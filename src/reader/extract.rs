//! Signal extraction over all the regions to scan, in parallel.
//!
//! Regions are processed in parallel on the current `rayon` thread pool, with
//! one reader per worker. A read overlapping two neighbouring regions is
//! processed once, by the first one. The result does not depend on the number
//! of threads: signals are sorted.

use std::io;

use rayon::prelude::*;

use super::alignment::{AlignmentError, AlignmentFile, RegionReader};
use super::filter::{FilterStats, ReadFilter};
use super::signals::{ClipOptions, Signal, SignalKind, read_signals};
use crate::regions::targets::{ScanRegion, TargetSet};

/// Number of median absolute deviations above the median insert size beyond
/// which a same-contig pair is a large insert.
const LARGE_INSERT_MADS: f64 = 6.0;

/// Scale of the median absolute deviation to the standard deviation of a
/// normal distribution.
const MAD_SCALE: f64 = 1.4826;

/// Signal extraction options.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ExtractOptions {
    /// Read filter.
    pub filter: ReadFilter,
    /// Soft-clip and poly(A/T) detection.
    pub clip: ClipOptions,
}

/// Template lengths above which proper pairs are only counted, not binned.
const MAX_BINNED_LENGTH: u64 = 10_000;

/// Template lengths of proper pairs, counted per length up to
/// [`MAX_BINNED_LENGTH`], plus an overflow count.
///
/// Memory is bounded whatever the number of pairs. The median and the MAD are
/// exact while the median is at most half the bound and the deviation rank of
/// the MAD falls below the overflow, which holds for any library of proper
/// pairs; otherwise overflowing lengths count as `MAX_BINNED_LENGTH + 1`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct InsertSizeHistogram {
    /// Pairs per template length (allocated on the first one).
    counts: Vec<u64>,
    /// Pairs above [`MAX_BINNED_LENGTH`].
    overflow: u64,
    /// All pairs.
    total: u64,
}

impl InsertSizeHistogram {
    /// Count one template length.
    pub fn add(&mut self, length: u64) {
        self.total += 1;
        if length > MAX_BINNED_LENGTH {
            self.overflow += 1;
            return;
        }
        #[allow(clippy::cast_possible_truncation)]
        if self.counts.is_empty() {
            self.counts = vec![0; MAX_BINNED_LENGTH as usize + 1];
        }
        #[allow(clippy::cast_possible_truncation)]
        let bin = length as usize;
        self.counts[bin] += 1;
    }

    /// Add the counts of another histogram.
    pub fn merge(&mut self, other: &Self) {
        if other.counts.is_empty() {
            self.overflow += other.overflow;
            self.total += other.total;
            return;
        }
        if self.counts.is_empty() {
            self.counts = vec![0; other.counts.len()];
        }
        for (count, other) in self.counts.iter_mut().zip(&other.counts) {
            *count += other;
        }
        self.overflow += other.overflow;
        self.total += other.total;
    }

    /// Number of pairs of a template length (overflow beyond the bound).
    fn count(&self, length: u64) -> u64 {
        if length > MAX_BINNED_LENGTH {
            return 0;
        }
        #[allow(clippy::cast_possible_truncation)]
        self.counts.get(length as usize).copied().unwrap_or(0)
    }

    /// Distribution of the counted lengths; `None` if empty.
    #[must_use]
    pub fn estimate(&self) -> Option<InsertSize> {
        if self.total == 0 {
            return None;
        }
        // Lower median, as `select_nth_unstable((n - 1) / 2)` on the values.
        let rank = (self.total - 1) / 2;
        let mut seen = 0;
        let mut median = MAX_BINNED_LENGTH + 1;
        for (length, &count) in (0..).zip(&self.counts) {
            seen += count;
            if seen > rank {
                median = length;
                break;
            }
        }
        // Deviations in increasing order: each distance d gathers the
        // lengths median - d and median + d; overflowing lengths come last.
        let mut seen = 0;
        let mut mad = (MAX_BINNED_LENGTH + 1).saturating_sub(median);
        for distance in 0..=median.max(MAX_BINNED_LENGTH.saturating_sub(median)) {
            seen += self.count(median + distance);
            if distance > 0 && distance <= median {
                seen += self.count(median - distance);
            }
            if seen > rank {
                mad = distance;
                break;
            }
        }
        #[allow(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            clippy::cast_precision_loss
        )]
        let spread = (LARGE_INSERT_MADS * MAD_SCALE * mad as f64).ceil() as u64;
        Some(InsertSize {
            pairs: self.total,
            median,
            mad,
            threshold: median + spread.max(1),
        })
    }
}

/// Insert size distribution of the proper pairs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct InsertSize {
    /// Number of proper pairs sampled (read 1 of each pair).
    pub pairs: u64,
    /// Median template length.
    pub median: u64,
    /// Median absolute deviation.
    pub mad: u64,
    /// Template length beyond which a same-contig pair is a large insert.
    pub threshold: u64,
}

/// Signals of all the regions, and statistics.
#[derive(Debug, Default)]
pub struct Extraction {
    /// Signals, sorted by contig, position, side and read.
    pub signals: Vec<Signal>,
    /// Read counts by filter outcome.
    pub stats: FilterStats,
    /// Soft-clips set aside for their low base quality.
    pub low_quality_clips: u64,
    /// Insert size distribution (`None` without proper pairs).
    pub insert_size: Option<InsertSize>,
}

impl Extraction {
    /// Number of signals of the given kind name (see [`SignalKind::name`]).
    #[must_use]
    pub fn count(&self, kind: &str) -> usize {
        self.signals
            .iter()
            .filter(|s| s.kind.name() == kind)
            .count()
    }

    /// Number of soft-clips with a poly(A/T) tail.
    #[must_use]
    pub fn tailed_clips(&self) -> usize {
        self.signals
            .iter()
            .filter(|s| matches!(s.kind, SignalKind::SoftClip { tail: Some(_), .. }))
            .count()
    }
}

/// Result of one region.
#[derive(Default)]
struct RegionResult {
    signals: Vec<Signal>,
    large_inserts: Vec<Signal>,
    insert_sizes: Vec<u64>,
    stats: FilterStats,
    low_quality_clips: u64,
}

/// Results of the regions processed by one worker, with the template lengths
/// counted in a histogram.
#[derive(Default)]
struct Accumulator {
    signals: Vec<Signal>,
    large_inserts: Vec<Signal>,
    insert_sizes: InsertSizeHistogram,
    stats: FilterStats,
    low_quality_clips: u64,
}

impl Accumulator {
    fn add(mut self, region: RegionResult) -> Self {
        self.signals.extend(region.signals);
        self.large_inserts.extend(region.large_inserts);
        for length in region.insert_sizes {
            self.insert_sizes.add(length);
        }
        self.stats.merge(&region.stats);
        self.low_quality_clips += region.low_quality_clips;
        self
    }

    fn merge(mut self, other: Self) -> Self {
        self.signals.extend(other.signals);
        self.large_inserts.extend(other.large_inserts);
        self.insert_sizes.merge(&other.insert_sizes);
        self.stats.merge(&other.stats);
        self.low_quality_clips += other.low_quality_clips;
        self
    }
}

/// Extract the signals of every region to scan.
///
/// # Errors
///
/// Fails if the alignment file cannot be read.
pub fn extract(
    file: &AlignmentFile,
    targets: &TargetSet,
    options: &ExtractOptions,
) -> Result<Extraction, AlignmentError> {
    let regions = targets.regions();
    let result = regions
        .par_iter()
        .enumerate()
        .map_init(
            || file.reader(),
            |reader, (i, region)| {
                let reader = reader
                    .as_mut()
                    .map_err(|e| io::Error::other(e.to_string()))?;
                // A read starting before the region and overlapping the
                // previous one belongs to the previous one.
                let previous_end = i
                    .checked_sub(1)
                    .map(|p| &regions[p])
                    .filter(|p| p.contig == region.contig)
                    .map(|p| p.end);
                scan_region(reader, file, region, previous_end, options)
            },
        )
        // The template lengths of each region go into one histogram per
        // split, not a list of all the pairs.
        .try_fold(Accumulator::default, |acc, region| {
            region.map(|region| acc.add(region))
        })
        .try_reduce(Accumulator::default, |a, b| Ok(a.merge(b)))
        .map_err(|source| AlignmentError::Io {
            path: file.path().to_path_buf(),
            source,
        })?;

    let mut extraction = Extraction {
        signals: result.signals,
        stats: result.stats,
        low_quality_clips: result.low_quality_clips,
        insert_size: result.insert_sizes.estimate(),
    };
    if let Some(insert_size) = extraction.insert_size {
        extraction
            .signals
            .extend(result.large_inserts.into_iter().filter(|s| {
                matches!(s.kind, SignalKind::LargeInsert { length } if length > insert_size.threshold)
            }));
    }
    extraction.signals.sort_unstable();
    extraction.signals.dedup();
    Ok(extraction)
}

/// Extract the signals of the reads of one region.
fn scan_region(
    reader: &mut RegionReader<'_>,
    file: &AlignmentFile,
    region: &ScanRegion,
    previous_end: Option<u64>,
    options: &ExtractOptions,
) -> io::Result<RegionResult> {
    let mut result = RegionResult::default();
    for record in reader.query(region)? {
        let record = record?;
        if let (Some(previous_end), Some(start)) =
            (previous_end, record.alignment_start().transpose()?)
            && (start.get() as u64 - 1) < region.start
            && (start.get() as u64 - 1) < previous_end
        {
            continue;
        }
        let exclusion = options.filter.check(record.as_ref())?;
        result.stats.add(exclusion);
        if exclusion.is_some() {
            continue;
        }
        let read = read_signals(
            record.as_ref(),
            file.header(),
            &options.clip,
            options.filter.min_mapq,
        )?;
        result.signals.extend(read.signals);
        result.large_inserts.extend(read.large_insert);
        result.insert_sizes.extend(read.insert_size);
        result.low_quality_clips += u64::from(read.low_quality_clips);
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reader::signals::{PolyTail, Side};
    use crate::regions::bed;
    use std::path::{Path, PathBuf};

    fn data(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data")
            .join(name)
    }

    fn run(bam: &str, threads: usize) -> Extraction {
        let reference = data("reference.fa");
        let file = AlignmentFile::open(&data(bam), Some(&reference)).unwrap();
        let probes = bed::read(&data("targets.bed")).unwrap();
        let targets = TargetSet::new(probes, 300, &file.contigs()).unwrap();
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| extract(&file, &targets, &ExtractOptions::default()))
            .unwrap()
    }

    /// Soft-clips by (position, side, tail).
    fn clips(extraction: &Extraction) -> Vec<(u64, Side, Option<u8>, usize)> {
        let mut groups: Vec<(u64, Side, Option<u8>, usize)> = Vec::new();
        for signal in &extraction.signals {
            if let SignalKind::SoftClip { tail, .. } = &signal.kind {
                let key = (signal.position, signal.side, tail.map(|t| t.base));
                match groups.last_mut() {
                    Some(g) if (g.0, g.1, g.2) == key => g.3 += 1,
                    _ => groups.push((key.0, key.1, key.2, 1)),
                }
            }
        }
        groups
    }

    #[test]
    fn positive_sample() {
        let extraction = run("positive.bam", 4);
        // Same counts as samtools (see tests/data/README.md for the truth).
        assert_eq!(extraction.stats.kept, 2005);
        // The other 2 duplicates are mates placed on chr2, outside the regions.
        assert_eq!(extraction.stats.duplicate, 4);
        assert_eq!(
            clips(&extraction),
            [
                (3000, Side::Left, None, 10),
                (10000, Side::Left, Some(b'A'), 21),
                (10012, Side::Right, None, 47),
            ]
        );
        // Every poly(A) clip holds at least 20 A from the junction.
        for signal in &extraction.signals {
            if let SignalKind::SoftClip {
                tail: Some(PolyTail { length, .. }),
                ..
            } = signal.kind
            {
                assert!(length >= 20, "{signal:?}");
            }
        }
        assert_eq!(extraction.count("mate_unmapped"), 14);
        assert_eq!(extraction.count("mate_low_mapq"), 123);
        assert_eq!(extraction.count("mate_elsewhere"), 0);
        assert_eq!(extraction.count("large_insert"), 0);
        // Mate signals point to the insertion: reads of the left flank end
        // before the right junction (10012), reads of the right flank start
        // after the left one (10000); the flanks overlap on the TSD.
        for signal in extraction
            .signals
            .iter()
            .filter(|s| s.kind.is_mate_signal())
        {
            assert_eq!(signal.contig, 0);
            match signal.side {
                Side::Right => assert!((9_550..=10_012).contains(&signal.position), "{signal:?}"),
                Side::Left => assert!((10_000..10_450).contains(&signal.position), "{signal:?}"),
            }
        }
        let insert_size = extraction.insert_size.unwrap();
        assert!((280..=320).contains(&insert_size.median), "{insert_size:?}");
    }

    #[test]
    fn negative_sample() {
        let extraction = run("negative.bam", 4);
        assert_eq!(clips(&extraction), [(3000, Side::Left, None, 10)]);
        assert!(!extraction.signals.iter().any(|s| s.kind.is_mate_signal()));
    }

    #[test]
    fn same_result_with_any_thread_count_and_format() {
        let reference = run("positive.bam", 1);
        for (bam, threads) in [
            ("positive.bam", 8),
            ("positive.cram", 1),
            ("positive.cram", 8),
        ] {
            let other = run(bam, threads);
            assert_eq!(other.signals, reference.signals, "{bam}, {threads} threads");
            assert_eq!(other.stats, reference.stats, "{bam}, {threads} threads");
            assert_eq!(other.insert_size, reference.insert_size);
        }
    }

    #[test]
    fn reads_overlapping_two_regions_are_processed_once() {
        let file = AlignmentFile::open(&data("positive.bam"), None).unwrap();
        let targets = |probes: &[(u64, u64)]| {
            let probes = probes
                .iter()
                .map(|&(start, end)| bed::Probe {
                    contig: "chr1".to_owned(),
                    start,
                    end,
                    name: None,
                })
                .collect();
            TargetSet::new(probes, 0, &file.contigs()).unwrap()
        };
        // Two regions 1 bp apart: most reads overlap both.
        let split = targets(&[(9_850, 10_005), (10_006, 10_150)]);
        assert_eq!(split.regions().len(), 2);
        let whole = targets(&[(9_850, 10_150)]);
        let options = ExtractOptions::default();
        let a = extract(&file, &split, &options).unwrap();
        let b = extract(&file, &whole, &options).unwrap();
        assert_eq!(a.stats, b.stats);
        assert_eq!(a.signals, b.signals);
    }

    fn histogram(lengths: &[u64]) -> InsertSizeHistogram {
        let mut histogram = InsertSizeHistogram::default();
        for &length in lengths {
            histogram.add(length);
        }
        histogram
    }

    /// Median and MAD of the values themselves (the estimate before #23).
    fn exact(mut lengths: Vec<u64>) -> (u64, u64) {
        let mid = (lengths.len() - 1) / 2;
        let median = *lengths.select_nth_unstable(mid).1;
        let mut deviations: Vec<u64> = lengths.iter().map(|&l| l.abs_diff(median)).collect();
        (median, *deviations.select_nth_unstable(mid).1)
    }

    #[test]
    fn insert_size_estimate() {
        assert_eq!(histogram(&[]).estimate(), None);
        let estimate = histogram(&[290, 300, 310, 300, 1_000]).estimate().unwrap();
        assert_eq!(estimate.pairs, 5);
        assert_eq!(estimate.median, 300);
        assert_eq!(estimate.mad, 10);
        assert_eq!(estimate.threshold, 300 + 89);
        // Lengths beyond the bound only count.
        let estimate = histogram(&[300, 300, 20_000]).estimate().unwrap();
        assert_eq!((estimate.pairs, estimate.median, estimate.mad), (3, 300, 0));
    }

    #[test]
    fn histogram_matches_the_exact_estimate() {
        // Deterministic pseudo-random samples (xorshift), with outliers and
        // lengths beyond the bound.
        let mut state = 0x2545_f491_4f6c_dd1d_u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for n in [1, 2, 3, 10, 101, 1_000, 10_000] {
            let lengths: Vec<u64> = (0..n)
                .map(|_| match next() % 20 {
                    0 => next() % 50_000,
                    1 => 1 + next() % 100,
                    _ => 200 + next() % 200,
                })
                .collect();
            // Split in two histograms and merge, as the workers do.
            let (a, b) = lengths.split_at(n / 3);
            let mut merged = histogram(a);
            merged.merge(&histogram(b));
            assert_eq!(merged, histogram(&lengths));
            let estimate = merged.estimate().unwrap();
            assert_eq!(
                (estimate.median, estimate.mad),
                exact(lengths.clone()),
                "{n} lengths"
            );
        }
    }
}
