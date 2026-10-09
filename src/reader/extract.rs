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

/// Insert size distribution of the proper pairs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct InsertSize {
    /// Number of proper pairs sampled (read 1 of each pair).
    pub pairs: usize,
    /// Median template length.
    pub median: u64,
    /// Median absolute deviation.
    pub mad: u64,
    /// Template length beyond which a same-contig pair is a large insert.
    pub threshold: u64,
}

impl InsertSize {
    /// Distribution of a sample of template lengths; `None` if empty.
    #[must_use]
    pub fn estimate(mut lengths: Vec<u64>) -> Option<Self> {
        let median = median(&mut lengths)?;
        let mut deviations: Vec<u64> = lengths.iter().map(|&l| l.abs_diff(median)).collect();
        let mad = median_of(&mut deviations).unwrap_or(0);
        #[allow(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            clippy::cast_precision_loss
        )]
        let spread = (LARGE_INSERT_MADS * MAD_SCALE * mad as f64).ceil() as u64;
        Some(Self {
            pairs: lengths.len(),
            median,
            mad,
            threshold: median + spread.max(1),
        })
    }
}

fn median(values: &mut [u64]) -> Option<u64> {
    median_of(values)
}

fn median_of(values: &mut [u64]) -> Option<u64> {
    if values.is_empty() {
        return None;
    }
    let mid = (values.len() - 1) / 2;
    Some(*values.select_nth_unstable(mid).1)
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
    let results = regions
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
        .collect::<io::Result<Vec<RegionResult>>>()
        .map_err(|source| AlignmentError::Io {
            path: file.path().to_path_buf(),
            source,
        })?;

    let mut extraction = Extraction::default();
    let mut large_inserts = Vec::new();
    let mut insert_sizes = Vec::new();
    for result in results {
        extraction.signals.extend(result.signals);
        large_inserts.extend(result.large_inserts);
        insert_sizes.extend(result.insert_sizes);
        extraction.stats.merge(&result.stats);
        extraction.low_quality_clips += result.low_quality_clips;
    }
    extraction.insert_size = InsertSize::estimate(insert_sizes);
    if let Some(insert_size) = extraction.insert_size {
        extraction.signals.extend(large_inserts.into_iter().filter(|s| {
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

    #[test]
    fn insert_size_estimate() {
        assert_eq!(InsertSize::estimate(vec![]), None);
        let estimate = InsertSize::estimate(vec![290, 300, 310, 300, 1_000]).unwrap();
        assert_eq!(estimate.pairs, 5);
        assert_eq!(estimate.median, 300);
        assert_eq!(estimate.mad, 10);
        assert_eq!(estimate.threshold, 300 + 89);
    }
}
