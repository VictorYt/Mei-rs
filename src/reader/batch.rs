//! Batched region queries, contig by contig, in parallel.
//!
//! A BAM query is cheap: each interval is queried on its own. A CRAM query
//! decodes whole containers (about 80 kb of an exome each), which would be
//! decoded once per interval they overlap: neighbouring intervals are queried
//! together, and the caller keeps the reads overlapping its intervals (see
//! [`overlap`]). Contigs are processed one at a time: the reference cache of
//! a CRAM file is loaded once per contig and holds one contig.

use std::io;
use std::ops::Range;

use noodles::sam;
use rayon::prelude::*;

use super::alignment::{AlignmentError, AlignmentFile, Format, RegionReader};
use crate::regions::targets::ScanRegion;

/// Largest gap between two intervals of a CRAM batch, in bp: larger than the
/// span of a CRAM container on an exome.
pub const CRAM_BATCH_GAP: u64 = 100_000;

/// Largest span of a CRAM batch, in bp, to keep enough batches per contig for
/// the workers.
pub const CRAM_BATCH_SPAN: u64 = 2_000_000;

/// Splits per worker of the parallel iterator over the batches of a contig.
/// `rayon` opens a reader (header parsing) per split: without a minimum
/// split length, it splits down to a few batches.
const SPLITS_PER_THREAD: usize = 8;

/// Consecutive intervals of one contig, read with one query.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Batch {
    /// Indices of the intervals.
    pub range: Range<usize>,
    /// Region to query: from the start of the first interval to the end of
    /// the last one.
    pub span: ScanRegion,
}

/// Group intervals, sorted by contig then start (0-based half-open, possibly
/// overlapping), into batches: one interval each for BAM; for CRAM,
/// neighbouring intervals of a contig, at most [`CRAM_BATCH_GAP`] apart and
/// spanning at most [`CRAM_BATCH_SPAN`].
#[must_use]
pub fn batches(intervals: &[ScanRegion], format: Format) -> Vec<Batch> {
    let mut batches: Vec<Batch> = Vec::new();
    for (i, interval) in intervals.iter().enumerate() {
        if format == Format::Cram
            && let Some(batch) = batches.last_mut()
            && batch.span.contig == interval.contig
            && interval.start.saturating_sub(batch.span.end) <= CRAM_BATCH_GAP
            && interval.end.max(batch.span.end) - batch.span.start <= CRAM_BATCH_SPAN
        {
            batch.range.end = i + 1;
            batch.span.end = batch.span.end.max(interval.end);
            continue;
        }
        batches.push(Batch {
            range: i..i + 1,
            span: interval.clone(),
        });
    }
    batches
}

/// Indices of the intervals a read overlaps, among intervals sorted by start
/// and by end, with the rule of the index queries: the 1-based inclusive
/// alignment interval of the read, a single base without alignment span.
///
/// # Errors
///
/// Fails if the position or the CIGAR of the read cannot be decoded.
pub fn overlap(
    intervals: &[ScanRegion],
    record: &dyn sam::alignment::Record,
) -> io::Result<Range<usize>> {
    let (Some(start), Some(end)) = (
        record.alignment_start().transpose()?,
        record.alignment_end().transpose()?,
    ) else {
        return Ok(0..0);
    };
    // 0-based half-open interval of the read.
    let (start, end) = (start.get() as u64 - 1, end.get() as u64);
    let first = intervals.partition_point(|r| r.end <= start);
    let last = first + intervals[first..].partition_point(|r| r.start < end);
    Ok(first..last)
}

/// Process the batches in parallel, contig by contig, with one reader per
/// split: `process` reads a batch, `add` folds its result into the
/// accumulator of the split, and `merge` combines the accumulators.
///
/// # Errors
///
/// Fails if the alignment file or the reference cannot be read.
pub fn fold<R, A>(
    file: &AlignmentFile,
    batches: &[Batch],
    process: impl Fn(&mut RegionReader<'_>, &Batch) -> io::Result<R> + Sync,
    add: impl Fn(A, R) -> A + Sync,
    merge: impl Fn(A, A) -> A + Sync,
) -> Result<A, AlignmentError>
where
    R: Send,
    A: Default + Send,
{
    let error = |source| AlignmentError::Io {
        path: file.path().to_path_buf(),
        source,
    };
    let mut result = A::default();
    let mut first = 0;
    while first < batches.len() {
        let contig = &batches[first].span.contig;
        let last = first
            + batches[first..]
                .iter()
                .take_while(|b| &b.span.contig == contig)
                .count();
        file.load_reference(contig).map_err(error)?;
        let min_len = (last - first) / (rayon::current_num_threads() * SPLITS_PER_THREAD);
        let contig_result = batches[first..last]
            .par_iter()
            .with_min_len(min_len.max(1))
            .map_init(
                || file.reader(),
                |reader, batch| {
                    let reader = reader
                        .as_mut()
                        .map_err(|e| io::Error::other(e.to_string()))?;
                    process(reader, batch)
                },
            )
            .try_fold(A::default, |acc, r| r.map(|r| add(acc, r)))
            .try_reduce(A::default, |a, b| Ok(merge(a, b)))
            .map_err(error)?;
        result = merge(result, contig_result);
        file.clear_reference_cache();
        first = last;
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn interval(contig: &str, start: u64, end: u64) -> ScanRegion {
        ScanRegion {
            contig: contig.to_owned(),
            start,
            end,
        }
    }

    #[test]
    fn batches_by_format() {
        let intervals = [
            interval("chr1", 0, 100),
            interval("chr1", 50_000, 50_100),
            interval("chr1", 200_000, 200_100),
            interval("chr2", 0, 100),
            interval("chr2", 1, 3),
        ];
        let ranges = |format| -> Vec<Range<usize>> {
            batches(&intervals, format)
                .into_iter()
                .map(|b| b.range)
                .collect()
        };
        assert_eq!(ranges(Format::Bam), [0..1, 1..2, 2..3, 3..4, 4..5]);
        // Gaps up to 100 kb, not across contigs; overlapping intervals.
        assert_eq!(ranges(Format::Cram), [0..2, 2..3, 3..5]);
        assert_eq!(
            batches(&intervals, Format::Cram)[0].span,
            interval("chr1", 0, 50_100)
        );

        // Span limit.
        let far: Vec<_> = (0..30)
            .map(|i| interval("chr1", i * 90_000, i * 90_000 + 100))
            .collect();
        let spans: Vec<_> = batches(&far, Format::Cram)
            .into_iter()
            .map(|b| b.span.end - b.span.start)
            .collect();
        assert!(spans.iter().all(|&s| s <= CRAM_BATCH_SPAN), "{spans:?}");
        assert_eq!(spans.len(), 2);
    }
}
