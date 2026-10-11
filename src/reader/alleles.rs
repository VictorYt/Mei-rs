//! Reads of each allele at the junctions of the candidates, for the
//! genotype (see [`crate::genotype`]).
//!
//! At a candidate site (its left and right junctions, or one of them), among
//! the reads that pass the read filter:
//!
//! - an **alternative** read is soft-clipped at a junction (its alignment
//!   ends at the left junction, or starts at the right one, within
//!   [`JUNCTION_SLACK`]) by at least `overhang` bases of median quality at
//!   least `min_quality`: it holds the inserted sequence;
//! - a **reference** read is aligned over the whole site (both junctions)
//!   with at least `overhang` bases on each side: it holds the site without
//!   the insertion. Over the whole site, so that an alternative read of one
//!   junction, which goes on over the TSD, is never a reference read of the
//!   other.
//!
//! The same overhang on both alleles keeps them comparable: a read that
//! crosses a junction by fewer bases counts for neither.
//!
//! Coordinates: junctions are boundaries (see [`crate::reader::signals`]).

use std::io;

use noodles::sam::alignment::Record;
use noodles::sam::alignment::record::cigar::op::Kind;

use super::alignment::{AlignmentError, AlignmentFile};
use super::batch;
use super::filter::ReadFilter;
use super::junctions::JUNCTION_SLACK;
use super::signals::clip_lengths;
use crate::regions::targets::ScanRegion;

/// A candidate site: its junctions (at least one).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Site {
    /// Contig name.
    pub contig: String,
    /// End of the upstream flank: alternative reads are clipped on their
    /// right end there.
    pub left: Option<u64>,
    /// Start of the downstream flank: alternative reads are clipped on their
    /// left end there.
    pub right: Option<u64>,
}

impl Site {
    /// First and last junction.
    fn span(&self) -> Option<(u64, u64)> {
        let junctions = [self.left, self.right];
        let lo = junctions.iter().flatten().min()?;
        let hi = junctions.iter().flatten().max()?;
        Some((*lo, *hi))
    }
}

/// Reads of each allele at a site.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SiteCounts {
    /// Alternative reads at the left junction.
    pub left: u32,
    /// Alternative reads at the right junction.
    pub right: u32,
    /// Reference reads over the whole site.
    pub reference: u32,
}

/// Count the reads of each allele at each site: one result per site, in
/// order (zero for a site without junction).
///
/// # Errors
///
/// Fails if the alignment file cannot be read.
pub fn allele_counts(
    file: &AlignmentFile,
    sites: &[Site],
    filter: ReadFilter,
    overhang: u64,
    min_quality: u8,
) -> Result<Vec<SiteCounts>, AlignmentError> {
    let mut order: Vec<(usize, u64, u64)> = sites
        .iter()
        .enumerate()
        .filter_map(|(i, site)| site.span().map(|(lo, hi)| (i, lo, hi)))
        .collect();
    order.sort_by(|a, b| (&sites[a.0].contig, a.1).cmp(&(&sites[b.0].contig, b.1)));
    // Reads overlapping the site and the base on each side.
    let windows: Vec<ScanRegion> = order
        .iter()
        .map(|&(i, lo, hi)| ScanRegion {
            contig: sites[i].contig.clone(),
            start: lo.saturating_sub(1),
            end: hi + 1,
        })
        .collect();
    let counts = batch::fold(
        file,
        &batch::batches(&windows, file.format()),
        |reader, batch| {
            let mut counts = vec![SiteCounts::default(); batch.range.len()];
            for record in reader.query(&batch.span)? {
                let record = record?;
                let covered = batch::overlap(&windows[batch.range.clone()], record.as_ref())?;
                if covered.is_empty() || filter.check(record.as_ref())?.is_some() {
                    continue;
                }
                let Some(read) = Read::new(record.as_ref(), overhang, min_quality)? else {
                    continue;
                };
                for k in covered {
                    let (i, lo, hi) = order[batch.range.start + k];
                    read.count(&sites[i], (lo, hi), overhang, &mut counts[k]);
                }
            }
            Ok((batch.range.start, counts))
        },
        |mut acc: Vec<(usize, Vec<SiteCounts>)>, batch| {
            acc.push(batch);
            acc
        },
        |mut a, b| {
            a.extend(b);
            a
        },
    )?;
    let mut result = vec![SiteCounts::default(); sites.len()];
    for (first, counts) in counts {
        for (offset, count) in counts.into_iter().enumerate() {
            result[order[first + offset].0] = count;
        }
    }
    Ok(result)
}

/// Alignment of a read, as the counts see it.
struct Read {
    /// 0-based alignment start.
    start: u64,
    /// 0-based alignment end (exclusive).
    end: u64,
    /// Soft-clip on the left end, long enough and of good enough quality.
    left_clip: bool,
    /// Same on the right end.
    right_clip: bool,
}

impl Read {
    fn new(record: &dyn Record, overhang: u64, min_quality: u8) -> io::Result<Option<Self>> {
        let (Some(start), Some(end)) = (
            record.alignment_start().transpose()?,
            record.alignment_end().transpose()?,
        ) else {
            return Ok(None);
        };
        let ops = record
            .cigar()
            .iter()
            .map(|op| op.map(|op| (op.kind(), op.len())))
            .collect::<io::Result<Vec<(Kind, usize)>>>()?;
        let (left, right) = clip_lengths(&ops);
        let qualities = record
            .quality_scores()
            .iter()
            .collect::<io::Result<Vec<u8>>>()?;
        let n = qualities.len();
        let good = |length: usize, range: std::ops::Range<usize>| {
            length as u64 >= overhang
                && qualities
                    .get(range)
                    .is_none_or(|q| median(q) >= min_quality)
        };
        Ok(Some(Self {
            start: start.get() as u64 - 1,
            end: end.get() as u64,
            left_clip: good(left, 0..left),
            right_clip: good(right, n.saturating_sub(right)..n),
        }))
    }

    /// Count the read at a site spanning `lo` to `hi`.
    fn count(&self, site: &Site, (lo, hi): (u64, u64), overhang: u64, counts: &mut SiteCounts) {
        let near = |junction: Option<u64>, position: u64| {
            junction.is_some_and(|j| j.abs_diff(position) <= JUNCTION_SLACK)
        };
        let left = self.right_clip && near(site.left, self.end);
        let right = self.left_clip && near(site.right, self.start);
        counts.left += u32::from(left);
        counts.right += u32::from(right);
        if !left && !right && self.start + overhang <= lo && self.end >= hi + overhang {
            counts.reference += 1;
        }
    }
}

/// Median of the qualities (lower median); 255 if empty.
fn median(qualities: &[u8]) -> u8 {
    if qualities.is_empty() {
        return u8::MAX;
    }
    let mut sorted = qualities.to_vec();
    sorted.sort_unstable();
    sorted[(sorted.len() - 1) / 2]
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn heterozygous_alu() {
        let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
        let file = AlignmentFile::open(&data.join("positive.bam"), None).unwrap();
        let site = Site {
            contig: "chr1".to_owned(),
            left: Some(10_012),
            right: Some(10_000),
        };
        let c = allele_counts(&file, &[site], ReadFilter::default(), 15, 20)
            .unwrap()
            .remove(0);
        // Heterozygous: about as many reads of each allele at each junction.
        for alternative in [c.left, c.right] {
            let fraction = f64::from(alternative) / f64::from(alternative + c.reference);
            assert!((0.3..0.7).contains(&fraction), "{c:?}");
        }
    }

    #[test]
    fn alternative_reads_of_one_junction_are_not_reference_reads() {
        // A TSD of 20 bp: the reads clipped at the right junction go on
        // over the left one.
        let site = Site {
            contig: "chr1".to_owned(),
            left: Some(1_020),
            right: Some(1_000),
        };
        let span = site.span().unwrap();
        let mut counts = SiteCounts::default();
        let clipped_right = Read {
            start: 1_000,
            end: 1_100,
            left_clip: true,
            right_clip: false,
        };
        clipped_right.count(&site, span, 15, &mut counts);
        let clipped_left = Read {
            start: 920,
            end: 1_021,
            left_clip: false,
            right_clip: true,
        };
        clipped_left.count(&site, span, 15, &mut counts);
        // Spanning the left junction only: neither.
        let partial = Read {
            start: 990,
            end: 1_090,
            left_clip: false,
            right_clip: false,
        };
        partial.count(&site, span, 15, &mut counts);
        let spanning = Read {
            start: 950,
            end: 1_050,
            left_clip: false,
            right_clip: false,
        };
        spanning.count(&site, span, 15, &mut counts);
        assert_eq!(
            counts,
            SiteCounts {
                left: 1,
                right: 1,
                reference: 1
            }
        );
    }
}
