//! Breakpoints of the candidates without junction soft-clips, looked for
//! again among all the soft-clips of their reads (those with MEI mates, see
//! `scan`).
//!
//! Without soft-clips of at least `--min-short-clip-len` bases at a junction,
//! a candidate is placed from its mate signals, within an interval of a few
//! hundred bases ([`crate::cluster::Candidate::mate_interval`]). The reads
//! that cross the junction are often still there: aligners clip them by a
//! few bases only, or match the first inserted bases with mismatches. The
//! interval is read again, and the boundary where most reads are clipped by
//! at least [`MIN_CLIP`] bases, at least [`MIN_CLIPS`] of them, is the
//! breakpoint (on ties, the one nearest the estimate).
//!
//! Reads pass the read filter of the extraction.
//!
//! Coordinates: boundaries (see [`crate::reader::signals`]); windows 0-based
//! half-open.

use std::collections::HashMap;

use noodles::sam::alignment::record::cigar::op::Kind;

use super::alignment::{AlignmentError, AlignmentFile};
use super::batch;
use super::filter::ReadFilter;
use super::signals::clip_lengths;
use crate::regions::targets::ScanRegion;

/// Minimum length of a soft-clip, in bp, to mark a junction.
pub const MIN_CLIP: usize = 3;

/// Minimum number of reads clipped at the same boundary.
pub const MIN_CLIPS: usize = 2;

/// A breakpoint to look for.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BreakpointRequest {
    /// Boundaries to look at: the first and the last, inclusive.
    pub first: u64,
    /// Last boundary looked at.
    pub last: u64,
    /// Contig name.
    pub contig: String,
    /// Estimated breakpoint, to break ties.
    pub estimate: u64,
}

impl BreakpointRequest {
    /// Region to read: every read clipped at one of the boundaries overlaps
    /// it.
    fn window(&self) -> ScanRegion {
        ScanRegion {
            contig: self.contig.clone(),
            start: self.first.saturating_sub(1),
            end: self.last + 1,
        }
    }
}

/// Breakpoint found: the boundary and the number of reads clipped there.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Breakpoint {
    /// Boundary where the reads are clipped.
    pub position: u64,
    /// Reads clipped there by at least [`MIN_CLIP`] bases.
    pub clips: usize,
}

/// Look for the breakpoint of each request (sorted by contig and first
/// boundary, as the batches require); one result per request, in order.
///
/// # Errors
///
/// Fails if the alignment file cannot be read.
pub fn find_breakpoints(
    file: &AlignmentFile,
    requests: &[BreakpointRequest],
    filter: &ReadFilter,
) -> Result<Vec<Option<Breakpoint>>, AlignmentError> {
    let windows: Vec<ScanRegion> = requests.iter().map(BreakpointRequest::window).collect();
    let longest = windows.iter().map(ScanRegion::len).max().unwrap_or(0);
    let found = batch::fold(
        file,
        &batch::batches(&windows, file.format()),
        |reader, batch| {
            let range = batch.range.clone();
            let mut counts: Vec<HashMap<u64, usize>> = vec![HashMap::new(); range.len()];
            for record in reader.query(&batch.span)? {
                let record = record?;
                let flags = record.flags()?;
                if flags.is_unmapped() || filter.check(&record)?.is_some() {
                    continue;
                }
                let ops = record
                    .cigar()
                    .iter()
                    .map(|op| op.map(|op| (op.kind(), op.len())))
                    .collect::<std::io::Result<Vec<(Kind, usize)>>>()?;
                let (left_clip, right_clip) = clip_lengths(&ops);
                if left_clip < MIN_CLIP && right_clip < MIN_CLIP {
                    continue;
                }
                let (Some(start), Some(end)) = (
                    record.alignment_start().transpose()?,
                    record.alignment_end().transpose()?,
                ) else {
                    continue;
                };
                let (start, end) = (start.get() as u64 - 1, end.get() as u64);
                let boundaries = [(left_clip, start), (right_clip, end)];
                // Requests sorted by first boundary, at most `longest` long.
                let requests = &requests[range.clone()];
                let after = requests.partition_point(|r| r.first <= end);
                for (k, request) in requests[..after].iter().enumerate().rev() {
                    if request.first + longest < start {
                        break;
                    }
                    for &(clip, boundary) in &boundaries {
                        if clip >= MIN_CLIP && (request.first..=request.last).contains(&boundary) {
                            *counts[k].entry(boundary).or_default() += 1;
                        }
                    }
                }
            }
            Ok(range
                .zip(counts)
                .map(|(i, counts)| (i, best(&counts, requests[i].estimate)))
                .collect::<Vec<_>>())
        },
        |mut acc: Vec<(usize, Option<Breakpoint>)>, batch| {
            acc.extend(batch);
            acc
        },
        |mut a, b| {
            a.extend(b);
            a
        },
    )?;
    let mut breakpoints = vec![None; requests.len()];
    for (i, breakpoint) in found {
        breakpoints[i] = breakpoint;
    }
    Ok(breakpoints)
}

/// The boundary with the most clipped reads, at least [`MIN_CLIPS`]; on
/// ties, the nearest to the estimate, then the first.
fn best(counts: &HashMap<u64, usize>, estimate: u64) -> Option<Breakpoint> {
    counts
        .iter()
        .filter(|&(_, &n)| n >= MIN_CLIPS)
        .min_by_key(|&(&position, &n)| {
            (std::cmp::Reverse(n), position.abs_diff(estimate), position)
        })
        .map(|(&position, &clips)| Breakpoint { position, clips })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::{Path, PathBuf};

    fn data(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data")
            .join(name)
    }

    #[test]
    fn most_clipped_boundary() {
        let counts: HashMap<u64, usize> = [(100, 2), (140, 3), (180, 3), (200, 1)].into();
        assert_eq!(
            best(&counts, 170),
            Some(Breakpoint {
                position: 180,
                clips: 3
            })
        );
        assert_eq!(
            best(&counts, 0),
            Some(Breakpoint {
                position: 140,
                clips: 3
            })
        );
        let one: HashMap<u64, usize> = [(100, 1)].into();
        assert_eq!(best(&one, 100), None);
    }

    #[test]
    fn breakpoint_of_the_alu() {
        // AluY of `positive`: right junction 10000 (reads clipped on their
        // left end), left junction 10012; more reads clipped at the first.
        let request = BreakpointRequest {
            contig: "chr1".to_owned(),
            first: 9_800,
            last: 10_100,
            estimate: 10_050,
        };
        let outside = BreakpointRequest {
            contig: "chr1".to_owned(),
            first: 12_000,
            last: 12_300,
            estimate: 12_000,
        };
        for (file, reference) in [
            ("positive.bam", None),
            ("positive.cram", Some(data("reference.fa"))),
        ] {
            let file = AlignmentFile::open(&data(file), reference.as_deref()).unwrap();
            let found = find_breakpoints(
                &file,
                &[request.clone(), outside.clone()],
                &ReadFilter::default(),
            )
            .unwrap();
            let alu = found[0].unwrap();
            assert!([10_000, 10_012].contains(&alu.position), "{alu:?}");
            assert!(alu.clips >= 20, "{alu:?}");
            assert_eq!(found[1], None);
        }
    }
}
