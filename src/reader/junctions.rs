//! Reads of the junctions of the candidates, read again for the local
//! assembly (see [`crate::assembly`]).
//!
//! The signals keep only the clipped part of a read; the assembly needs the
//! whole reads that cross a junction. For each candidate, one window around
//! its signals is read again, and three sets of reads are kept, each at most
//! [`MAX_READS`]:
//!
//! - **left**: reads soft-clipped on their right end, whose alignment ends
//!   within [`JUNCTION_SLACK`] of the left junction (end of the upstream
//!   flank);
//! - **right**: reads soft-clipped on their left end, whose alignment starts
//!   within [`JUNCTION_SLACK`] of the right junction (start of the downstream
//!   flank);
//! - **inside**: unmapped reads placed in the window (at their anchored
//!   mate): their sequence lies in the inserted element, in either
//!   orientation.
//!
//! Anchored reads pass the read filter of the extraction. Bases below
//! [`MIN_BASE_QUALITY`] are masked (`N`).
//!
//! Coordinates: junctions are boundaries, windows 0-based half-open (see
//! [`crate::reader::signals`]).

use std::io;

use noodles::sam::alignment::Record;
use noodles::sam::alignment::record::cigar::op::Kind;

use super::alignment::{AlignmentError, AlignmentFile};
use super::batch;
use super::filter::ReadFilter;
use super::signals::clip_lengths;
use crate::regions::targets::ScanRegion;

/// Maximum number of reads kept per junction (and inside the element).
pub const MAX_READS: usize = 300;

/// Distance, in bp, between the end of an alignment and a junction for a
/// clipped read to cross it: aligners place a clip a few bases away from
/// the junction when the inserted bases match the reference.
pub const JUNCTION_SLACK: u64 = 5;

/// Bases below this quality are masked.
pub const MIN_BASE_QUALITY: u8 = 10;

/// Reads to gather for one candidate.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct JunctionRequest {
    /// Window to read: the junctions and the anchors of the mates.
    pub window: ScanRegion,
    /// End of the upstream flank, if any.
    pub left: Option<u64>,
    /// Start of the downstream flank, if any.
    pub right: Option<u64>,
}

/// Reads gathered for one candidate (upper case, masked).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct JunctionReads {
    /// Reads crossing the left junction.
    pub left: Vec<Vec<u8>>,
    /// Reads crossing the right junction.
    pub right: Vec<Vec<u8>>,
    /// Unmapped reads placed in the window.
    pub inside: Vec<Vec<u8>>,
}

/// Gather the reads of each request (sorted by contig and window start, as
/// the batches require); one result per request, in order.
///
/// # Errors
///
/// Fails if the alignment file cannot be read.
pub fn fetch_junction_reads(
    file: &AlignmentFile,
    requests: &[JunctionRequest],
    filter: &ReadFilter,
) -> Result<Vec<JunctionReads>, AlignmentError> {
    let windows: Vec<ScanRegion> = requests.iter().map(|r| r.window.clone()).collect();
    let found = batch::fold(
        file,
        &batch::batches(&windows, file.format()),
        |reader, batch| {
            let mut found: Vec<(usize, JunctionReads)> = batch
                .range
                .clone()
                .map(|w| (w, JunctionReads::default()))
                .collect();
            for record in reader.query(&batch.span)? {
                let record = record?;
                let Some(read) = classify(record.as_ref(), *filter)? else {
                    continue;
                };
                for (w, reads) in &mut found {
                    add(&requests[*w], &read, reads);
                }
            }
            Ok(found)
        },
        |mut acc: Vec<(usize, JunctionReads)>, batch| {
            acc.extend(batch);
            acc
        },
        |mut a, b| {
            a.extend(b);
            a
        },
    )?;
    let mut reads = vec![JunctionReads::default(); requests.len()];
    for (w, r) in found {
        reads[w] = r;
    }
    Ok(reads)
}

/// A read of a window, as the assembly sees it.
#[derive(Debug)]
enum WindowRead {
    /// Anchored read: 0-based alignment start and end (exclusive), soft-clip
    /// lengths on the left and right, masked sequence.
    Anchored {
        start: u64,
        end: u64,
        left_clip: usize,
        right_clip: usize,
        sequence: Vec<u8>,
    },
    /// Unmapped read placed at its mate (0-based position).
    Unmapped { position: u64, sequence: Vec<u8> },
}

/// A soft-clipped anchored read or a placed unmapped read; `None` for the
/// others.
fn classify(record: &dyn Record, filter: ReadFilter) -> io::Result<Option<WindowRead>> {
    let flags = record.flags()?;
    if flags.is_unmapped() {
        if flags.is_mate_unmapped() || flags.is_secondary() || flags.is_supplementary() {
            return Ok(None);
        }
        let Some(position) = record.alignment_start().transpose()? else {
            return Ok(None);
        };
        return Ok(Some(WindowRead::Unmapped {
            position: position.get() as u64 - 1,
            sequence: masked(record)?,
        }));
    }
    if filter.check(record)?.is_some() {
        return Ok(None);
    }
    let ops = record
        .cigar()
        .iter()
        .map(|op| op.map(|op| (op.kind(), op.len())))
        .collect::<io::Result<Vec<(Kind, usize)>>>()?;
    let (left_clip, right_clip) = clip_lengths(&ops);
    if left_clip == 0 && right_clip == 0 {
        return Ok(None);
    }
    let (Some(start), Some(end)) = (
        record.alignment_start().transpose()?,
        record.alignment_end().transpose()?,
    ) else {
        return Ok(None);
    };
    Ok(Some(WindowRead::Anchored {
        start: start.get() as u64 - 1,
        end: end.get() as u64,
        left_clip,
        right_clip,
        sequence: masked(record)?,
    }))
}

/// Add a read to the reads of a request, if it belongs to them.
fn add(request: &JunctionRequest, read: &WindowRead, reads: &mut JunctionReads) {
    let near = |junction: Option<u64>, position: u64| {
        junction.is_some_and(|j| j.abs_diff(position) <= JUNCTION_SLACK)
    };
    let push = |set: &mut Vec<Vec<u8>>, sequence: &[u8]| {
        if set.len() < MAX_READS {
            set.push(sequence.to_vec());
        }
    };
    match read {
        WindowRead::Anchored {
            start,
            end,
            left_clip,
            right_clip,
            sequence,
        } => {
            if *right_clip > 0 && near(request.left, *end) {
                push(&mut reads.left, sequence);
            }
            if *left_clip > 0 && near(request.right, *start) {
                push(&mut reads.right, sequence);
            }
        }
        WindowRead::Unmapped { position, sequence } => {
            if (request.window.start..request.window.end).contains(position) {
                push(&mut reads.inside, sequence);
            }
        }
    }
}

/// Sequence of a read, upper case, with the bases below
/// [`MIN_BASE_QUALITY`] masked.
fn masked(record: &dyn Record) -> io::Result<Vec<u8>> {
    let mut sequence: Vec<u8> = record.sequence().iter().collect();
    sequence.make_ascii_uppercase();
    let qualities = record
        .quality_scores()
        .iter()
        .collect::<io::Result<Vec<u8>>>()?;
    if qualities.len() == sequence.len() {
        for (base, &q) in sequence.iter_mut().zip(&qualities) {
            if q < MIN_BASE_QUALITY {
                *base = b'N';
            }
        }
    }
    Ok(sequence)
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
    fn reads_of_the_alu_junctions() {
        // AluY: TSD chr1:10001-10012 (left junction 10012, right 10000).
        let request = JunctionRequest {
            window: ScanRegion {
                contig: "chr1".to_owned(),
                start: 9_600,
                end: 10_400,
            },
            left: Some(10_012),
            right: Some(10_000),
        };
        for (file, reference) in [
            ("positive.bam", None),
            ("positive.cram", Some(data("reference.fa"))),
        ] {
            let file = AlignmentFile::open(&data(file), reference.as_deref()).unwrap();
            let reads = fetch_junction_reads(
                &file,
                std::slice::from_ref(&request),
                &ReadFilter::default(),
            )
            .unwrap()
            .remove(0);
            assert!(reads.left.len() >= 20, "{}", reads.left.len());
            assert!(reads.right.len() >= 20, "{}", reads.right.len());
            // Reads of the poly(A) tail and the Alu are mapped to chr2 or
            // unmapped and placed at their mates.
            assert_ne!(reads.inside.len(), 0);
            // The left reads hold the end of the TSD and the Alu head.
            assert!(
                reads
                    .left
                    .iter()
                    .any(|r| r.windows(14).any(|w| w == b"CCCCCGGGCCGGGC"))
            );
        }
    }
}
