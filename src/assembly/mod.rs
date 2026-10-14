//! Phase 4: local assembly of the junctions, to the base.
//!
//! The reads that cross each junction of a candidate (see
//! [`crate::reader::junctions`]) are assembled from the reference flank
//! into the inserted element ([`graph`]): the left junction forward, from
//! the upstream flank; the right junction on the reverse complement, from
//! the downstream flank. The unmapped reads placed near the candidate (in
//! the element) are added in both orientations. Each contig is then compared
//! with the reference flank it starts from:
//!
//! - the **junction** is where the contig leaves the reference (a single
//!   mismatch followed by [`SNP_CHECK`] matching bases is a SNP or an error,
//!   not the junction); the rest of the contig is the start (left) or the
//!   end (right) of the inserted sequence, which must hold [`MIN_INSERT`]
//!   bases;
//! - with both junctions, the **TSD** is the reference between them when the
//!   flanks overlap by 1 to [`MAX_TSD`] bp, and there is none otherwise;
//! - the **poly(A/T) tail** is found at the inserted end next to a junction
//!   as in the clips; its length is the pure run of A (or T) from the
//!   junction, when the contig holds its end.
//!
//! When the inserted bases next to a junction match the reference beyond it
//! (microhomology), the contig matches the reference over them and the
//! junction is placed after them (before them on the right). The alternative
//! haplotype is the same whichever junction they are given to, so the
//! length of the TSD and of the tail cannot tell them apart.
//!
//! Coordinates: junctions are boundaries (see [`crate::reader::signals`]);
//! [`ReferenceWindow::start`] is 0-based.

pub mod graph;

use crate::cluster::tracts::ReferenceWindow;
use crate::reader::junctions::JunctionReads;
use crate::reader::signals::{ClipOptions, PolyTail, poly_tail};
use crate::typer::consensus::{ConsensusSet, Family};
use crate::typer::mates::without_runs;
use crate::typer::type_sequence;
use graph::{Graph, K, MIN_COUNT, reverse_complement};

/// Minimum number of inserted bases in a contig for its junction to count.
pub const MIN_INSERT: usize = 10;

/// Matching bases required after a mismatch for the contig to stay on the
/// reference flank.
pub const SNP_CHECK: usize = 8;

/// Distance, in bp, between the junction estimated from the clips and the
/// end of the seed taken from the reference: the closest seed held by
/// [`MIN_COUNT`] reads, from [`MIN_SEED_GAP`] to [`MAX_SEED_GAP`].
pub const MIN_SEED_GAP: u64 = 3;

/// See [`MIN_SEED_GAP`].
pub const MAX_SEED_GAP: u64 = 40;

/// Longest TSD.
pub const MAX_TSD: u64 = 50;

/// An assembled junction.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct JunctionContig {
    /// Junction (boundary): end of the upstream flank for the left one,
    /// start of the downstream flank for the right one.
    pub junction: u64,
    /// Contig, in reference orientation: flank then inserted bases (left),
    /// or inserted bases then flank (right).
    pub contig: Vec<u8>,
    /// Bases of the flank in the contig.
    pub flank: usize,
}

impl JunctionContig {
    /// Inserted bases of a left contig.
    #[must_use]
    pub fn left_insert(&self) -> &[u8] {
        &self.contig[self.flank..]
    }

    /// Inserted bases of a right contig.
    #[must_use]
    pub fn right_insert(&self) -> &[u8] {
        &self.contig[..self.contig.len() - self.flank]
    }
}

/// Result of the assembly of a candidate.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Assembly {
    /// Left junction, if assembled.
    pub left: Option<JunctionContig>,
    /// Right junction, if assembled.
    pub right: Option<JunctionContig>,
    /// With both junctions: the TSD (reference between them), `None` when the
    /// flanks do not overlap by 1 to [`MAX_TSD`] bp.
    pub tsd: Option<Vec<u8>>,
    /// Poly(A/T) tail at an assembled junction (its length: the pure run
    /// from the junction), when the contig holds its whole length.
    pub tail: Option<PolyTail>,
}

impl Assembly {
    /// Whether both junctions are assembled.
    #[must_use]
    pub fn both(&self) -> bool {
        self.left.is_some() && self.right.is_some()
    }

    /// `both`, `left`, `right` or `none`.
    #[must_use]
    pub fn assembled(&self) -> &'static str {
        match (&self.left, &self.right) {
            (Some(_), Some(_)) => "both",
            (Some(_), None) => "left",
            (None, Some(_)) => "right",
            (None, None) => "none",
        }
    }
}

/// Assemble the junctions of a candidate from its reads, the reference
/// around it and its junctions estimated from the clips.
#[must_use]
pub fn assemble(
    reads: &JunctionReads,
    reference: &ReferenceWindow<'_>,
    left: Option<u64>,
    right: Option<u64>,
    clip: &ClipOptions,
) -> Assembly {
    let inside: Vec<Vec<u8>> = reads
        .inside
        .iter()
        .flat_map(|r| [r.clone(), reverse_complement(r)])
        .collect();
    let left = left.and_then(|junction| {
        let mut walk: Vec<&[u8]> = reads.left.iter().map(Vec::as_slice).collect();
        walk.extend(inside.iter().map(Vec::as_slice));
        assemble_left(&Graph::new(walk), reference, junction)
    });
    let right = right.and_then(|junction| {
        let reverse: Vec<Vec<u8>> = reads.right.iter().map(|r| reverse_complement(r)).collect();
        let mut walk: Vec<&[u8]> = reverse.iter().map(Vec::as_slice).collect();
        walk.extend(inside.iter().map(Vec::as_slice));
        assemble_right(&Graph::new(walk), reference, junction)
    });

    let tsd = match (&left, &right) {
        (Some(l), Some(r)) if l.junction > r.junction && l.junction - r.junction <= MAX_TSD => {
            reference_bases(reference, r.junction, l.junction).map(<[u8]>::to_vec)
        }
        _ => None,
    };
    // Tails next to the junctions: at the start of the left insert, at the
    // end of the right one. A tail is found as in the clips (sequencing
    // errors tolerated), and its length is the pure run from the junction
    // (the contig is a consensus); only when the contig goes beyond it.
    let tail_of = |insert: &[u8], at_end: bool| {
        let tail = poly_tail(insert, at_end, clip)?;
        let run = if at_end {
            insert.iter().rev().take_while(|&&b| b == tail.base).count()
        } else {
            insert.iter().take_while(|&&b| b == tail.base).count()
        };
        (run < insert.len()).then_some(PolyTail {
            base: tail.base,
            length: u32::try_from(run).ok()?,
        })
    };
    let tails = [
        left.as_ref().and_then(|l| tail_of(l.left_insert(), false)),
        right.as_ref().and_then(|r| tail_of(r.right_insert(), true)),
    ];
    let tail = tails.into_iter().flatten().max_by_key(|t| t.length);
    Assembly {
        left,
        right,
        tsd,
        tail,
    }
}

/// Assemble a candidate placed by its mates, which has no junction from its
/// clips (#66): the breakpoint (a boundary where reads are clipped by a few
/// bases) does not say which flank the clipped reads leave, so it is tried as
/// the left and as the right junction, each alone. A side is kept if its
/// inserted bases hold a poly(A/T) tail next to the junction or are an
/// element of the `family` of the candidate (short exome clips often reach
/// only the tail and a few bases of the element); the sides kept are then
/// assembled together (TSD).
#[must_use]
pub fn assemble_at_breakpoint(
    reads: &JunctionReads,
    reference: &ReferenceWindow<'_>,
    breakpoint: u64,
    family: Family,
    set: &ConsensusSet,
    clip: &ClipOptions,
) -> Assembly {
    let is_element = |insert: &[u8]| {
        type_sequence(set, &without_runs(insert))
            .iter()
            .any(|hit| set.sequences()[hit.consensus].family == family)
    };
    let from_left = assemble(reads, reference, Some(breakpoint), None, clip);
    let from_right = assemble(reads, reference, None, Some(breakpoint), clip);
    let keep_left = from_left
        .left
        .as_ref()
        .is_some_and(|l| from_left.tail.is_some() || is_element(l.left_insert()));
    let keep_right = from_right
        .right
        .as_ref()
        .is_some_and(|r| from_right.tail.is_some() || is_element(r.right_insert()));
    match (keep_left, keep_right) {
        (true, true) => assemble(reads, reference, Some(breakpoint), Some(breakpoint), clip),
        (true, false) => from_left,
        (false, true) => from_right,
        (false, false) => Assembly::default(),
    }
}

/// Reference bases of a 0-based half-open interval, if in the window.
fn reference_bases<'a>(reference: &ReferenceWindow<'a>, start: u64, end: u64) -> Option<&'a [u8]> {
    let from = usize::try_from(start.checked_sub(reference.start)?).ok()?;
    let to = usize::try_from(end.checked_sub(reference.start)?).ok()?;
    reference.sequence.get(from..to)
}

/// Assemble the left junction: from a seed of the upstream flank, forward.
fn assemble_left(
    graph: &Graph<'_>,
    reference: &ReferenceWindow<'_>,
    junction: u64,
) -> Option<JunctionContig> {
    let k = K as u64;
    let (seed_start, seed) = (MIN_SEED_GAP..=MAX_SEED_GAP).find_map(|gap| {
        let end = junction.checked_sub(gap)?;
        let seed = reference_bases(reference, end.checked_sub(k)?, end)?;
        (graph.count(seed) >= MIN_COUNT).then_some((end - k, seed))
    })?;
    let contig = graph.extend(seed);
    let flank_reference = reference_bases(
        reference,
        seed_start,
        reference.start + reference.sequence.len() as u64,
    )?;
    let flank = flank_match(contig.iter().copied(), flank_reference.iter().copied())?;
    (contig.len() - flank >= MIN_INSERT).then(|| JunctionContig {
        junction: seed_start + flank as u64,
        contig,
        flank,
    })
}

/// Assemble the right junction: from a seed of the downstream flank, on the
/// reverse complement.
fn assemble_right(
    graph: &Graph<'_>,
    reference: &ReferenceWindow<'_>,
    junction: u64,
) -> Option<JunctionContig> {
    let k = K as u64;
    let (seed_end, seed) = (MIN_SEED_GAP..=MAX_SEED_GAP).find_map(|gap| {
        let start = junction + gap;
        let seed = reverse_complement(reference_bases(reference, start, start + k)?);
        (graph.count(&seed) >= MIN_COUNT).then_some((start + k, seed))
    })?;
    let contig = graph.extend(&seed);
    let flank_reference =
        reverse_complement(reference_bases(reference, reference.start, seed_end)?);
    let flank = flank_match(contig.iter().copied(), flank_reference.iter().copied())?;
    (contig.len() - flank >= MIN_INSERT).then(|| JunctionContig {
        junction: seed_end - flank as u64,
        contig: reverse_complement(&contig),
        flank,
    })
}

/// Number of leading bases of a contig on the reference flank: up to the
/// first mismatch not followed by [`SNP_CHECK`] matching bases. `None` if
/// the contig never leaves the flank (or the window ends first).
fn flank_match(
    contig: impl Iterator<Item = u8>,
    reference: impl Iterator<Item = u8>,
) -> Option<usize> {
    let pairs: Vec<(u8, u8)> = contig.zip(reference).collect();
    let mut i = 0;
    while i < pairs.len() {
        let (c, r) = pairs[i];
        if c != r {
            let next = &pairs[i + 1..pairs.len().min(i + 1 + SNP_CHECK)];
            if next.len() < SNP_CHECK || next.iter().any(|(c, r)| c != r) {
                return Some(i);
            }
        }
        i += 1;
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reads of `length` bases tiling `sequence` every base.
    fn tiles(sequence: &[u8], length: usize) -> Vec<Vec<u8>> {
        (0..=sequence.len() - length)
            .map(|i| sequence[i..i + length].to_vec())
            .collect()
    }

    /// A reference of 400 bases (no repeat of 21 bases), an insertion with
    /// a TSD of `tsd` bases ending at 200, and its junction reads.
    fn insertion(element: &[u8], tsd: usize) -> (Vec<u8>, JunctionReads) {
        let mut state = 7u64;
        let reference: Vec<u8> = (0..400)
            .map(|_| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1);
                b"ACGT"[(state >> 33) as usize % 4]
            })
            .collect();
        let (right, left) = (200 - tsd, 200);
        let alt = [&reference[..left], element, &reference[right..]].concat();
        let reads = tiles(&alt, 70);
        // Reads crossing each junction, as the reader gathers them.
        let left_reads = reads[left - 60..left - 10].to_vec();
        let right_end = left + element.len();
        let right_reads = reads[right_end - 60..right_end - 10].to_vec();
        (
            reference,
            JunctionReads {
                left: left_reads,
                right: right_reads,
                inside: Vec::new(),
            },
        )
    }

    const ELEMENT: &[u8] =
        b"GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAG";

    #[test]
    fn junctions_tsd_and_tail() {
        let element = [ELEMENT, &[b'A'; 30][..]].concat();
        let (reference, reads) = insertion(&element, 12);
        let window = ReferenceWindow {
            start: 1_000,
            sequence: &reference,
        };
        // Junctions estimated 2 bp off.
        let a = assemble(
            &reads,
            &window,
            Some(1_202),
            Some(1_186),
            &ClipOptions::default(),
        );
        let (left, right) = (a.left.as_ref().unwrap(), a.right.as_ref().unwrap());
        assert_eq!((left.junction, right.junction), (1_200, 1_188));
        assert_eq!(a.tsd.as_deref(), Some(&reference[188..200]));
        assert_eq!(
            a.tail,
            Some(PolyTail {
                base: b'A',
                length: 30
            })
        );
        assert!(left.left_insert().starts_with(&ELEMENT[..20]));
        assert!(right.right_insert().ends_with(&[b'A'; 30]));
        assert_eq!(a.assembled(), "both");
    }

    #[test]
    fn no_reads_no_junction() {
        let (reference, mut reads) = insertion(ELEMENT, 12);
        reads.right.clear();
        let window = ReferenceWindow {
            start: 1_000,
            sequence: &reference,
        };
        let a = assemble(
            &reads,
            &window,
            Some(1_200),
            Some(1_188),
            &ClipOptions::default(),
        );
        assert!(a.left.is_some() && a.right.is_none());
        assert_eq!((a.tsd.as_deref(), a.assembled()), (None, "left"));
    }

    #[test]
    fn breakpoint_without_clip_junctions() {
        // Only the reads of the left junction exist; the breakpoint is tried
        // on both sides and the left contig is an Alu.
        let (reference, mut reads) = insertion(ELEMENT, 12);
        reads.right.clear();
        let window = ReferenceWindow {
            start: 1_000,
            sequence: &reference,
        };
        let set = ConsensusSet::load(None).unwrap();
        let clip = ClipOptions::default();
        let a = assemble_at_breakpoint(&reads, &window, 1_200, Family::Alu, &set, &clip);
        assert_eq!(a.assembled(), "left");
        assert_eq!(a.left.as_ref().map(|l| l.junction), Some(1_200));
        // Not an element of another family (and no tail): no contig kept.
        let a = assemble_at_breakpoint(&reads, &window, 1_200, Family::Line1, &set, &clip);
        assert_eq!(a.assembled(), "none");
    }

    #[test]
    fn flank_matches() {
        let m = |c: &[u8], r: &[u8]| flank_match(c.iter().copied(), r.iter().copied());
        assert_eq!(m(b"ACGTAC", b"ACGTTT"), Some(4));
        // A SNP followed by enough matching bases stays on the flank.
        assert_eq!(m(b"ACGTACGTACGTAG", b"ACCTACGTACGTCC"), Some(12));
        // The contig never leaves the reference.
        assert_eq!(m(b"ACGT", b"ACGTAA"), None);
    }
}
