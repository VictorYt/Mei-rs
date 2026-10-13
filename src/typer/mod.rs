//! Phase 3: typing of the inserted element (Alu, LINE-1, SVA) from the
//! soft-clipped sequences of a candidate.
//!
//! Each junction soft-clip, without the pure run of its poly(A/T) tail, is
//! seeded against the k-mer index of the consensus sequences (both strands,
//! [`consensus`]); the consensus strands with the most seeds are aligned
//! locally ([`align`]) in a window around the seeds. A clip is typed by the
//! alignments covering at least [`MIN_COVERAGE`] of it with at least
//! [`MIN_IDENTITY`].
//!
//! A candidate takes the family with the highest summed score over its
//! clips, the subfamily (consensus) with the highest summed score in that
//! family, the strand of the element (`+` when the clips match the forward
//! consensus, in reference orientation) and the part of the consensus
//! covered. Close subfamilies (e.g. `AluY` and `AluYa5`) differ by a few
//! bases: short clips often cannot tell them apart, and the first one in the
//! consensus file wins ties.
//!
//! The alignment is scalar: a clip of at most a read length is aligned
//! against a few windows of a few hundred bases, and only the junction clips
//! of candidates are typed.

pub mod align;
pub mod consensus;
pub mod mates;

use std::collections::HashMap;

use crate::cluster::{Candidate, JUNCTION_TOLERANCE};
use crate::reader::signals::{Side, SignalKind};
use consensus::{ConsensusSet, Family, K, kmers};

/// Minimum length of a clip, without the pure run of its tail, to be typed.
pub const MIN_TYPED_LENGTH: usize = 15;

/// Number of consensus strands aligned per clip (those with the most seeds).
pub const ALIGNED_STRANDS: usize = 8;

/// Bases added on each side of the seeds to the aligned consensus window.
pub const WINDOW_MARGIN: usize = 30;

/// Minimum fraction of the clip (without its tail run) in the alignment.
pub const MIN_COVERAGE: f64 = 0.7;

/// Minimum identity of the alignment.
pub const MIN_IDENTITY: f64 = 0.9;

/// Alignment of one clip to one consensus.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ClipHit {
    /// Index of the consensus.
    pub consensus: usize,
    /// The clip matches the reverse complement of the consensus.
    pub reverse: bool,
    /// Alignment score.
    pub score: i32,
    /// Identity of the alignment.
    pub identity: f64,
    /// Aligned part of the forward consensus, 0-based half-open.
    pub element_start: usize,
    /// End of the aligned part of the forward consensus (exclusive).
    pub element_end: usize,
}

/// Element typed from the clips of a candidate.
#[derive(Clone, Debug, PartialEq)]
pub struct Typing {
    /// Family.
    pub family: Family,
    /// Best consensus of the family (subfamily).
    pub subfamily: String,
    /// The element lies on the reverse strand of the reference.
    pub reverse: bool,
    /// Part of the subfamily consensus covered by its clips, 0-based
    /// half-open.
    pub element_start: usize,
    /// End of the covered part (exclusive).
    pub element_end: usize,
    /// Mean identity of the clips aligned to the subfamily consensus.
    pub identity: f64,
    /// Junction clips typed to the family.
    pub typed_clips: usize,
    /// Of these, short clips (see
    /// [`crate::reader::signals::ClipOptions::is_short`]).
    pub short_typed_clips: usize,
    /// Typed from the mates, no junction clip being typed
    /// ([`mates::type_mates`]).
    pub by_mates: bool,
}

/// Alignments of a sequence (reference orientation) that pass
/// [`MIN_COVERAGE`] and [`MIN_IDENTITY`], at most one per consensus strand,
/// best score first.
#[must_use]
pub fn type_sequence(set: &ConsensusSet, sequence: &[u8]) -> Vec<ClipHit> {
    type_part(set, sequence, None)
}

/// [`type_sequence`], with alignments of at least `min_aligned` bases of the
/// sequence accepted instead of [`MIN_COVERAGE`] when given: for a sequence
/// of which only a part is expected to be the element.
#[must_use]
pub fn type_part(set: &ConsensusSet, sequence: &[u8], min_aligned: Option<usize>) -> Vec<ClipHit> {
    if sequence.len() < MIN_TYPED_LENGTH.max(K) {
        return Vec::new();
    }
    // Seed diagonals (consensus position - clip position) per strand.
    let mut seeds: HashMap<(u32, bool), Vec<i64>> = HashMap::new();
    for (position, code) in kmers(sequence) {
        for occurrence in set.occurrences(code) {
            seeds
                .entry((occurrence.consensus, occurrence.reverse))
                .or_default()
                .push(i64::from(occurrence.position) - i64::try_from(position).unwrap_or(0));
        }
    }
    // Rank the strands by their seeds on one diagonal band; two seeds at
    // least, unless the clip is too short to hold two.
    let min_seeds = if sequence.len() >= 2 * K + 3 { 2 } else { 1 };
    let mut strands: Vec<(Strand, Band)> = seeds
        .into_iter()
        .map(|(strand, mut diagonals)| (strand, best_band(&mut diagonals)))
        .filter(|(_, (count, _, _))| *count >= min_seeds)
        .collect();
    strands.sort_by(|a, b| b.1.0.cmp(&a.1.0).then(a.0.cmp(&b.0)));
    strands.truncate(ALIGNED_STRANDS);

    let length = i64::try_from(sequence.len()).unwrap_or(i64::MAX);
    let margin = i64::try_from(WINDOW_MARGIN).unwrap_or(0);
    let mut hits: Vec<ClipHit> = strands
        .into_iter()
        .filter_map(|((consensus, reverse), (_, first, last))| {
            let consensus = consensus as usize;
            let strand = set.strand(consensus, reverse);
            let strand_length = i64::try_from(strand.len()).unwrap_or(i64::MAX);
            let lo = usize::try_from((first - margin).max(0)).unwrap_or(0);
            let hi = usize::try_from((last + length + margin).min(strand_length)).unwrap_or(0);
            let alignment = align::local(sequence, strand.get(lo..hi)?)?;
            let aligned = alignment.query_end - alignment.query_start;
            #[allow(clippy::cast_precision_loss)]
            let covered = min_aligned.map_or_else(
                || aligned as f64 / sequence.len() as f64 >= MIN_COVERAGE,
                |min| aligned >= min,
            );
            if !covered || alignment.identity() < MIN_IDENTITY {
                return None;
            }
            let (start, end) = (lo + alignment.target_start, lo + alignment.target_end);
            let (element_start, element_end) = if reverse {
                (strand.len() - end, strand.len() - start)
            } else {
                (start, end)
            };
            Some(ClipHit {
                consensus,
                reverse,
                score: alignment.score,
                identity: alignment.identity(),
                element_start,
                element_end,
            })
        })
        .collect();
    hits.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then(a.consensus.cmp(&b.consensus))
            .then(a.reverse.cmp(&b.reverse))
    });
    // One hit per consensus: its best strand.
    let mut seen = Vec::new();
    hits.retain(|h| {
        let new = !seen.contains(&h.consensus);
        seen.push(h.consensus);
        new
    });
    hits
}

/// A consensus strand: consensus index and reverse strand.
type Strand = (u32, bool);

/// Seeds of a diagonal band: count, first and last diagonal.
type Band = (usize, i64, i64);

/// Width of the diagonal band in which seeds are counted (indels shift the
/// diagonal).
pub const BAND: i64 = 20;

/// The band of [`BAND`] diagonals holding the most seeds: seed count, first
/// and last diagonal of the band. Sorts `diagonals`.
fn best_band(diagonals: &mut [i64]) -> Band {
    diagonals.sort_unstable();
    let mut best = (0, 0, 0);
    let mut first = 0;
    for last in 0..diagonals.len() {
        while diagonals[last] - diagonals[first] > BAND {
            first += 1;
        }
        let count = last - first + 1;
        if count > best.0 {
            best = (count, diagonals[first], diagonals[last]);
        }
    }
    best
}

/// Whether a soft-clip is evidence of a mobile element on its own: it has a
/// poly(A/T) tail, or its sequence is typed to an MEI consensus
/// ([`type_sequence`]). Short clips must be (see
/// [`crate::reader::signals::ClipOptions::is_short`]): a random 15 bp
/// sequence is typed about once in a hundred times with the default set,
/// and short clips are far more frequent than long ones. Other signals are
/// not soft-clips: `false`.
#[must_use]
pub fn is_mei_clip(set: &ConsensusSet, kind: &SignalKind) -> bool {
    match kind {
        SignalKind::SoftClip { tail: Some(_), .. } => true,
        SignalKind::SoftClip { sequence, .. } => !type_sequence(set, sequence).is_empty(),
        _ => false,
    }
}

/// Type the element inserted at a candidate from its junction soft-clips;
/// `None` if no clip is typed.
#[must_use]
pub fn type_candidate(set: &ConsensusSet, candidate: &Candidate) -> Option<Typing> {
    let near = |junction: Option<u64>, position: u64| {
        junction.is_some_and(|j| j.abs_diff(position) <= JUNCTION_TOLERANCE)
    };
    // Hits of each typed clip, and whether it is short.
    let typed: Vec<(Vec<ClipHit>, bool)> = candidate
        .signals
        .iter()
        .filter_map(|signal| {
            let SignalKind::SoftClip {
                sequence,
                tail,
                short,
                ..
            } = &signal.kind
            else {
                return None;
            };
            let junction = match signal.side {
                Side::Right => candidate.left_junction,
                Side::Left => candidate.right_junction,
            };
            if !near(junction, signal.position) {
                return None;
            }
            // Without the pure run of the tail, at the junction end.
            let run = tail.map_or(0, |t| match signal.side {
                Side::Left => sequence.iter().rev().take_while(|&&b| b == t.base).count(),
                Side::Right => sequence.iter().take_while(|&&b| b == t.base).count(),
            });
            let rest = match signal.side {
                Side::Left => &sequence[..sequence.len() - run],
                Side::Right => &sequence[run..],
            };
            let hits = type_sequence(set, rest);
            (!hits.is_empty()).then_some((hits, *short))
        })
        .collect();
    if typed.is_empty() {
        return None;
    }

    let family_of = |hit: &ClipHit| set.sequences()[hit.consensus].family;
    // Best hit of each family per clip, summed per family, with the number
    // of clips and of short clips.
    let mut families: HashMap<Family, (i64, usize, usize)> = HashMap::new();
    for (hits, short) in &typed {
        let mut best: HashMap<Family, i32> = HashMap::new();
        for hit in hits {
            let score = best.entry(family_of(hit)).or_insert(hit.score);
            *score = (*score).max(hit.score);
        }
        for (family, score) in best {
            let entry = families.entry(family).or_default();
            entry.0 += i64::from(score);
            entry.1 += 1;
            entry.2 += usize::from(*short);
        }
    }
    let (&family, &(_, typed_clips, short_typed_clips)) = families
        .iter()
        .max_by(|a, b| a.1.0.cmp(&b.1.0).then(b.0.cmp(a.0)))?;

    // Subfamily: the consensus of the family with the highest summed score.
    let mut scores: HashMap<usize, i64> = HashMap::new();
    for hit in typed
        .iter()
        .flat_map(|(hits, _)| hits)
        .filter(|h| family_of(h) == family)
    {
        *scores.entry(hit.consensus).or_default() += i64::from(hit.score);
    }
    let (&best, _) = scores
        .iter()
        .max_by(|a, b| a.1.cmp(b.1).then(b.0.cmp(a.0)))?;
    let on_best: Vec<&ClipHit> = typed
        .iter()
        .flat_map(|(hits, _)| hits)
        .filter(|h| h.consensus == best)
        .collect();
    let reverse_score: i64 = on_best
        .iter()
        .map(|h| if h.reverse { 1 } else { -1 } * i64::from(h.score))
        .sum();
    #[allow(clippy::cast_precision_loss)]
    let identity = on_best.iter().map(|h| h.identity).sum::<f64>() / on_best.len() as f64;
    Some(Typing {
        family,
        subfamily: set.sequences()[best].name.clone(),
        reverse: reverse_score > 0,
        element_start: on_best.iter().map(|h| h.element_start).min()?,
        element_end: on_best.iter().map(|h| h.element_end).max()?,
        identity,
        typed_clips,
        short_typed_clips,
        by_mates: false,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cluster::{ClusterOptions, cluster};
    use crate::reader::signals::{PolyTail, Signal};
    use consensus::reverse_complement;

    fn set() -> ConsensusSet {
        ConsensusSet::load(None).unwrap()
    }

    fn consensus(set: &ConsensusSet, name: &str) -> Vec<u8> {
        set.sequences()
            .iter()
            .find(|c| c.name == name)
            .unwrap()
            .sequence
            .clone()
    }

    fn clip(position: u64, side: Side, read: &str, sequence: &[u8], tail: Option<u8>) -> Signal {
        let mut signal = crate::cluster::tests::clip(position, side, read, None);
        let length = tail.map_or(0, |base| match side {
            Side::Left => sequence.iter().rev().take_while(|&&b| b == base).count(),
            Side::Right => sequence.iter().take_while(|&&b| b == base).count(),
        });
        signal.kind = SignalKind::SoftClip {
            sequence: sequence.to_vec(),
            median_quality: 30,
            tail: tail.map(|base| PolyTail {
                base,
                length: u32::try_from(length).unwrap(),
            }),
            supplementary: false,
            short: false,
        };
        signal
    }

    fn candidate(signals: Vec<Signal>) -> Candidate {
        let mut signals = signals;
        signals.sort();
        cluster(
            signals,
            &ClusterOptions {
                require_tail: false,
                ..ClusterOptions::default()
            },
        )
        .remove(0)
    }

    #[test]
    fn alu_clip_on_both_strands() {
        let set = set();
        let alu = consensus(&set, "AluY");
        let piece = &alu[100..160];
        let hits = type_sequence(&set, piece);
        let best = hits[0];
        assert_eq!(set.sequences()[best.consensus].family, Family::Alu);
        assert!(!best.reverse);
        assert_eq!((best.element_start, best.element_end), (100, 160));
        assert!((best.identity - 1.0).abs() < 1e-9);

        let hits = type_sequence(&set, &reverse_complement(piece));
        assert!(hits[0].reverse);
        assert_eq!((hits[0].element_start, hits[0].element_end), (100, 160));
    }

    #[test]
    fn line1_and_sva_clips() {
        let set = set();
        let l1 = consensus(&set, "L1HS_3end");
        let hits = type_sequence(&set, &l1[400..450]);
        assert_eq!(set.sequences()[hits[0].consensus].family, Family::Line1);
        let sva = consensus(&set, "SVA_F");
        let hits = type_sequence(&set, &sva[700..760]);
        assert_eq!(set.sequences()[hits[0].consensus].family, Family::Sva);
    }

    #[test]
    fn line1_orf2_clips() {
        // The middle of LINE-1 (ORF2), not in the L1HS 5' and 3' end models:
        // a 5'-truncated LINE-1 can start there. L1HS differs from the L1P1
        // ORF2 model by a few percent.
        let set = set();
        let mut piece = consensus(&set, "L1P1_orf2")[1400..1460].to_vec();
        for i in (10..60).step_by(20) {
            piece[i] = if piece[i] == b'A' { b'G' } else { b'A' };
        }
        let hits = type_sequence(&set, &reverse_complement(&piece));
        let best = hits[0];
        assert_eq!(set.sequences()[best.consensus].family, Family::Line1);
        assert_eq!(set.sequences()[best.consensus].name, "L1P1_orf2");
        assert!(best.reverse);
        assert_eq!((best.element_start, best.element_end), (1400, 1460));
    }

    #[test]
    fn mismatches_are_tolerated_not_random_sequence() {
        let set = set();
        let mut piece = consensus(&set, "AluY")[30..90].to_vec();
        for i in (5..60).step_by(12) {
            piece[i] = if piece[i] == b'A' { b'C' } else { b'A' };
        }
        let hits = type_sequence(&set, &piece);
        assert_eq!(set.sequences()[hits[0].consensus].family, Family::Alu);
        assert!(hits[0].identity < 1.0 && hits[0].identity >= MIN_IDENTITY);

        let random = b"TGCATGCTAGCTAGGATCGATCGTAGCTAGCTTTAGGCGATCGAT";
        assert_eq!(type_sequence(&set, random), []);
        // Too short, or only a tail.
        assert_eq!(type_sequence(&set, &consensus(&set, "AluY")[30..42]), []);
        assert_eq!(type_sequence(&set, &[b'A'; 40]), []);
    }

    #[test]
    fn partly_aligned_sequence() {
        // A mate aligned in part to a reference Alu: the reference read over
        // a read length is 34 Alu bases, then unrelated sequence.
        let set = set();
        let mut sequence = consensus(&set, "AluY")[240..274].to_vec();
        sequence.extend_from_slice(b"TGCATGCTAGCTAGGATCGATCGTAGCTAGCTTTAGGCGA");
        assert_eq!(type_sequence(&set, &sequence), []);
        let hits = type_part(&set, &sequence, Some(30));
        assert_eq!(set.sequences()[hits[0].consensus].family, Family::Alu);
        assert_eq!(type_part(&set, &sequence, Some(40)), []);
    }

    #[test]
    fn mei_clips() {
        let set = set();
        let soft_clip = |sequence: &[u8], tail| SignalKind::SoftClip {
            sequence: sequence.to_vec(),
            median_quality: 30,
            tail,
            supplementary: false,
            short: false,
        };
        let alu = consensus(&set, "AluY");
        assert!(is_mei_clip(&set, &soft_clip(&alu[120..136], None)));
        assert!(!is_mei_clip(&set, &soft_clip(b"TGCATGCTAGCTAGGA", None)));
        let tail = PolyTail {
            base: b'A',
            length: 12,
        };
        assert!(is_mei_clip(
            &set,
            &soft_clip(b"AAAAAAAAAAAAGTC", Some(tail))
        ));
        assert!(!is_mei_clip(&set, &SignalKind::MateUnmapped));
    }

    #[test]
    fn band() {
        let mut d = vec![50, 3, 0, 1, 2, 70, 60, 2];
        assert_eq!(best_band(&mut d), (5, 0, 3));
        assert_eq!(best_band(&mut []), (0, 0, 0));
    }

    #[test]
    fn candidate_typing() {
        let set = set();
        let alu = consensus(&set, "AluY");
        // Plus-strand Alu: clips of the left junction start with the Alu
        // head, clips of the right junction end with its 3' end and a tail.
        let mut tailed = alu[250..290].to_vec();
        tailed.extend_from_slice(&[b'A'; 20]);
        let c = candidate(vec![
            clip(1_012, Side::Right, "l1", &alu[..50], None),
            clip(1_012, Side::Right, "l2", &alu[..35], None),
            clip(1_000, Side::Left, "r1", &tailed, Some(b'A')),
            clip(1_000, Side::Left, "r2", &tailed[10..], Some(b'A')),
        ]);
        let typing = type_candidate(&set, &c).unwrap();
        assert_eq!(typing.family, Family::Alu);
        assert!(!typing.reverse);
        assert_eq!(typing.typed_clips, 4);
        assert_eq!(typing.subfamily, "AluY");
        // Up to the A run of the consensus, set aside with the tail.
        assert_eq!((typing.element_start, typing.element_end), (0, 281));

        // The same insertion on the minus strand.
        let c = candidate(vec![
            clip(
                1_012,
                Side::Right,
                "l1",
                &reverse_complement(&tailed),
                Some(b'T'),
            ),
            clip(
                1_000,
                Side::Left,
                "r1",
                &reverse_complement(&alu[..50]),
                None,
            ),
        ]);
        let typing = type_candidate(&set, &c).unwrap();
        assert_eq!(typing.family, Family::Alu);
        assert!(typing.reverse);
        assert_eq!(typing.typed_clips, 2);
    }

    #[test]
    fn untyped_candidate() {
        let set = set();
        let c = candidate(vec![
            clip(
                1_000,
                Side::Left,
                "r1",
                b"TGCATGCTAGCTAGGATCGATCGTAGCTAG",
                None,
            ),
            clip(
                1_000,
                Side::Left,
                "r2",
                b"GCATGCTAGCTAGGATCGATCGTAGCTAG",
                None,
            ),
        ]);
        assert_eq!(type_candidate(&set, &c), None);
    }
}
