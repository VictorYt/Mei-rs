//! MEI mates: mates whose sequence is a mobile element.
//!
//! A mate signal says that the mate of an anchored read lies in the inserted
//! sequence; it is evidence of a mobile element only if the mate is one.
//! Up to [`MAX_CHECKED_MATES`] mates of each typed candidate with enough
//! junction soft-clips (the others cannot pass) are sampled
//! ([`sample_mates`]), and their sequence is found:
//!
//! - unmapped mates: their records are placed at the anchor; they are read
//!   again at the anchors of the candidate (see [`crate::reader::mates`]);
//! - mapped mates (low MAPQ, another contig, large insert): the mate aligned
//!   at its locus, so the reference there stands for its sequence; it is read
//!   from `--reference`, one query per mate ([`fetch_mapped_mates`]; without
//!   a reference, mapped mates are not checked).
//!
//! A mate is an MEI mate when its best alignment ([`super::type_sequence`])
//! is in the family of the candidate, or in any family for an untyped one.
//!
//! The checks are kept in the candidate ([`Candidate::mate_checks`]), not in
//! the signals: millions of signals, a few thousand checked candidates.
//!
//! Coordinates: [`Locus::start`] is 0-based; FASTA regions are 1-based
//! inclusive.

use std::io;
use std::path::Path;

use noodles::core::{Position, Region};
use noodles::fasta;

use super::consensus::ConsensusSet;
use super::type_sequence;
use crate::cluster::{Candidate, Filter};
use crate::reader::mates::{MateRequest, name_hash};
use crate::reader::signals::{Locus, SignalKind};
use crate::regions::targets::Contig;

/// Maximum number of mates checked per candidate, sampled evenly along it:
/// the fraction of MEI mates does not need them all, and deep samples have
/// candidates with thousands of mate signals.
pub const MAX_CHECKED_MATES: usize = 40;

/// Bases of reference read at the locus of a mapped mate when the length of
/// the anchored read is shorter.
pub const MIN_MATE_LENGTH: u64 = 50;

/// A mate sampled to be checked.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MateCheck {
    /// Index of the mate signal in [`Candidate::signals`].
    pub signal: usize,
    /// Sequence of the mate, once found.
    pub sequence: Option<Vec<u8>>,
    /// Whether it is a mobile element of the candidate family, once checked.
    pub mei: Option<bool>,
}

/// Whether the mates of a candidate are checked: it is typed and has enough
/// junction soft-clips (`no_junction` is decided at clustering), so it can
/// pass.
#[must_use]
pub fn checked(candidate: &Candidate) -> bool {
    candidate.typing.is_some() && !candidate.filters.contains(&Filter::NoJunction)
}

/// Sample at most [`MAX_CHECKED_MATES`] mates of a checked candidate, evenly
/// spaced among its unmapped and mapped mate signals, into
/// [`Candidate::mate_checks`].
pub fn sample_mates(candidate: &mut Candidate) {
    candidate.mate_checks.clear();
    if !checked(candidate) {
        return;
    }
    let eligible: Vec<usize> = (0..candidate.signals.len())
        .filter(|&i| {
            let kind = &candidate.signals[i].kind;
            *kind == SignalKind::MateUnmapped || kind.mate_locus().is_some()
        })
        .collect();
    let n = eligible.len();
    // The first rank of each of MAX_CHECKED_MATES equal slices of the ranks.
    let slice = |rank: usize| rank * MAX_CHECKED_MATES / n;
    candidate.mate_checks = eligible
        .iter()
        .enumerate()
        .filter(|&(rank, _)| n <= MAX_CHECKED_MATES || rank == 0 || slice(rank) != slice(rank - 1))
        .map(|(_, &signal)| MateCheck {
            signal,
            sequence: None,
            mei: None,
        })
        .collect();
}

/// The unmapped mates to read: one request group per candidate (with its
/// index), each request with the index of its mate check.
#[must_use]
pub fn unmapped_requests(candidates: &[Candidate]) -> Vec<(usize, Vec<(usize, MateRequest)>)> {
    candidates
        .iter()
        .enumerate()
        .filter_map(|(c, candidate)| {
            let requests: Vec<(usize, MateRequest)> = candidate
                .mate_checks
                .iter()
                .enumerate()
                .filter_map(|(k, check)| {
                    let signal = &candidate.signals[check.signal];
                    (signal.kind == SignalKind::MateUnmapped).then(|| {
                        (
                            k,
                            MateRequest {
                                contig: signal.contig,
                                anchor_start: signal.start,
                                name: name_hash(signal.read.as_bytes()),
                            },
                        )
                    })
                })
                .collect();
            (!requests.is_empty()).then_some((c, requests))
        })
        .collect()
}

/// Read the reference at the locus of the sampled mapped mates, as their
/// sequence (a read length from the mate start, the length of the anchored
/// alignment).
///
/// # Errors
///
/// Fails if the reference cannot be opened or read.
pub fn fetch_mapped_mates(
    candidates: &mut [Candidate],
    contigs: &[Contig],
    reference: &Path,
) -> io::Result<()> {
    let mut reader = fasta::io::indexed_reader::Builder::default().build_from_path(reference)?;
    // (candidate, check, locus, 0-based end), sorted by locus for sequential
    // reads.
    let mut loci: Vec<(usize, usize, Locus, u64)> = Vec::new();
    for (c, candidate) in candidates.iter().enumerate() {
        for (k, check) in candidate.mate_checks.iter().enumerate() {
            let signal = &candidate.signals[check.signal];
            let Some(locus) = signal.kind.mate_locus() else {
                continue;
            };
            let Some(length) = contigs.get(locus.contig).map(|c| c.length) else {
                continue;
            };
            let read_length = (signal.end - signal.start).max(MIN_MATE_LENGTH);
            let end = (locus.start + read_length).min(length);
            if locus.start < end {
                loci.push((c, k, locus, end));
            }
        }
    }
    loci.sort_unstable_by_key(|&(_, _, locus, _)| locus);
    for (c, k, locus, end) in loci {
        let (Some(first), Some(last)) = (
            usize::try_from(locus.start + 1)
                .ok()
                .and_then(|p| Position::try_from(p).ok()),
            usize::try_from(end)
                .ok()
                .and_then(|p| Position::try_from(p).ok()),
        ) else {
            continue;
        };
        let region = Region::new(contigs[locus.contig].name.as_str(), first..=last);
        let record = reader.query(&region)?;
        let mut sequence = record.sequence().as_ref().to_vec();
        sequence.make_ascii_uppercase();
        candidates[c].mate_checks[k].sequence = Some(sequence);
    }
    Ok(())
}

/// Type the sampled mates of a candidate with a sequence: `checked_mates` and
/// `mei_mates` (typed to its family, or to any family if it is untyped).
pub fn type_mates(set: &ConsensusSet, candidate: &mut Candidate) {
    let family = candidate.typing.as_ref().map(|t| t.family);
    let (mut checked, mut mei) = (0, 0);
    for check in &mut candidate.mate_checks {
        let Some(sequence) = check.sequence.as_deref() else {
            continue;
        };
        checked += 1;
        let best = type_sequence(set, sequence)
            .first()
            .map(|hit| set.sequences()[hit.consensus].family);
        let is_mei = best.is_some_and(|f| family.is_none_or(|family| f == family));
        check.mei = Some(is_mei);
        mei += usize::from(is_mei);
    }
    candidate.checked_mates = checked;
    candidate.mei_mates = mei;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cluster::tests::mate;
    use crate::cluster::{ClusterOptions, cluster};
    use crate::reader::signals::Side;
    use crate::typer::consensus::{Family, reverse_complement};

    fn alu_typing() -> crate::typer::Typing {
        crate::typer::Typing {
            family: Family::Alu,
            subfamily: "AluY".to_owned(),
            reverse: false,
            element_start: 0,
            element_end: 50,
            identity: 1.0,
            typed_clips: 2,
        }
    }

    fn unchecked_options() -> ClusterOptions {
        ClusterOptions {
            min_clips: 0,
            ..ClusterOptions::default()
        }
    }

    #[test]
    fn sampling() {
        let signals: Vec<_> = (0..100)
            .map(|i| {
                mate(
                    1_000 + i,
                    Side::Right,
                    &format!("m{i:03}"),
                    SignalKind::MateUnmapped,
                )
            })
            .collect();
        let typed = |mut c: Candidate| {
            c.typing = Some(alu_typing());
            c
        };
        let mut candidate = typed(cluster(signals.clone(), &unchecked_options()).remove(0));
        sample_mates(&mut candidate);
        let kept: Vec<usize> = candidate.mate_checks.iter().map(|m| m.signal).collect();
        assert_eq!(kept.len(), MAX_CHECKED_MATES);
        assert_eq!(kept[..4], [0, 3, 5, 8]);
        assert_eq!(kept[MAX_CHECKED_MATES - 1], 98);
        let requests = unmapped_requests(std::slice::from_ref(&candidate));
        assert_eq!(requests.len(), 1);
        assert_eq!(requests[0].1.len(), MAX_CHECKED_MATES);
        assert_eq!(requests[0].1[1].1.anchor_start, candidate.signals[3].start);

        // Few mates: all kept. Not checked: `no_junction`, or untyped.
        let mut few = typed(cluster(signals[..10].to_vec(), &unchecked_options()).remove(0));
        sample_mates(&mut few);
        assert_eq!(few.mate_checks.len(), 10);
        let mut no_junction =
            typed(cluster(signals[..10].to_vec(), &ClusterOptions::default()).remove(0));
        sample_mates(&mut no_junction);
        assert_eq!(no_junction.mate_checks, []);
        let mut untyped = cluster(signals[..10].to_vec(), &unchecked_options()).remove(0);
        sample_mates(&mut untyped);
        assert_eq!(untyped.mate_checks, []);
    }

    #[test]
    fn mei_mates_are_counted() {
        let set = ConsensusSet::load(None).unwrap();
        let alu = set.sequences()[0].sequence.clone();
        let sva = set
            .sequences()
            .iter()
            .find(|c| c.family == Family::Sva)
            .unwrap()
            .sequence
            .clone();
        let signals: Vec<_> = ["a", "b", "c", "d", "e"]
            .iter()
            .map(|read| mate(1_000, Side::Right, read, SignalKind::MateUnmapped))
            .collect();
        let mut candidate = cluster(signals, &unchecked_options()).remove(0);
        let sequences: [Option<Vec<u8>>; 5] = [
            Some(alu[20..96].to_vec()),
            Some(reverse_complement(&alu[150..226])),
            Some(b"TGCATGCTAGCTAGGATCGATCGTAGCTAGCTTTAGGCGATCGATTGCATGCTAGCTAGG".to_vec()),
            Some(sva[300..376].to_vec()),
            None,
        ];
        candidate.mate_checks = sequences
            .into_iter()
            .enumerate()
            .map(|(signal, sequence)| MateCheck {
                signal,
                sequence,
                mei: None,
            })
            .collect();

        // Untyped: any family counts.
        type_mates(&set, &mut candidate);
        assert_eq!((candidate.checked_mates, candidate.mei_mates), (4, 3));

        // Typed Alu: the SVA mate does not count.
        candidate.typing = Some(alu_typing());
        type_mates(&set, &mut candidate);
        assert_eq!((candidate.checked_mates, candidate.mei_mates), (4, 2));
        let flags: Vec<_> = candidate.mate_checks.iter().map(|m| m.mei).collect();
        assert_eq!(
            flags,
            [Some(true), Some(true), Some(false), Some(false), None]
        );
    }
}
