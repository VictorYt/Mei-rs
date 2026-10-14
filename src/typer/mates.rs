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
//!   a reference, mapped mates are not checked). Only the aligned part of the
//!   mate stands for it: mates near the 3' end of the element are often
//!   clipped (`34M42S`), their clipped part being the poly(A) tail or the
//!   insertion flank, not the reference after the copy. The reference is
//!   read over the span of the mate alignment, from its CIGAR (`MC` tag),
//!   else over a read length (aligners often leave `MC` out of pairs on two
//!   contigs).
//!
//! Mapped mates found in the reference next to their candidate
//! ([`LOCAL_WINDOW`]) are not MEI mates: the pair fits the reference without
//! an insertion. Next to a reference copy of an element, the mates of the
//! reads in that copy map with a low MAPQ to other copies, and the pairs look
//! like mates in an inserted element.
//!
//! Mates are typed without their runs of A or T ([`MIN_MATE_RUN`]): the
//! poly(A) tail does not align to the consensus sequences, and A-rich
//! sequences would align to the A-rich parts of the consensus.
//!
//! A mate in the poly(A) tail of the element has too few element bases to be
//! typed: an untyped mate with a run of [`TAIL_MATE_RUN`] A or T is not
//! counted, as an MEI mate or as another (#71).
//!
//! A reference copy of the family of the candidate near the site
//! ([`flag_reference_copies`]) makes its mate evidence unreliable (#71).
//!
//! A mate is an MEI mate when its best alignment ([`super::type_sequence`])
//! is in the family of the candidate, or in any family for an untyped one.
//!
//! Untyped candidates with a poly(A/T) tail clip have their mates checked
//! too: on short reads, the clips at the tail junction are often the tail
//! alone, with too few element bases to be typed. At least
//! [`MIN_TYPED_MATES`] mates of one family type such a candidate
//! ([`Typing::by_mates`]).
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

use std::collections::{HashMap, HashSet};

use super::consensus::{ConsensusSet, Family};
use super::{ClipHit, Typing, type_sequence};
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

/// Minimum number of mates of one family to type a candidate whose clips are
/// not typed.
pub const MIN_TYPED_MATES: usize = 2;

/// Minimum number of mate signals of a candidate without junction clips to
/// have its mates checked.
pub const MATE_ONLY_MIN_MATES: usize = 3;

/// Minimum supporting reads per read of local depth of a candidate without
/// junction clips to have its mates checked: most of them are background
/// (HGSVC3 exomes: 59 of the 61 true ones above, a seventh of the others).
pub const MATE_ONLY_MIN_RATIO: f64 = 0.2;

/// A mate sampled to be checked.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MateCheck {
    /// Index of the mate signal in [`Candidate::signals`].
    pub signal: usize,
    /// Sequence of the mate, once found.
    pub sequence: Option<Vec<u8>>,
    /// Whether it is a mobile element of the candidate family, once checked.
    pub mei: Option<bool>,
    /// A mapped mate found in the reference next to the candidate: the pair
    /// fits the reference, the mate is not an MEI mate.
    pub local: bool,
}

/// Whether the mates of a candidate are checked: it is typed and has enough
/// junction soft-clips (`no_junction` is decided at clustering), so it can
/// pass; or it is untyped and has a soft-clip with a poly(A/T) tail, so its
/// mates may type it; or, typed or not, it has too few junction clips but
/// enough mate signals for its depth ([`MATE_ONLY_MIN_MATES`],
/// [`MATE_ONLY_MIN_RATIO`]), so its mates may make it an insertion (#59,
/// #70).
#[must_use]
pub fn checked(candidate: &Candidate) -> bool {
    let no_junction = candidate.filters.contains(&Filter::NoJunction);
    let mate_only = no_junction
        && candidate.unmapped_mates + candidate.discordant_mates >= MATE_ONLY_MIN_MATES
        && candidate
            .signal_ratio()
            .is_some_and(|ratio| ratio >= MATE_ONLY_MIN_RATIO);
    match candidate.typing {
        Some(_) => !no_junction || mate_only,
        None => candidate.polya_clips + candidate.polyt_clips > 0 || mate_only,
    }
}

/// Whether a candidate can pass once its mates are checked: typed, by its
/// clips or by its mates, with enough junction soft-clips.
#[must_use]
pub fn may_pass(candidate: &Candidate) -> bool {
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
            local: false,
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
/// sequence: over the span of the mate alignment (`MC` tag), else a read
/// length from the mate start (the length of the anchored alignment); then
/// mark the local ones ([`MateCheck::local`]).
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
            let Some(length) = contigs.get(locus.contig()).map(|c| c.length) else {
                continue;
            };
            let read_length = locus.span.map_or_else(
                || (signal.end - signal.start).max(MIN_MATE_LENGTH),
                |span| u64::from(span.get()),
            );
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
        let region = Region::new(contigs[locus.contig()].name.as_str(), first..=last);
        let record = reader.query(&region)?;
        let mut sequence = record.sequence().as_ref().to_vec();
        sequence.make_ascii_uppercase();
        let check = &mut candidates[c].mate_checks[k];
        check.sequence = Some(sequence);
    }
    mark_local_mates(candidates, contigs, &mut reader)?;
    Ok(())
}

/// Bases of reference read on each side of the anchored reads of a
/// candidate, where a mate of a pair without insertion lies (fragment
/// length).
pub const LOCAL_WINDOW: u64 = 600;

/// Minimum length of the runs of A or T removed from a mate before typing.
pub const MIN_MATE_RUN: usize = 10;

/// The sequence of a mate without its runs of at least [`MIN_MATE_RUN`] A
/// or T: the poly(A) tail of the element, or that of the reference copy a
/// mate maps to, does not align to the consensus and would lower its
/// coverage.
#[must_use]
pub fn without_runs(sequence: &[u8]) -> Vec<u8> {
    sequence
        .chunk_by(|a, b| a == b)
        .filter(|run| !(matches!(run[0], b'A' | b'T') && run.len() >= MIN_MATE_RUN))
        .flatten()
        .copied()
        .collect()
}

/// Run of A or T, in bases, that makes an untyped mate a mate in the poly(A)
/// tail of the element: with too few element bases to be typed, it is not
/// counted (#71).
pub const TAIL_MATE_RUN: usize = 20;

/// Longest run of A or of T in a sequence.
fn longest_run(sequence: &[u8]) -> usize {
    sequence
        .chunk_by(|a, b| a == b)
        .filter(|run| matches!(run[0], b'A' | b'T'))
        .map(<[u8]>::len)
        .max()
        .unwrap_or(0)
}

/// Length of the k-mers compared between a mapped mate and the reference
/// next to the candidate.
pub const LOCAL_K: usize = 20;

/// Fraction of the k-mers of a mapped mate found next to the candidate for
/// it to be local.
pub const LOCAL_MIN_SHARED: f64 = 0.5;

/// Mark the mapped mates found in the reference next to their candidate
/// (within [`LOCAL_WINDOW`]) as local ([`MateCheck::local`]): the pair fits
/// the reference without an insertion. Next to a reference copy of an
/// element, the mates of the reads in that copy map with a low MAPQ to other
/// copies, and these pairs look like mates in an inserted element.
fn mark_local_mates(
    candidates: &mut [Candidate],
    contigs: &[Contig],
    reader: &mut fasta::io::IndexedReader<fasta::io::BufReader<std::fs::File>>,
) -> io::Result<()> {
    for candidate in candidates.iter_mut() {
        let mapped: Vec<usize> = candidate
            .mate_checks
            .iter()
            .enumerate()
            .filter(|(_, check)| {
                check.sequence.is_some()
                    && candidate.signals[check.signal].kind.mate_locus().is_some()
            })
            .map(|(k, _)| k)
            .collect();
        let Some(first) = mapped.first() else {
            continue;
        };
        let contig = candidate.signals[candidate.mate_checks[*first].signal].contig;
        let Some(contig) = contigs.get(contig) else {
            continue;
        };
        let start = candidate.signals.iter().map(|s| s.start).min().unwrap_or(0);
        let end = candidate.signals.iter().map(|s| s.end).max().unwrap_or(0);
        let (Some(first), Some(last)) = (
            usize::try_from(start.saturating_sub(LOCAL_WINDOW) + 1)
                .ok()
                .and_then(|p| Position::try_from(p).ok()),
            usize::try_from((end + LOCAL_WINDOW).min(contig.length))
                .ok()
                .and_then(|p| Position::try_from(p).ok()),
        ) else {
            continue;
        };
        let record = reader.query(&Region::new(contig.name.as_str(), first..=last))?;
        let window: HashSet<u64> = canonical_kmers(record.sequence().as_ref()).collect();
        for k in mapped {
            let check = &mut candidate.mate_checks[k];
            let Some(sequence) = check.sequence.as_deref() else {
                continue;
            };
            let (shared, total) = canonical_kmers(sequence)
                .fold((0usize, 0usize), |(s, t), kmer| {
                    (s + usize::from(window.contains(&kmer)), t + 1)
                });
            #[allow(clippy::cast_precision_loss)]
            let local = total > 0 && shared as f64 >= LOCAL_MIN_SHARED * total as f64;
            check.local = local;
        }
    }
    Ok(())
}

/// The [`LOCAL_K`]-mers of a sequence, each as the smaller code of its two
/// strands; k-mers with other bases than ACGT are skipped.
fn canonical_kmers(sequence: &[u8]) -> impl Iterator<Item = u64> + '_ {
    sequence.windows(LOCAL_K).filter_map(|kmer| {
        let (mut forward, mut reverse) = (0u64, 0u64);
        for (i, &b) in kmer.iter().enumerate() {
            let bits = match b.to_ascii_uppercase() {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => return None,
            };
            forward = (forward << 2) | bits;
            reverse |= (3 - bits) << (2 * i);
        }
        Some(forward.min(reverse))
    })
}

/// Bases of reference read on each side of a candidate to look for a copy of
/// its family (one fragment length; HGSVC3 exomes: 300 bp gave 128 true and
/// 42 false PASS, 600 bp 126 and 40, the fraction of true ones with a copy
/// near being small).
pub const COPY_WINDOW: u64 = 600;

/// Length of the pieces of the reference typed to find a copy, and their
/// step.
const COPY_PIECE: usize = 100;
const COPY_STEP: usize = 50;

/// Flag the candidates with a reference copy of their family near the site
/// ([`Candidate::reference_copy`]): the mates of the reads in that copy map
/// with a low MAPQ to other copies, and those not caught by the local check
/// ([`MateCheck::local`], diverged copies) look like mates in an inserted
/// element. Only candidates typed and placed by mate evidence, or with mates
/// in the poly(A) tail, are looked at: those typed by their mates, those
/// without enough junction clips (`no_junction`) and those with
/// [`Candidate::tail_mates`]; the mates in the tail are then discounted
/// ([`discount_tail_mates`]).
///
/// # Errors
///
/// Fails if the reference cannot be opened or read.
pub fn flag_reference_copies(
    candidates: &mut [Candidate],
    contigs: &[Contig],
    reference: &Path,
    set: &ConsensusSet,
) -> io::Result<()> {
    let mut reader = fasta::io::indexed_reader::Builder::default().build_from_path(reference)?;
    for candidate in candidates.iter_mut() {
        let Some(typing) = candidate.typing.as_ref() else {
            continue;
        };
        if !(typing.by_mates
            || candidate.filters.contains(&Filter::NoJunction)
            || candidate.tail_mates > 0)
        {
            continue;
        }
        let Some(contig) = contigs.get(candidate.contig) else {
            continue;
        };
        let family = typing.family;
        let start = candidate.breakpoint.saturating_sub(COPY_WINDOW);
        let end = (candidate.breakpoint + COPY_WINDOW).min(contig.length);
        let (Some(first), Some(last)) = (
            usize::try_from(start + 1)
                .ok()
                .and_then(|p| Position::try_from(p).ok()),
            usize::try_from(end)
                .ok()
                .and_then(|p| Position::try_from(p).ok()),
        ) else {
            continue;
        };
        if first > last {
            continue;
        }
        let record = reader.query(&Region::new(contig.name.as_str(), first..=last))?;
        let mut window = record.sequence().as_ref().to_vec();
        window.make_ascii_uppercase();
        candidate.reference_copy = (0..window.len().saturating_sub(COPY_PIECE - 1))
            .step_by(COPY_STEP)
            .any(|at| {
                type_sequence(set, &window[at..at + COPY_PIECE])
                    .iter()
                    .any(|hit| set.sequences()[hit.consensus].family == family)
            });
    }
    candidates.iter_mut().for_each(discount_tail_mates);
    Ok(())
}

/// Leave the mates in the poly(A) tail of the element out of the checked
/// mates ([`Candidate::tail_mates`], #71), unless a reference copy of the
/// family lies near the site: the mates of the reads in the poly(A) tail of
/// that copy look like them.
pub fn discount_tail_mates(candidate: &mut Candidate) {
    if !candidate.reference_copy {
        candidate.checked_mates -= candidate.tail_mates;
        candidate.tail_mates = 0;
    }
}

/// Type the sampled mates of a candidate with a sequence: `checked_mates` and
/// `mei_mates` (typed to its family, or to any family if it is untyped). An
/// untyped candidate is first typed by its mates (see the module documentation).
pub fn type_mates(set: &ConsensusSet, candidate: &mut Candidate) {
    let hits: Vec<Option<Option<ClipHit>>> = candidate
        .mate_checks
        .iter()
        .map(|check| {
            check.sequence.as_deref().map(|sequence| {
                if check.local {
                    return None;
                }
                type_sequence(set, &without_runs(sequence)).first().copied()
            })
        })
        .collect();
    if candidate.typing.is_none() {
        let typed: Vec<ClipHit> = hits.iter().flatten().flatten().copied().collect();
        candidate.typing = mate_typing(set, candidate, &typed);
    }
    let family = candidate.typing.as_ref().map(|t| t.family);
    let (mut checked, mut mei, mut identity, mut tail) = (0, 0, 0.0, 0);
    for (check, hit) in candidate.mate_checks.iter_mut().zip(hits) {
        let Some(hit) = hit else {
            continue;
        };
        // A mate in the tail of the element: no family to type; counted as a
        // non-MEI mate until [`discount_tail_mates`].
        if hit.is_none()
            && check
                .sequence
                .as_deref()
                .is_some_and(|sequence| longest_run(sequence) >= TAIL_MATE_RUN)
        {
            tail += 1;
        }
        checked += 1;
        let mei_hit =
            hit.filter(|hit| family.is_none_or(|f| set.sequences()[hit.consensus].family == f));
        check.mei = Some(mei_hit.is_some());
        if let Some(hit) = mei_hit {
            mei += 1;
            identity += hit.identity;
        }
    }
    candidate.checked_mates = checked;
    candidate.tail_mates = tail;
    candidate.mei_mates = mei;
    #[allow(clippy::cast_precision_loss)]
    let mean = (mei > 0).then(|| identity / mei as f64);
    candidate.mate_identity = mean;
}

/// Element typed from the best hits of the mates of an untyped candidate:
/// the family of the most mates (at least [`MIN_TYPED_MATES`]; ties go to
/// the higher summed score). Mapped mates stand for a copy of the element in
/// the reference, often of an older subfamily: the subfamily is the
/// consensus with more than half of the summed score of the family, else
/// only the family is known (its name). The element strand comes from the
/// tails (poly(T): minus strand).
fn mate_typing(set: &ConsensusSet, candidate: &Candidate, hits: &[ClipHit]) -> Option<Typing> {
    let family_of = |hit: &ClipHit| set.sequences()[hit.consensus].family;
    let mut families: HashMap<Family, (usize, i64)> = HashMap::new();
    for hit in hits {
        let entry = families.entry(family_of(hit)).or_default();
        entry.0 += 1;
        entry.1 += i64::from(hit.score);
    }
    let (&family, &(mates, family_score)) = families
        .iter()
        .max_by(|a, b| a.1.cmp(b.1).then(b.0.cmp(a.0)))?;
    if mates < MIN_TYPED_MATES {
        return None;
    }
    let in_family: Vec<&ClipHit> = hits.iter().filter(|h| family_of(h) == family).collect();
    let mut scores: HashMap<usize, i64> = HashMap::new();
    for hit in &in_family {
        *scores.entry(hit.consensus).or_default() += i64::from(hit.score);
    }
    let (&best, &best_score) = scores
        .iter()
        .max_by(|a, b| a.1.cmp(b.1).then(b.0.cmp(a.0)))?;
    let on_best: Vec<&&ClipHit> = in_family.iter().filter(|h| h.consensus == best).collect();
    #[allow(clippy::cast_precision_loss)]
    let identity = on_best.iter().map(|h| h.identity).sum::<f64>() / on_best.len() as f64;
    Some(Typing {
        family,
        subfamily: if 2 * best_score > family_score {
            set.sequences()[best].name.clone()
        } else {
            family.name().to_owned()
        },
        reverse: candidate.polyt_clips > candidate.polya_clips,
        element_start: on_best.iter().map(|h| h.element_start).min()?,
        element_end: on_best.iter().map(|h| h.element_end).max()?,
        identity,
        typed_clips: 0,
        short_typed_clips: 0,
        by_mates: true,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cluster::tests::{clip, mate};
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
            short_typed_clips: 0,
            by_mates: false,
        }
    }

    fn unchecked_options() -> ClusterOptions {
        ClusterOptions {
            min_clips: 0,
            ..ClusterOptions::default()
        }
    }

    #[test]
    fn local_mates_are_marked() {
        // Two mapped mates of a candidate at chr1:5000 of the fixture
        // reference: one maps next to it (chr1:5300, the pair fits the
        // reference), one on chr2.
        let reference =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/reference.fa");
        let contigs = [
            Contig {
                name: "chr1".to_owned(),
                length: 20_000,
            },
            Contig {
                name: "chr2".to_owned(),
                length: 10_000,
            },
        ];
        let locus = |contig, start| Locus {
            contig,
            span: std::num::NonZeroU32::new(60),
            start,
        };
        let signals = vec![
            mate(
                5_000,
                Side::Right,
                "near",
                SignalKind::MateLowMapq {
                    mapq: 0,
                    mate: Some(locus(0, 5_300)),
                },
            ),
            mate(
                5_010,
                Side::Right,
                "far",
                SignalKind::MateElsewhere {
                    mate: locus(1, 3_000),
                },
            ),
        ];
        let mut candidates = cluster(signals, &unchecked_options());
        candidates[0].mate_checks = (0..2)
            .map(|signal| MateCheck {
                signal,
                sequence: None,
                mei: None,
                local: false,
            })
            .collect();
        fetch_mapped_mates(&mut candidates, &contigs, &reference).unwrap();
        let checks = &candidates[0].mate_checks;
        let near = checks
            .iter()
            .find(|c| candidates[0].signals[c.signal].read == "near");
        let far = checks
            .iter()
            .find(|c| candidates[0].signals[c.signal].read == "far");
        assert!(near.unwrap().local);
        assert!(!far.unwrap().local);
        assert_eq!(far.unwrap().sequence.as_ref().map(Vec::len), Some(60));
    }

    #[test]
    fn local_mates_are_not_mei_mates() {
        let set = ConsensusSet::load(None).unwrap();
        let alu = set.sequences()[0].sequence.clone();
        let signals: Vec<_> = ["a", "b", "c"]
            .iter()
            .map(|read| mate(1_000, Side::Right, read, SignalKind::MateUnmapped))
            .collect();
        let mut candidate = cluster(signals, &unchecked_options()).remove(0);
        candidate.typing = Some(alu_typing());
        candidate.mate_checks = (0..3)
            .map(|signal| MateCheck {
                signal,
                sequence: Some(alu[20..96].to_vec()),
                mei: None,
                local: signal == 2,
            })
            .collect();
        type_mates(&set, &mut candidate);
        assert_eq!((candidate.checked_mates, candidate.mei_mates), (3, 2));
        assert_eq!(candidate.mate_checks[2].mei, Some(false));
    }

    #[test]
    fn tail_mates_are_discounted_without_a_reference_copy() {
        let set = ConsensusSet::load(None).unwrap();
        let alu = set.sequences()[0].sequence.clone();
        let signals: Vec<_> = ["a", "b", "c", "d"]
            .iter()
            .map(|read| mate(1_000, Side::Right, read, SignalKind::MateUnmapped))
            .collect();
        let mut candidate = cluster(signals, &unchecked_options()).remove(0);
        candidate.typing = Some(alu_typing());
        // Three MEI mates and one mate in the poly(A) tail (untyped).
        candidate.mate_checks = (0..4)
            .map(|signal| MateCheck {
                signal,
                sequence: Some(if signal < 3 {
                    alu[20..96].to_vec()
                } else {
                    vec![b'A'; 60]
                }),
                mei: None,
                local: false,
            })
            .collect();
        type_mates(&set, &mut candidate);
        assert_eq!(
            (
                candidate.checked_mates,
                candidate.mei_mates,
                candidate.tail_mates
            ),
            (4, 3, 1)
        );
        // A reference copy of the family nearby: the tail mate stays counted.
        let mut near_copy = candidate.clone();
        near_copy.reference_copy = true;
        discount_tail_mates(&mut near_copy);
        assert_eq!(near_copy.checked_mates, 4);
        discount_tail_mates(&mut candidate);
        assert_eq!((candidate.checked_mates, candidate.tail_mates), (3, 0));
    }

    #[test]
    fn reference_copies_are_flagged() {
        // chr2:7001-7281 of the fixture reference is a copy of the Alu.
        let reference =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/reference.fa");
        let contigs = [
            Contig {
                name: "chr1".to_owned(),
                length: 20_000,
            },
            Contig {
                name: "chr2".to_owned(),
                length: 10_000,
            },
        ];
        let set = ConsensusSet::load(None).unwrap();
        let candidate_at = |contig: usize, position: u64| {
            let mut signal = mate(position, Side::Right, "m", SignalKind::MateUnmapped);
            signal.contig = contig;
            let mut candidate = cluster(vec![signal], &unchecked_options()).remove(0);
            candidate.typing = Some(alu_typing());
            candidate.filters.push(Filter::NoJunction);
            candidate
        };
        // 100 bp before the copy, and on chr1 far from any copy.
        let mut candidates = [candidate_at(1, 6_900), candidate_at(0, 5_000)];
        flag_reference_copies(&mut candidates, &contigs, &reference, &set).unwrap();
        assert!(candidates[0].reference_copy);
        assert!(!candidates[1].reference_copy);
    }

    #[test]
    fn mates_without_their_tail() {
        // A mate in the tail of the element: its A run would halve the
        // coverage of the element end.
        let set = ConsensusSet::load(None).unwrap();
        let alu = &set.sequences()[0].sequence;
        let mut sequence = vec![b'A'; 40];
        sequence.extend_from_slice(&alu[200..260]);
        assert_eq!(without_runs(&sequence).len(), sequence.len() - 40);
        assert_ne!(type_sequence(&set, &without_runs(&sequence)), []);
        assert_eq!(without_runs(b"ACGTAAAAAAAAACGT"), b"ACGTAAAAAAAAACGT");
    }

    #[test]
    fn canonical_kmers_of_both_strands() {
        let sequence = b"ACGTTGCATGCATTGACCAGTA";
        let reverse = reverse_complement(sequence);
        let forward: Vec<u64> = canonical_kmers(sequence).collect();
        let mut backward: Vec<u64> = canonical_kmers(&reverse).collect();
        backward.reverse();
        assert_eq!(forward.len(), sequence.len() - LOCAL_K + 1);
        assert_eq!(forward, backward);
        assert_eq!(canonical_kmers(b"ACGTNACGTACGTACGTACGTAC").count(), 0);
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

        // Few mates: all kept. Not checked: `no_junction`, or untyped
        // without a tail.
        let mut few = typed(cluster(signals[..10].to_vec(), &unchecked_options()).remove(0));
        sample_mates(&mut few);
        assert_eq!(few.mate_checks.len(), 10);
        let mut no_junction =
            typed(cluster(signals[..10].to_vec(), &ClusterOptions::default()).remove(0));
        sample_mates(&mut no_junction);
        assert_eq!(no_junction.mate_checks, []);
        // Unless its mates are many for its depth: they may make it an
        // insertion (#70).
        no_junction.depth = Some(20);
        sample_mates(&mut no_junction);
        assert_eq!(no_junction.mate_checks.len(), 10);
        let mut untyped = cluster(signals[..10].to_vec(), &unchecked_options()).remove(0);
        sample_mates(&mut untyped);
        assert_eq!(untyped.mate_checks, []);
        // Untyped with a tail clip, even without enough junction clips: its
        // mates may type it.
        let mut tailed = signals[..10].to_vec();
        tailed.push(clip(1_005, Side::Right, "t", Some(b'A')));
        tailed.sort();
        let mut tailed = cluster(tailed, &ClusterOptions::default()).remove(0);
        assert!(tailed.filters.contains(&Filter::NoJunction));
        sample_mates(&mut tailed);
        assert_eq!(tailed.mate_checks.len(), 10);
        assert!(!may_pass(&tailed));
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
        let random = b"TGCATGCTAGCTAGGATCGATCGTAGCTAGCTTTAGGCGATCGATTGCATGCTAGCTAGG".to_vec();
        let checks = |sequences: Vec<Option<Vec<u8>>>| -> Vec<MateCheck> {
            sequences
                .into_iter()
                .enumerate()
                .map(|(signal, sequence)| MateCheck {
                    signal,
                    sequence,
                    mei: None,
                    local: false,
                })
                .collect()
        };

        // Untyped, one mate per family: not typed, any family counts.
        candidate.mate_checks = checks(vec![
            Some(alu[20..96].to_vec()),
            Some(random.clone()),
            Some(sva[300..376].to_vec()),
            None,
        ]);
        type_mates(&set, &mut candidate);
        assert_eq!(candidate.typing, None);
        assert_eq!((candidate.checked_mates, candidate.mei_mates), (3, 2));

        // Untyped, two Alu mates: typed Alu by its mates; the SVA mate does
        // not count.
        candidate.mate_checks = checks(vec![
            Some(alu[20..96].to_vec()),
            Some(reverse_complement(&alu[150..226])),
            Some(random),
            Some(sva[300..376].to_vec()),
            None,
        ]);
        type_mates(&set, &mut candidate);
        let typing = candidate.typing.clone().unwrap();
        assert!(typing.by_mates);
        assert_eq!(
            (typing.family, typing.subfamily.as_str(), typing.typed_clips),
            (Family::Alu, "AluY", 0)
        );
        assert_eq!((typing.element_start, typing.element_end), (20, 226));
        assert_eq!((candidate.checked_mates, candidate.mei_mates), (4, 2));

        // Typed Alu by its clips: the SVA mate does not count.
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
