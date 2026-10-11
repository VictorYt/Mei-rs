//! Length variants of reference tracts (homopolymers and short tandem
//! repeats) mistaken for insertions.
//!
//! When the reads carry a longer A tract (or `(CAAAA)n`, `(AAATAAT)n`...) than
//! the reference, typically at the 3' end of an Alu of the reference, the
//! aligner soft-clips them on both sides of the tract: the clips have a
//! poly(A/T) tail and the overlap of the two junctions looks like a TSD. Such
//! a clip is the reference itself with more repeat units: beyond the tract,
//! its far end is the reference sequence that flanks the tract.
//!
//! A candidate is flagged [`Filter::TandemRepeat`] when
//! - at least [`MIN_EXPLAINED`] of its junction soft-clips, and at least half
//!   of the informative ones, end in the reference where the clip would end
//!   if the read were the reference with a longer tract (the far end of the
//!   clip, homopolymer-compressed, found within [`ANCHOR_SLACK`] bp of a clip
//!   length away from the junction, on the clip side); or
//! - its tailed clips hold nothing but the tail and the reference has a
//!   homopolymer of the tail base of at least [`MIN_HOMOPOLYMER`] bp at their
//!   junction (expansion of a reference homopolymer).
//!
//! A real insertion next to a reference Alu can be flagged too: its clips are
//! Alu sequence, found in the neighbouring Alu (2 of 60 true candidates on
//! the HGSVC3 exomes of #26).
//!
//! Coordinates: [`ReferenceWindow::start`] is the 0-based position of the
//! first base of the window, signal positions are boundaries (see
//! [`crate::reader::signals`]).

use super::{Candidate, Filter, JUNCTION_TOLERANCE};
use crate::reader::signals::{Side, SignalKind};

/// Reference bases fetched on each side of the breakpoint: junctions up to
/// 50 bp apart, clips up to a read length, plus the slack.
pub const WINDOW: u64 = 250;

/// Length of the homopolymer-compressed far end of a clip that must be found
/// in the reference.
pub const ANCHOR_LENGTH: usize = 8;

/// Distinct bases required in the far end of a clip for it to be
/// informative (low-complexity ends such as `ATATATAT` match anywhere).
pub const MIN_ANCHOR_BASES: usize = 3;

/// Maximum distance, in bp, between where the far end of a clip is found in
/// the reference and where it would be if the read were the reference.
pub const ANCHOR_SLACK: i64 = 30;

/// Compressed bases that may be skipped at the read end (sequencing errors);
/// each one skipped widens the slack by 3 bp.
pub const READ_END_TRIM: usize = 2;

/// Minimum number of junction soft-clips explained by the reference.
pub const MIN_EXPLAINED: usize = 2;

/// Minimum length of a reference homopolymer at the junction of tail-only
/// clips.
pub const MIN_HOMOPOLYMER: usize = 10;

/// Maximum distance, in bp, between that homopolymer and the junction.
pub const HOMOPOLYMER_REACH: usize = 3;

/// Reference bases around a candidate.
#[derive(Clone, Copy, Debug)]
pub struct ReferenceWindow<'a> {
    /// 0-based position of the first base.
    pub start: u64,
    /// Bases, upper case.
    pub sequence: &'a [u8],
}

/// Reference interval to fetch for a candidate: 0-based, half-open, clipped
/// to the contig.
#[must_use]
pub fn window(candidate: &Candidate, contig_length: u64) -> (u64, u64) {
    (
        candidate.breakpoint.saturating_sub(WINDOW),
        (candidate.breakpoint + WINDOW).min(contig_length),
    )
}

/// Flag the candidate [`Filter::TandemRepeat`] if its soft-clips are a length
/// variant of a reference tract.
pub fn apply(candidate: &mut Candidate, reference: &ReferenceWindow<'_>) {
    if explained_by_reference(candidate, reference) {
        candidate.filters.push(Filter::TandemRepeat);
    }
}

/// A run of identical bases: base and 0-based half-open interval in the
/// sequence it comes from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Run {
    base: u8,
    start: usize,
    end: usize,
}

/// Homopolymer runs of a sequence.
fn runs(sequence: &[u8]) -> Vec<Run> {
    let mut runs: Vec<Run> = Vec::new();
    for (i, &base) in sequence.iter().enumerate() {
        match runs.last_mut() {
            Some(run) if run.base == base => run.end = i + 1,
            _ => runs.push(Run {
                base,
                start: i,
                end: i + 1,
            }),
        }
    }
    runs
}

/// Whether the soft-clips of a candidate are explained by the reference.
#[must_use]
pub fn explained_by_reference(candidate: &Candidate, reference: &ReferenceWindow<'_>) -> bool {
    let reference_runs = runs(reference.sequence);
    let near = |junction: Option<u64>, position: u64| {
        junction.is_some_and(|j| j.abs_diff(position) <= JUNCTION_TOLERANCE)
    };
    let mut informative = 0;
    let mut explained = 0;
    let mut tailed = Vec::new();
    let mut tailed_informative = false;
    for signal in &candidate.signals {
        let SignalKind::SoftClip { sequence, tail, .. } = &signal.kind else {
            continue;
        };
        // Clips on the right end of reads make the left junction.
        let junction = match signal.side {
            Side::Right => candidate.left_junction,
            Side::Left => candidate.right_junction,
        };
        if !near(junction, signal.position) {
            continue;
        }
        let boundary = offset(signal.position, reference.start);
        let clip_length = i64::try_from(sequence.len()).unwrap_or(i64::MAX);
        let expected = match signal.side {
            Side::Left => boundary - clip_length,
            Side::Right => boundary + clip_length,
        };
        // The tail lies at the junction end of the clip. Only its pure run is
        // set aside to anchor the clip: a tail of `(CAAAA)n` passes the tail
        // test, but it is the repeat whose far end is informative.
        let pure = tail.map_or(0, |t| {
            junction_end(sequence, signal.side)
                .take_while(|&&b| b == t.base)
                .count()
        });
        if let Some(found) = anchored(
            without_junction_end(sequence, signal.side, pure),
            signal.side,
            expected,
            &reference_runs,
        ) {
            informative += 1;
            explained += usize::from(found);
        }
        if let Some(t) = tail {
            tailed.push((boundary, t.base));
            // Whether the clip holds more than the tail found by the tail
            // test (homopolymer rule).
            let rest = without_junction_end(
                sequence,
                signal.side,
                (t.length as usize).min(sequence.len()),
            );
            tailed_informative |= anchored(rest, signal.side, expected, &reference_runs).is_some();
        }
    }
    let anchored_clips = explained >= MIN_EXPLAINED && 2 * explained >= informative;
    let homopolymer = !tailed.is_empty()
        && !tailed_informative
        && tailed
            .iter()
            .any(|&(boundary, base)| homopolymer_at(&reference_runs, boundary, base));
    anchored_clips || homopolymer
}

/// Bases of a clip from its junction end: the end of a clip before the read,
/// the start of one after it.
fn junction_end(sequence: &[u8], side: Side) -> Box<dyn Iterator<Item = &u8> + '_> {
    match side {
        Side::Left => Box::new(sequence.iter().rev()),
        Side::Right => Box::new(sequence.iter()),
    }
}

/// A clip without `length` bases at its junction end.
fn without_junction_end(sequence: &[u8], side: Side, length: usize) -> &[u8] {
    match side {
        Side::Left => &sequence[..sequence.len() - length],
        Side::Right => &sequence[length..],
    }
}

/// Signed offset of a boundary in the window.
fn offset(position: u64, start: u64) -> i64 {
    i64::try_from(position).unwrap_or(i64::MAX) - i64::try_from(start).unwrap_or(i64::MAX)
}

/// Whether the far end of a clip (without its tail) is found in the
/// reference where expected: `None` if the clip has no informative far end.
fn anchored(rest: &[u8], side: Side, expected: i64, reference: &[Run]) -> Option<bool> {
    let compressed: Vec<u8> = runs(rest).iter().map(|r| r.base).collect();
    let mut informative = false;
    for trim in 0..=READ_END_TRIM {
        // The far end: start of a clip before the read, end of one after it.
        let Some(anchor) = (match side {
            Side::Left => compressed.get(trim..trim + ANCHOR_LENGTH),
            Side::Right => compressed
                .len()
                .checked_sub(trim + ANCHOR_LENGTH)
                .map(|from| &compressed[from..compressed.len() - trim]),
        }) else {
            break;
        };
        let mut bases = anchor.to_vec();
        bases.sort_unstable();
        bases.dedup();
        if bases.len() < MIN_ANCHOR_BASES {
            continue;
        }
        informative = true;
        let slack = ANCHOR_SLACK + 3 * i64::try_from(trim).unwrap_or(0);
        let found = reference
            .windows(ANCHOR_LENGTH)
            .filter(|w| w.iter().map(|r| r.base).eq(anchor.iter().copied()))
            .any(|w| {
                let at = match side {
                    Side::Left => w[0].start,
                    Side::Right => w[ANCHOR_LENGTH - 1].end,
                };
                (i64::try_from(at).unwrap_or(i64::MAX) - expected).abs() <= slack
            });
        if found {
            return Some(true);
        }
    }
    informative.then_some(false)
}

/// Whether a homopolymer of `base` of at least [`MIN_HOMOPOLYMER`] bp lies
/// within [`HOMOPOLYMER_REACH`] bp of the boundary.
fn homopolymer_at(reference: &[Run], boundary: i64, base: u8) -> bool {
    let reach = i64::try_from(HOMOPOLYMER_REACH).unwrap_or(0);
    reference.iter().any(|r| {
        r.base == base
            && r.end - r.start >= MIN_HOMOPOLYMER
            && i64::try_from(r.start).unwrap_or(i64::MAX) <= boundary + reach
            && i64::try_from(r.end).unwrap_or(i64::MAX) >= boundary - reach
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cluster::{ClusterOptions, cluster};
    use crate::reader::signals::{PolyTail, Signal};

    /// A soft-clip with its sequence; the tail is the run of `A`/`T` at the
    /// junction end, if any.
    fn clip(position: u64, side: Side, read: &str, sequence: &str) -> Signal {
        let bytes = sequence.as_bytes();
        let end = match side {
            Side::Left => bytes.iter().rev().copied().collect::<Vec<_>>(),
            Side::Right => bytes.to_vec(),
        };
        let tail = b"AT".iter().copied().find_map(|base| {
            let length = end.iter().take_while(|&&b| b == base).count();
            (length >= 10).then(|| PolyTail {
                base,
                length: u32::try_from(length).unwrap(),
            })
        });
        let mut signal = crate::cluster::tests::clip(position, side, read, None);
        signal.kind = SignalKind::SoftClip {
            sequence: bytes.to_vec(),
            median_quality: 30,
            tail,
            supplementary: false,
            short: false,
        };
        signal
    }

    fn candidate(signals: Vec<Signal>) -> Candidate {
        let mut signals = signals;
        signals.sort();
        cluster(signals, &ClusterOptions::default()).remove(0)
    }

    /// `chr1:40089250-40089340` of `GRCh38`: 3' end of a reference Alu, A9 tract
    /// at bases 40089287-40089295 (1-based).
    const ALU_END: &[u8] =
        b"CCCTCCAGCCTGGGCAACAAGAGCGAAACTCCATCTCAAAAAAAAAGATCTCATTTGAATATCACAAGGAACTACTGCTGGAATCACTGTG";
    const ALU_END_START: u64 = 40_089_249;

    fn alu_end() -> ReferenceWindow<'static> {
        ReferenceWindow {
            start: ALU_END_START,
            sequence: ALU_END,
        }
    }

    #[test]
    fn poly_a_length_variant() {
        // Reads with a longer A tract (HG00512): clipped before the tract
        // (left junction 40089295) and after the bases before it (right
        // junction 40089286).
        let c = candidate(vec![
            clip(40_089_286, Side::Left, "l1", "TCCATCTCAAAAAAAAAAAAAAAAA"),
            clip(40_089_286, Side::Left, "l2", "ACTCCATCTCAAAAAAAAAAAA"),
            clip(40_089_295, Side::Right, "r1", "AAAAAAAAAAAAAAAAGATCTCATT"),
            clip(
                40_089_295,
                Side::Right,
                "r2",
                "AAAAAAAAAAAAAAAAAGATCTCATTTG",
            ),
        ]);
        assert_eq!(c.tsd_length(), Some(9));
        assert!(explained_by_reference(&c, &alu_end()));
    }

    #[test]
    fn repeat_tail_length_variant() {
        // NA12878, chr2:201350686: `(CAAAA)n` at the end of a reference Alu.
        // The tail test takes the whole repeat as a tail (A23, A17), so
        // only the pure run of the tail is set aside to anchor the clips.
        let window = ReferenceWindow {
            start: 201_350_639,
            sequence: b"CACTGCACTCTAGCCTGGGCAACAGAGTGAGACTCCATCTCAAAAACAAAACAAAACAAAACAATACTTACTTTGGGGAAACAACAGGTTT",
        };
        let with_tail = |mut signal: Signal, length: u32| {
            if let SignalKind::SoftClip { tail, .. } = &mut signal.kind {
                *tail = Some(PolyTail { base: b'A', length });
            }
            signal
        };
        let c = candidate(vec![
            with_tail(
                clip(
                    201_350_686,
                    Side::Left,
                    "l1",
                    "TCCATCTCAAAAACAAAACAAAACAAAACAA",
                ),
                23,
            ),
            with_tail(
                clip(
                    201_350_686,
                    Side::Left,
                    "l2",
                    "TCCATCTCAAAAACAAAACAAAACAAAACAA",
                ),
                23,
            ),
            with_tail(
                clip(201_350_686, Side::Left, "l3", "TCTCAAAAACAAAACAAAACAAAACAA"),
                23,
            ),
            with_tail(
                clip(201_350_703, Side::Right, "r1", "AAAACAAAACAAAACAATACTTACTT"),
                17,
            ),
            with_tail(
                clip(201_350_703, Side::Right, "r2", "AAAACAAAACAAAACAATACTTA"),
                17,
            ),
        ]);
        assert_eq!(c.tsd_length(), Some(17));
        assert!(explained_by_reference(&c, &window));
    }

    #[test]
    fn insertion_is_not_explained() {
        // An Alu inserted at the same place: the clips end in Alu sequence
        // that is not in the reference there.
        let c = candidate(vec![
            clip(
                40_089_286,
                Side::Left,
                "l1",
                "GGCCGGGCGCGGTGGCTCAAAAAAAAAAAAAA",
            ),
            clip(
                40_089_286,
                Side::Left,
                "l2",
                "TGGCTCACGCCTGTAAAAAAAAAAAAAAAAA",
            ),
            clip(40_089_295, Side::Right, "r1", "GGCCGGGCGCGGTGGCTCACGCC"),
            clip(
                40_089_295,
                Side::Right,
                "r2",
                "GGCCGGGCGCGGTGGCTCACGCCTGTAATC",
            ),
        ]);
        assert!(!explained_by_reference(&c, &alu_end()));
    }

    #[test]
    fn far_end_must_be_where_expected() {
        // The far end of the clips is reference sequence, but 40 bp away
        // from where a longer tract would put it.
        let c = candidate(vec![
            clip(40_089_295, Side::Right, "r1", "AAAAAAAAAAAAAGGCAACAAGAGCG"),
            clip(40_089_295, Side::Right, "r2", "AAAAAAAAAAAAAAGGCAACAAGAGCG"),
        ]);
        assert!(!explained_by_reference(&c, &alu_end()));
    }

    #[test]
    fn low_complexity_ends_are_not_informative() {
        // `TATATATA` once compressed: found anywhere, so it says nothing,
        // even where a longer tract would put it.
        let reference = b"GGCATTCGCTTATATATATATATTGCCAGTACGTGCATGCCAT";
        let window = ReferenceWindow {
            start: 1_000,
            sequence: reference,
        };
        let signals = vec![
            clip(1_030, Side::Left, "l1", "TATATATATAAAAAAAAAAAAA"),
            clip(1_030, Side::Left, "l2", "TATATATATAAAAAAAAAAAAA"),
            clip(1_030, Side::Left, "l3", "TTATATATATAAAAAAAAAAAA"),
        ];
        assert!(!explained_by_reference(&candidate(signals), &window));

        // The same position with a clip of three bases is informative.
        let signals = vec![
            clip(1_030, Side::Left, "l1", "GCTTATATATAAAAAAAAAAAAA"),
            clip(1_030, Side::Left, "l2", "CGCTTATATATAAAAAAAAAAAA"),
        ];
        assert!(explained_by_reference(&candidate(signals), &window));
    }

    #[test]
    fn reference_homopolymer() {
        // chr5:160412537: A13 in the reference, tail-only clips at its end.
        let reference = b"TGGGCTTTCAATGTTCGTTTTGCCTTTAAAAAAAAAAAAAAGAAAGAAAGAAAGAAAAAGTAAAC";
        let window = ReferenceWindow {
            start: 160_412_496,
            sequence: reference,
        };
        let signals = vec![
            clip(160_412_537, Side::Right, "r1", "AAAAAAAAAAAAAAAAAAAAAAA"),
            clip(160_412_537, Side::Right, "r2", "AAAAAAAAAAAAAAAAAAAAAAAAA"),
            clip(
                160_412_537,
                Side::Right,
                "r3",
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            ),
        ];
        assert!(explained_by_reference(&candidate(signals.clone()), &window));

        // The same clips away from any homopolymer: kept.
        let window = ReferenceWindow {
            start: 160_412_496,
            sequence: b"TGGGCTTTCAATGTTCGTTTTGCCTTTGCAGCCTGACAGAGAAAGAAAGAAAGAAAAAGTAAAC",
        };
        assert!(!explained_by_reference(&candidate(signals), &window));
    }

    #[test]
    fn flag() {
        let mut c = candidate(vec![
            clip(40_089_286, Side::Left, "l1", "TCCATCTCAAAAAAAAAAAAAAAAA"),
            clip(40_089_295, Side::Right, "r1", "AAAAAAAAAAAAAAAAGATCTCATT"),
            clip(
                40_089_295,
                Side::Right,
                "r2",
                "AAAAAAAAAAAAAAAAAGATCTCATTTG",
            ),
        ]);
        assert!(c.is_pass());
        apply(&mut c, &alu_end());
        assert_eq!(c.filters, [Filter::TandemRepeat]);
        assert_eq!(c.confidence(), None);
    }

    #[test]
    fn window_is_clipped_to_the_contig() {
        let mut c = candidate(vec![clip(100, Side::Right, "r", "AAAAAAAAAAAAAAAAAAAAA")]);
        assert_eq!(window(&c, 10_000), (0, 350));
        c.breakpoint = 9_900;
        assert_eq!(window(&c, 10_000), (9_650, 10_000));
    }

    #[test]
    fn runs_of_bases() {
        let r = runs(b"AAGCCC");
        assert_eq!(
            r,
            [
                Run {
                    base: b'A',
                    start: 0,
                    end: 2
                },
                Run {
                    base: b'G',
                    start: 2,
                    end: 3
                },
                Run {
                    base: b'C',
                    start: 3,
                    end: 6
                },
            ]
        );
    }
}
