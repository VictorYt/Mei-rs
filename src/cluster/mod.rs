//! Phase 2: spatial clustering of the signals into candidate insertion sites,
//! breakpoint estimation, capture artefact filters ([`filters`]) and length
//! variants of reference tracts ([`tracts`]).
//!
//! Positions are *boundaries* (see [`crate::reader::signals`]): a junction at
//! boundary `b` means the insertion lies between bases `b - 1` and `b`
//! (0-based), i.e. right after base `b` in 1-based coordinates.

pub mod filters;
pub mod tracts;

use std::collections::{BTreeMap, HashSet};

use crate::reader::signals::{Side, Signal, SignalKind};

/// Clustering options.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ClusterOptions {
    /// Maximum distance between two consecutive signals of a cluster.
    pub window: u64,
    /// Minimum number of distinct supporting reads.
    pub min_support: usize,
    /// Minimum number of soft-clips at the junctions; each junction needs as
    /// many to give a TSD.
    pub min_clips: usize,
    /// Require a soft-clip with a poly(A/T) tail (`no_tail`). Off by default:
    /// typed clips (`untyped`) are the evidence of a retrotransposon, and the
    /// tail sets the confidence level.
    pub require_tail: bool,
}

impl Default for ClusterOptions {
    fn default() -> Self {
        Self {
            window: 100,
            min_support: 3,
            min_clips: 2,
            require_tail: false,
        }
    }
}

/// Distance within which a soft-clip supports a junction.
pub const JUNCTION_TOLERANCE: u64 = 2;

/// Span, in bp, beyond the inner edge of a group of mate signals in which
/// the insertion can lie when only one side of it has reads: a fragment
/// length minus a read length on exomes.
pub const MATE_SPAN: u64 = 300;

/// Overlap, in bp, tolerated between the upstream and downstream groups of
/// mate signals of one insertion.
pub const MATE_OVERLAP: u64 = 20;

/// Reason a candidate is filtered out.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Filter {
    /// Fewer distinct supporting reads than `--min-support`.
    LowSupport,
    /// Fewer soft-clips at the junctions than `--min-clips` (short ones
    /// count only with an MEI mate): mate signals alone do not resolve the
    /// insertion to the base, unless reads clipped by a few bases and young
    /// MEI mates do (see [`filters::mate_only_insertion`]).
    NoJunction,
    /// No soft-clip with a poly(A/T) tail (only with `--require-tail`).
    NoTail,
    /// Fewer junction soft-clips typed to an MEI consensus than
    /// `--min-typed-clips` (short ones count only with an MEI mate), and not
    /// one typed clip with a poly(A/T) tail and at least as many MEI mates:
    /// nothing shows that the inserted sequence is a mobile element.
    Untyped,
    /// Of the mates checked against the consensus sequences, fewer than
    /// `--min-mei-mate-fraction` are mobile elements: the mate signals are
    /// mostly noise.
    NonMeiMates,
    /// The soft-clips are a length variant of a reference homopolymer or
    /// tandem repeat, not an insertion (see [`tracts`]).
    TandemRepeat,
    /// Every soft-clip stops at a probe edge, without a poly(A/T) tail or a
    /// TSD (capture artefact).
    ProbeEdge,
    /// Too few supporting reads for the local depth.
    LowRatio,
    /// Quality score below `--min-qual` (see [`crate::genotype::score`]).
    LowQual,
}

impl Filter {
    /// Name, as written in the outputs (and later in the VCF `FILTER`).
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::LowSupport => "low_support",
            Self::NoJunction => "no_junction",
            Self::NoTail => "no_tail",
            Self::Untyped => "untyped",
            Self::NonMeiMates => "non_mei_mates",
            Self::TandemRepeat => "tandem_repeat",
            Self::ProbeEdge => "probe_edge",
            Self::LowRatio => "low_ratio",
            Self::LowQual => "low_qual",
        }
    }
}

/// Confidence level of a candidate that passes the filters.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Confidence {
    /// Typed soft-clips with a poly(A/T) tail and a TSD: complete evidence
    /// of a retrotransposition.
    High,
    /// Typed soft-clips with a poly(A/T) tail, without a TSD.
    Medium,
    /// Typed soft-clips without a poly(A/T) tail (often the tailed junction
    /// is not covered by clips long enough).
    Low,
}

impl Confidence {
    /// Name, as written in the outputs.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::High => "high",
            Self::Medium => "medium",
            Self::Low => "low",
        }
    }
}

/// A candidate insertion site: a cluster of signals.
#[derive(Clone, Debug, PartialEq)]
pub struct Candidate {
    /// Contig (index in the header).
    pub contig: usize,
    /// First signal position (boundary).
    pub start: u64,
    /// Last signal position (boundary).
    pub end: u64,
    /// Estimated insertion point (boundary): the most supported soft-clip
    /// junction, or else estimated from the mate signals (see
    /// [`Self::mate_interval`]).
    pub breakpoint: u64,
    /// Without enough junction soft-clips (`--min-clips`) but with mate
    /// signals, the interval (boundaries, inclusive) where the insertion lies
    /// according to them; `None` with junction soft-clips, or once the
    /// breakpoint is found among the soft-clips of the reads.
    pub mate_interval: Option<(u64, u64)>,
    /// Reads soft-clipped by at least
    /// [`crate::reader::breakpoints::MIN_CLIP`] bases at the breakpoint,
    /// counted when it is looked for again among all clips (candidates
    /// without junction soft-clips, see [`crate::reader::breakpoints`]).
    pub breakpoint_clips: usize,
    /// End of the upstream flank: most supported junction of the clips on
    /// the right end of reads, if any; the assembled one once assembled.
    pub left_junction: Option<u64>,
    /// Start of the downstream flank: most supported junction of the clips
    /// on the left end of reads, if any; the assembled one once assembled.
    pub right_junction: Option<u64>,
    /// Soft-clips supporting the left junction (within
    /// [`JUNCTION_TOLERANCE`]).
    pub left_junction_clips: usize,
    /// Soft-clips supporting the right junction (within
    /// [`JUNCTION_TOLERANCE`]).
    pub right_junction_clips: usize,
    /// Of the soft-clips supporting either junction, short ones (see
    /// [`crate::reader::signals::ClipOptions::is_short`]).
    pub short_junction_clips: usize,
    /// Length of the target site duplication: assembled, or suggested by
    /// the junctions.
    tsd: Option<u64>,
    /// Signals of the cluster, sorted.
    pub signals: Vec<Signal>,
    /// Number of distinct supporting reads (by name).
    pub support: usize,
    /// Soft-clips on the left end of reads (insertion before the read).
    pub left_clips: usize,
    /// Soft-clips on the right end of reads (insertion after the read).
    pub right_clips: usize,
    /// Reads whose mate is unmapped.
    pub unmapped_mates: usize,
    /// Reads whose mate maps with a low MAPQ, to another contig or far away.
    pub discordant_mates: usize,
    /// Soft-clips with a poly(A) tail.
    pub polya_clips: usize,
    /// Soft-clips with a poly(T) tail.
    pub polyt_clips: usize,
    /// Number of distinct fragment start positions (strand-aware), a proxy
    /// for the number of distinct (non-duplicate) fragments.
    pub distinct_starts: usize,
    /// Local depth at the breakpoint (filtered reads), once computed.
    pub depth: Option<u64>,
    /// Element typed from the junction soft-clips, once typed.
    pub typing: Option<crate::typer::Typing>,
    /// Mates sampled to be checked against the consensus sequences (see
    /// [`crate::typer::mates`]).
    pub mate_checks: Vec<crate::typer::mates::MateCheck>,
    /// Mates whose sequence was checked.
    pub checked_mates: usize,
    /// Of these, mates typed to the family of the candidate (to any family
    /// if it is untyped).
    pub mei_mates: usize,
    /// Mean identity of the MEI mates to the consensus, once checked.
    pub mate_identity: Option<f64>,
    /// A reference copy of the family of the candidate lies near its site:
    /// its mate evidence is not used ([`crate::typer::mates::flag_reference_copies`], #71).
    pub reference_copy: bool,
    /// Checked mates in the poly(A) tail of the element (untyped, a run of A or T), counted in
    /// [`Self::checked_mates`] until [`crate::typer::mates::discount_tail_mates`] (#71).
    pub tail_mates: usize,
    /// Local assembly of the junctions, once done (see [`crate::assembly`]).
    pub assembly: Option<crate::assembly::Assembly>,
    /// Genotype, once the reads of each allele are counted at the junctions
    /// (see [`crate::genotype`]).
    pub genotype: Option<crate::genotype::Genotype>,
    /// Quality score, once the filters are applied (see
    /// [`crate::genotype::score`]).
    pub qual: Option<u32>,
    /// Filters the candidate fails; empty if it passes.
    pub filters: Vec<Filter>,
}

impl Candidate {
    /// Build a candidate from a non-empty, sorted cluster of signals.
    fn new(signals: Vec<Signal>, options: &ClusterOptions) -> Self {
        let clips = |side: Side| {
            signals
                .iter()
                .filter(move |s| s.side == side && matches!(s.kind, SignalKind::SoftClip { .. }))
        };
        let left_junction = mode(clips(Side::Right).map(|s| s.position));
        let right_junction = mode(clips(Side::Left).map(|s| s.position));
        // Clips of a junction: all of them, and the short ones.
        let near = |junction: Option<u64>, side: Side| {
            let at = clips(side)
                .filter(|s| junction.is_some_and(|j| s.position.abs_diff(j) <= JUNCTION_TOLERANCE));
            at.fold((0, 0), |(all, short), s| {
                let is_short = matches!(s.kind, SignalKind::SoftClip { short: true, .. });
                (all + 1, short + usize::from(is_short))
            })
        };
        let (left_junction_clips, left_short) = near(left_junction, Side::Right);
        let (right_junction_clips, right_short) = near(right_junction, Side::Left);
        let short_junction_clips = left_short + right_short;
        let tsd = match (left_junction, right_junction) {
            (Some(left), Some(right))
                if left_junction_clips >= options.min_clips
                    && right_junction_clips >= options.min_clips =>
            {
                left.checked_sub(right)
                    .filter(|overlap| (1..=50).contains(overlap))
            }
            _ => None,
        };
        let (breakpoint, mate_interval) = place(
            &signals,
            left_junction_clips + right_junction_clips >= options.min_clips,
        );

        let count = |f: &dyn Fn(&Signal) -> bool| signals.iter().filter(|s| f(s)).count();
        let tails = |base: u8| {
            count(
                &|s| matches!(&s.kind, SignalKind::SoftClip { tail: Some(t), .. } if t.base == base),
            )
        };
        let support = signals
            .iter()
            .map(|s| s.read.as_str())
            .collect::<HashSet<_>>()
            .len();
        let distinct_starts = signals
            .iter()
            .map(|s| (s.reverse, if s.reverse { s.end } else { s.start }))
            .collect::<HashSet<_>>()
            .len();

        let mut candidate = Self {
            contig: signals[0].contig,
            start: signals.iter().map(|s| s.position).min().unwrap_or(0),
            end: signals.iter().map(|s| s.position).max().unwrap_or(0),
            breakpoint,
            mate_interval,
            breakpoint_clips: 0,
            left_junction,
            right_junction,
            left_junction_clips,
            right_junction_clips,
            short_junction_clips,
            tsd,
            support,
            left_clips: clips(Side::Left).count(),
            right_clips: clips(Side::Right).count(),
            unmapped_mates: count(&|s| s.kind == SignalKind::MateUnmapped),
            discordant_mates: count(&|s| {
                s.kind.is_mate_signal() && s.kind != SignalKind::MateUnmapped
            }),
            polya_clips: tails(b'A'),
            polyt_clips: tails(b'T'),
            distinct_starts,
            depth: None,
            typing: None,
            mate_checks: Vec::new(),
            checked_mates: 0,
            mei_mates: 0,
            mate_identity: None,
            reference_copy: false,
            tail_mates: 0,
            assembly: None,
            genotype: None,
            qual: None,
            filters: Vec::new(),
            signals,
        };
        if candidate.support < options.min_support {
            candidate.filters.push(Filter::LowSupport);
        }
        if left_junction_clips + right_junction_clips < options.min_clips {
            candidate.filters.push(Filter::NoJunction);
        }
        if options.require_tail && candidate.polya_clips + candidate.polyt_clips == 0 {
            candidate.filters.push(Filter::NoTail);
        }
        candidate
    }

    /// Whether the candidate passes every filter.
    #[must_use]
    pub fn is_pass(&self) -> bool {
        self.filters.is_empty()
    }

    /// Whether signals support the insertion from both flanks.
    #[must_use]
    pub fn both_sides(&self) -> bool {
        let has = |side| self.signals.iter().any(|s| s.side == side);
        has(Side::Left) && has(Side::Right)
    }

    /// Supporting reads per read of local depth, once the depth is known.
    #[must_use]
    pub fn signal_ratio(&self) -> Option<f64> {
        #[allow(clippy::cast_precision_loss)]
        self.depth
            .map(|depth| self.support as f64 / depth.max(1) as f64)
    }

    /// Length of the target site duplication: with both junctions assembled,
    /// the assembled one (see [`Self::set_assembly`]); otherwise suggested by
    /// the junctions, the overlap of the two flanks, when each junction has
    /// at least `min_clips` soft-clips and they overlap by 1 to 50 bp.
    #[must_use]
    pub fn tsd_length(&self) -> Option<u64> {
        self.tsd
    }

    /// Keep the assembly of the junctions: the assembled junctions replace
    /// the estimated ones, and with both, the assembled TSD replaces the
    /// estimate (an assembly that finds no TSD removes it).
    pub fn set_assembly(&mut self, assembly: crate::assembly::Assembly) {
        if let Some(left) = &assembly.left {
            self.left_junction = Some(left.junction);
        }
        if let Some(right) = &assembly.right {
            self.right_junction = Some(right.junction);
        }
        if assembly.both() {
            self.tsd = assembly.tsd.as_ref().map(|t| t.len() as u64);
        }
        self.assembly = Some(assembly);
    }

    /// Whether the candidate failed `--min-support` (it is never written).
    #[must_use]
    pub fn is_low_support(&self) -> bool {
        self.filters.contains(&Filter::LowSupport)
    }

    /// Confidence level, for a candidate that passes the filters: `low` when
    /// typed by its mates only (no clip typed), else from the tail and TSD.
    #[must_use]
    pub fn confidence(&self) -> Option<Confidence> {
        if !self.is_pass() {
            None
        } else if self.polya_clips + self.polyt_clips == 0
            || self.typing.as_ref().is_some_and(|t| t.by_mates)
        {
            Some(Confidence::Low)
        } else if self.tsd.is_some() {
            Some(Confidence::High)
        } else if self.polya_clips + self.polyt_clips > 0 {
            Some(Confidence::Medium)
        } else {
            Some(Confidence::Low)
        }
    }
}

/// Most frequent value, the smallest on ties; `None` if empty.
fn mode(values: impl Iterator<Item = u64>) -> Option<u64> {
    let mut counts: BTreeMap<u64, usize> = BTreeMap::new();
    for value in values {
        *counts.entry(value).or_default() += 1;
    }
    let max = *counts.values().max()?;
    counts.into_iter().find(|&(_, n)| n == max).map(|(v, _)| v)
}

/// Breakpoint of a cluster: the most frequent soft-clip position when the
/// junctions have enough clips (`junctions`) or there is no mate signal;
/// otherwise a lone clip is as often noise, and the mate signals place the
/// insertion, within an interval ([`mate_breakpoint`]).
fn place(signals: &[Signal], junctions: bool) -> (u64, Option<(u64, u64)>) {
    let clip_breakpoint = mode(
        signals
            .iter()
            .filter(|s| matches!(s.kind, SignalKind::SoftClip { .. }))
            .map(|s| s.position),
    );
    let has_mates = signals.iter().any(|s| s.kind.is_mate_signal());
    match clip_breakpoint {
        Some(position) if junctions || !has_mates => (position, None),
        _ => {
            let (position, interval) = mate_breakpoint(signals);
            (position, Some(interval))
        }
    }
}

/// Breakpoint from mate signals only, and the interval where the insertion
/// lies. Reads of the upstream flank point downstream (their mates are in
/// the insertion, after them): the insertion lies after the furthest of
/// them; reads of the downstream flank before the nearest of them.
///
/// - Both groups, the upstream one ending before the downstream one starts
///   (within [`MATE_OVERLAP`]): between their inner edges, at the middle.
/// - Otherwise (one group, or groups that overlap: some signals are
///   background), the inner edge of the larger group, the insertion lying
///   within [`MATE_SPAN`] beyond it.
///
/// The interval reaches [`MATE_OVERLAP`] beyond the inner edges: a read can
/// end a few bases past the junction when the inserted bases match the
/// reference.
fn mate_breakpoint(signals: &[Signal]) -> (u64, (u64, u64)) {
    let positions = |side| {
        signals
            .iter()
            .filter(move |s| s.side == side && s.kind.is_mate_signal())
            .map(|s| s.position)
    };
    let upstream = positions(Side::Right).max();
    let downstream = positions(Side::Left).min();
    let (n_up, n_down) = (
        positions(Side::Right).count(),
        positions(Side::Left).count(),
    );
    match (upstream, downstream) {
        (Some(up), Some(down)) if up <= down + MATE_OVERLAP => {
            let (first, last) = (up.min(down), up.max(down));
            (
                first.midpoint(last),
                (first.saturating_sub(MATE_OVERLAP), last + MATE_OVERLAP),
            )
        }
        (Some(up), _) if n_up >= n_down => (up, (up.saturating_sub(MATE_OVERLAP), up + MATE_SPAN)),
        (_, Some(down)) => (down, (down.saturating_sub(MATE_SPAN), down + MATE_OVERLAP)),
        _ => {
            let position = signals[0].position;
            (position, (position, position))
        }
    }
}

/// Group sorted signals into candidates: consecutive signals of a contig at
/// most `window` bases apart belong to the same candidate. Candidates with
/// too little support are flagged [`Filter::LowSupport`].
#[must_use]
pub fn cluster(signals: Vec<Signal>, options: &ClusterOptions) -> Vec<Candidate> {
    let mut candidates = Vec::new();
    let mut current: Vec<Signal> = Vec::new();
    for signal in signals {
        if let Some(last) = current.last()
            && (last.contig != signal.contig || signal.position - last.position > options.window)
        {
            candidates.push(Candidate::new(std::mem::take(&mut current), options));
        }
        current.push(signal);
    }
    if !current.is_empty() {
        candidates.push(Candidate::new(current, options));
    }
    candidates
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::reader::signals::PolyTail;

    pub(crate) fn clip(position: u64, side: Side, read: &str, tail: Option<u8>) -> Signal {
        Signal {
            contig: 0,
            position,
            side,
            read: read.to_owned(),
            first_segment: true,
            kind: SignalKind::SoftClip {
                sequence: b"ACGT".to_vec(),
                median_quality: 30,
                tail: tail.map(|base| PolyTail { base, length: 20 }),
                supplementary: false,
                short: false,
            },
            reverse: side == Side::Left,
            start: if side == Side::Left {
                position
            } else {
                position - 80
            },
            end: if side == Side::Left {
                position + 80
            } else {
                position
            },
            mapq: Some(60),
        }
    }

    pub(crate) fn mate(position: u64, side: Side, read: &str, kind: SignalKind) -> Signal {
        Signal {
            kind,
            first_segment: false,
            ..clip(position, side, read, None)
        }
    }

    fn sorted(mut signals: Vec<Signal>) -> Vec<Signal> {
        signals.sort();
        signals
    }

    #[test]
    fn groups_by_window_and_contig() {
        let mut far = clip(1_000, Side::Left, "c", None);
        far.contig = 1;
        let signals = sorted(vec![
            clip(1_000, Side::Left, "a", None),
            clip(1_100, Side::Left, "b", None),
            clip(1_201, Side::Left, "c", None),
            far,
        ]);
        let candidates = cluster(signals, &ClusterOptions::default());
        let spans: Vec<_> = candidates
            .iter()
            .map(|c| (c.contig, c.start, c.end, c.signals.len()))
            .collect();
        assert_eq!(
            spans,
            [
                (0, 1_000, 1_100, 2),
                (0, 1_201, 1_201, 1),
                (1, 1_000, 1_000, 1)
            ]
        );
    }

    #[test]
    fn insertion_with_tsd() {
        let mut signals = Vec::new();
        for i in 0..5 {
            signals.push(clip(1_012, Side::Right, &format!("r{i}"), None));
            signals.push(clip(1_000, Side::Left, &format!("l{i}"), Some(b'A')));
        }
        signals.push(clip(1_003, Side::Left, "x", None));
        signals.push(mate(950, Side::Right, "m1", SignalKind::MateUnmapped));
        signals.push(mate(
            1_050,
            Side::Left,
            "m2",
            SignalKind::MateLowMapq {
                mapq: 0,
                mate: None,
            },
        ));
        // A read with both a clip and a mate signal counts once.
        signals.push(mate(1_000, Side::Left, "l0", SignalKind::MateUnmapped));
        let candidates = cluster(sorted(signals), &ClusterOptions::default());
        assert_eq!(candidates.len(), 1);
        let c = &candidates[0];
        assert_eq!(c.breakpoint, 1_000);
        assert_eq!(c.left_junction, Some(1_012));
        assert_eq!(c.right_junction, Some(1_000));
        assert_eq!(c.tsd_length(), Some(12));
        assert_eq!(c.support, 13);
        assert_eq!((c.left_clips, c.right_clips), (6, 5));
        assert_eq!((c.unmapped_mates, c.discordant_mates), (2, 1));
        assert_eq!((c.polya_clips, c.polyt_clips), (5, 0));
        assert!(c.both_sides());
        assert!(c.is_pass());
        assert_eq!((c.start, c.end), (950, 1_050));
    }

    #[test]
    fn breakpoint_from_mates_only() {
        let signals = sorted(vec![
            mate(900, Side::Right, "a", SignalKind::MateUnmapped),
            mate(980, Side::Right, "b", SignalKind::MateUnmapped),
            mate(1_020, Side::Left, "c", SignalKind::MateUnmapped),
        ]);
        let c = &cluster(signals, &ClusterOptions::default())[0];
        assert_eq!(c.breakpoint, 1_000);
        assert_eq!((c.left_junction, c.right_junction), (None, None));
        assert_eq!(c.tsd_length(), None);

        let one_side = sorted(vec![
            mate(900, Side::Right, "a", SignalKind::MateUnmapped),
            mate(980, Side::Right, "b", SignalKind::MateUnmapped),
        ]);
        let c = &cluster(one_side, &ClusterOptions::default())[0];
        assert_eq!(c.breakpoint, 980);
        assert!(!c.both_sides());
        // Beyond the inner edge of the upstream reads, within a fragment.
        assert_eq!(c.mate_interval, Some((980 - MATE_OVERLAP, 980 + MATE_SPAN)));
    }

    #[test]
    fn breakpoint_from_overlapping_mate_groups() {
        // The downstream reads start before the upstream ones end (some are
        // background): the inner edge of the larger group.
        let signals = sorted(vec![
            mate(1_000, Side::Left, "x", SignalKind::MateUnmapped),
            mate(1_100, Side::Right, "a", SignalKind::MateUnmapped),
            mate(1_150, Side::Right, "b", SignalKind::MateUnmapped),
            mate(1_170, Side::Right, "c", SignalKind::MateUnmapped),
        ]);
        let c = &cluster(signals, &ClusterOptions::default())[0];
        assert_eq!(c.breakpoint, 1_170);
        assert_eq!(c.mate_interval, Some((1_150, 1_470)));
        // Downstream reads only: before the nearest one.
        let signals = sorted(vec![
            mate(1_100, Side::Left, "a", SignalKind::MateUnmapped),
            mate(1_140, Side::Left, "b", SignalKind::MateUnmapped),
            mate(1_190, Side::Right, "x", SignalKind::MateUnmapped),
            mate(1_200, Side::Left, "c", SignalKind::MateUnmapped),
        ]);
        let c = &cluster(signals, &ClusterOptions::default())[0];
        assert_eq!(c.breakpoint, 1_100);
        assert_eq!(c.mate_interval, Some((800, 1_120)));
        // Consistent groups: between them, with the margin.
        let signals = sorted(vec![
            mate(980, Side::Right, "a", SignalKind::MateUnmapped),
            mate(1_020, Side::Left, "b", SignalKind::MateUnmapped),
        ]);
        let c = &cluster(signals, &ClusterOptions::default())[0];
        assert_eq!(c.mate_interval, Some((960, 1_040)));
    }

    #[test]
    fn lone_clip_with_mates() {
        // One clip is not enough junction clips: the mates place it.
        let signals = sorted(vec![
            clip(1_120, Side::Left, "z", None),
            mate(1_000, Side::Right, "a", SignalKind::MateUnmapped),
            mate(1_050, Side::Right, "b", SignalKind::MateUnmapped),
        ]);
        let candidates = cluster(signals, &ClusterOptions::default());
        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].breakpoint, 1_050);
        assert!(candidates[0].mate_interval.is_some());
        // Without mates, the clip does.
        let c = &cluster(
            sorted(vec![clip(1_120, Side::Left, "z", None)]),
            &ClusterOptions::default(),
        )[0];
        assert_eq!((c.breakpoint, c.mate_interval), (1_120, None));
    }

    #[test]
    fn low_support() {
        let signals = sorted(vec![
            clip(1_000, Side::Left, "a", None),
            mate(1_000, Side::Left, "a", SignalKind::MateUnmapped),
            clip(1_000, Side::Left, "b", None),
        ]);
        let lenient = ClusterOptions {
            require_tail: false,
            ..ClusterOptions::default()
        };
        let c = &cluster(signals.clone(), &lenient)[0];
        assert_eq!(c.support, 2);
        assert_eq!(c.filters, [Filter::LowSupport]);
        let options = ClusterOptions {
            min_support: 2,
            ..lenient
        };
        assert!(cluster(signals, &options)[0].is_pass());
    }

    #[test]
    fn distinct_starts_and_ratio() {
        let mut duplicate = clip(1_000, Side::Left, "b", None);
        duplicate.read = "b".to_owned();
        let signals = sorted(vec![
            clip(1_000, Side::Left, "a", None),
            duplicate,
            clip(1_012, Side::Right, "c", None),
        ]);
        let mut c = cluster(signals, &ClusterOptions::default()).remove(0);
        assert_eq!(c.distinct_starts, 2);
        assert_eq!(c.signal_ratio(), None);
        c.depth = Some(30);
        assert_eq!(c.signal_ratio(), Some(0.1));
        c.depth = Some(0);
        assert_eq!(c.signal_ratio(), Some(3.0));
    }

    #[test]
    fn no_junction_and_tsd_need_enough_clips() {
        let options = ClusterOptions {
            require_tail: false,
            ..ClusterOptions::default()
        };
        // Mates only.
        let mates = sorted(
            (0..5)
                .map(|i| {
                    mate(
                        1_000,
                        Side::Right,
                        &format!("m{i}"),
                        SignalKind::MateUnmapped,
                    )
                })
                .collect(),
        );
        let c = &cluster(mates.clone(), &options)[0];
        assert_eq!(c.filters, [Filter::NoJunction]);
        assert_eq!(c.confidence(), None);

        // One clip is not enough; two (within the tolerance) are.
        let mut one = mates.clone();
        one.push(clip(1_000, Side::Right, "c1", None));
        let c = &cluster(sorted(one), &options)[0];
        assert_eq!(c.filters, [Filter::NoJunction]);
        let mut two = mates;
        two.push(clip(1_000, Side::Right, "c1", None));
        two.push(clip(1_002, Side::Right, "c2", None));
        let c = &cluster(sorted(two), &options)[0];
        assert!(c.is_pass());
        assert_eq!(c.left_junction_clips, 2);
        assert_eq!(c.confidence(), Some(Confidence::Low));

        // A clip 3 bp away does not support the junction.
        let far = sorted(vec![
            clip(1_000, Side::Right, "a", None),
            clip(1_003, Side::Right, "b", None),
            clip(1_020, Side::Left, "c", None),
        ]);
        let c = &cluster(far, &options)[0];
        assert_eq!((c.left_junction_clips, c.right_junction_clips), (1, 1));
        assert!(c.is_pass());

        // TSD: each junction needs `min_clips` clips.
        let mut signals = vec![
            clip(1_012, Side::Right, "r0", None),
            clip(1_012, Side::Right, "r1", None),
            clip(1_000, Side::Left, "l0", Some(b'A')),
        ];
        let c = &cluster(sorted(signals.clone()), &options)[0];
        assert_eq!(c.tsd_length(), None);
        signals.push(clip(1_000, Side::Left, "l1", Some(b'A')));
        let c = &cluster(sorted(signals), &options)[0];
        assert_eq!(c.tsd_length(), Some(12));
    }

    #[test]
    fn confidence_levels() {
        let options = ClusterOptions::default();
        let mut signals = Vec::new();
        for i in 0..3 {
            signals.push(clip(1_012, Side::Right, &format!("r{i}"), None));
            signals.push(clip(1_000, Side::Left, &format!("l{i}"), Some(b'A')));
        }
        let high = &cluster(sorted(signals.clone()), &options)[0];
        assert_eq!(high.confidence(), Some(Confidence::High));

        // A tail without a TSD: medium, with or without mates.
        let mut one_junction: Vec<Signal> = (0..3)
            .map(|i| clip(1_000, Side::Left, &format!("l{i}"), Some(b'T')))
            .collect();
        let medium = &cluster(sorted(one_junction.clone()), &options)[0];
        assert_eq!(medium.tsd_length(), None);
        assert_eq!(medium.confidence(), Some(Confidence::Medium));
        one_junction.push(mate(1_050, Side::Left, "m", SignalKind::MateUnmapped));
        let medium = &cluster(sorted(one_junction), &options)[0];
        assert_eq!(medium.confidence(), Some(Confidence::Medium));

        // No tail: low, even with a TSD.
        let untailed: Vec<Signal> = (0..3)
            .flat_map(|i| {
                [
                    clip(1_012, Side::Right, &format!("r{i}"), None),
                    clip(1_000, Side::Left, &format!("l{i}"), None),
                ]
            })
            .collect();
        let low = &cluster(sorted(untailed), &options)[0];
        assert_eq!(low.tsd_length(), Some(12));
        assert_eq!(low.confidence(), Some(Confidence::Low));
    }

    #[test]
    fn no_tail() {
        let untailed = || {
            sorted(
                (0..3)
                    .flat_map(|i| {
                        [
                            clip(1_012, Side::Right, &format!("r{i}"), None),
                            clip(1_000, Side::Left, &format!("l{i}"), None),
                        ]
                    })
                    .chain([mate(950, Side::Right, "m", SignalKind::MateUnmapped)])
                    .collect(),
            )
        };
        // --require-tail.
        let options = ClusterOptions {
            require_tail: true,
            ..ClusterOptions::default()
        };
        let c = &cluster(untailed(), &options)[0];
        assert_eq!(c.filters, [Filter::NoTail]);
        assert_eq!(c.confidence(), None);

        // By default the tail is not required, and the TSD alone is not
        // `high`.
        let c = &cluster(untailed(), &ClusterOptions::default())[0];
        assert_eq!(c.tsd_length(), Some(12));
        assert_eq!(c.confidence(), Some(Confidence::Low));
    }

    #[test]
    fn mode_prefers_the_smallest_on_ties() {
        assert_eq!(mode([5, 3, 5, 3, 7].into_iter()), Some(3));
        assert_eq!(mode([5, 5, 3].into_iter()), Some(5));
        assert_eq!(mode(std::iter::empty()), None);
    }
}
