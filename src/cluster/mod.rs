//! Phase 2: spatial clustering of the signals into candidate insertion sites,
//! breakpoint estimation and capture artefact filters ([`filters`]).
//!
//! Positions are *boundaries* (see [`crate::reader::signals`]): a junction at
//! boundary `b` means the insertion lies between bases `b - 1` and `b`
//! (0-based), i.e. right after base `b` in 1-based coordinates.

pub mod filters;

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
    /// Require a soft-clip with a poly(A/T) tail: until element typing, the
    /// only evidence that the inserted sequence is a retrotransposon.
    pub require_tail: bool,
}

impl Default for ClusterOptions {
    fn default() -> Self {
        Self {
            window: 100,
            min_support: 3,
            min_clips: 2,
            require_tail: true,
        }
    }
}

/// Distance within which a soft-clip supports a junction.
pub const JUNCTION_TOLERANCE: u64 = 2;

/// Reason a candidate is filtered out.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Filter {
    /// Fewer distinct supporting reads than `--min-support`.
    LowSupport,
    /// Fewer soft-clips at the junctions than `--min-clips`: mate signals
    /// alone do not resolve the insertion to the base.
    NoJunction,
    /// No soft-clip with a poly(A/T) tail.
    NoTail,
    /// Every soft-clip stops at a probe edge, without a poly(A/T) tail or a
    /// TSD (capture artefact).
    ProbeEdge,
    /// Too few supporting reads for the local depth.
    LowRatio,
}

impl Filter {
    /// Name, as written in the outputs (and later in the VCF `FILTER`).
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::LowSupport => "low_support",
            Self::NoJunction => "no_junction",
            Self::NoTail => "no_tail",
            Self::ProbeEdge => "probe_edge",
            Self::LowRatio => "low_ratio",
        }
    }
}

/// Confidence level of a candidate that passes the filters.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Confidence {
    /// Soft-clips with a poly(A/T) tail and a TSD: complete evidence of a
    /// retrotransposition.
    High,
    /// Soft-clips and mate signals, without a TSD (or without a tail when
    /// tails are not required).
    Medium,
    /// Soft-clips only, without a TSD.
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
    /// junction, or else estimated from the mate signals.
    pub breakpoint: u64,
    /// End of the upstream flank: most supported junction of the clips on
    /// the right end of reads, if any.
    pub left_junction: Option<u64>,
    /// Start of the downstream flank: most supported junction of the clips
    /// on the left end of reads, if any.
    pub right_junction: Option<u64>,
    /// Soft-clips supporting the left junction (within
    /// [`JUNCTION_TOLERANCE`]).
    pub left_junction_clips: usize,
    /// Soft-clips supporting the right junction (within
    /// [`JUNCTION_TOLERANCE`]).
    pub right_junction_clips: usize,
    /// Length of the target site duplication suggested by the junctions.
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
        let near = |junction: Option<u64>, side: Side| {
            junction.map_or(0, |j| {
                clips(side)
                    .filter(|s| s.position.abs_diff(j) <= JUNCTION_TOLERANCE)
                    .count()
            })
        };
        let left_junction_clips = near(left_junction, Side::Right);
        let right_junction_clips = near(right_junction, Side::Left);
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
        let breakpoint = mode(
            signals
                .iter()
                .filter(|s| matches!(s.kind, SignalKind::SoftClip { .. }))
                .map(|s| s.position),
        )
        .unwrap_or_else(|| mate_breakpoint(&signals));

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
            left_junction,
            right_junction,
            left_junction_clips,
            right_junction_clips,
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

    /// Length of the target site duplication suggested by the junctions: the
    /// overlap of the two flanks, when each junction has at least
    /// `min_clips` soft-clips and they overlap by 1 to 50 bp.
    #[must_use]
    pub fn tsd_length(&self) -> Option<u64> {
        self.tsd
    }

    /// Whether the candidate failed `--min-support` (it is never written).
    #[must_use]
    pub fn is_low_support(&self) -> bool {
        self.filters.contains(&Filter::LowSupport)
    }

    /// Confidence level, for a candidate that passes the filters.
    #[must_use]
    pub fn confidence(&self) -> Option<Confidence> {
        if !self.is_pass() {
            None
        } else if self.tsd.is_some() && self.polya_clips + self.polyt_clips > 0 {
            Some(Confidence::High)
        } else if self.unmapped_mates + self.discordant_mates > 0 {
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

/// Breakpoint from mate signals only: between the furthest read of the
/// upstream flank and the nearest read of the downstream flank.
fn mate_breakpoint(signals: &[Signal]) -> u64 {
    let upstream = signals
        .iter()
        .filter(|s| s.side == Side::Right)
        .map(|s| s.position)
        .max();
    let downstream = signals
        .iter()
        .filter(|s| s.side == Side::Left)
        .map(|s| s.position)
        .min();
    match (upstream, downstream) {
        (Some(up), Some(down)) => up.midpoint(down),
        (Some(position), None) | (None, Some(position)) => position,
        (None, None) => signals[0].position,
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
            SignalKind::MateLowMapq { mapq: 0 },
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
        assert_eq!(c.confidence(), Some(Confidence::Medium));

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

        // A tail without a TSD: medium with mates, low without.
        let mut one_junction: Vec<Signal> = (0..3)
            .map(|i| clip(1_000, Side::Left, &format!("l{i}"), Some(b'T')))
            .collect();
        let low = &cluster(sorted(one_junction.clone()), &options)[0];
        assert_eq!(low.tsd_length(), None);
        assert_eq!(low.confidence(), Some(Confidence::Low));
        one_junction.push(mate(1_050, Side::Left, "m", SignalKind::MateUnmapped));
        let medium = &cluster(sorted(one_junction), &options)[0];
        assert_eq!(medium.confidence(), Some(Confidence::Medium));
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
        let c = &cluster(untailed(), &ClusterOptions::default())[0];
        assert_eq!(c.filters, [Filter::NoTail]);
        assert_eq!(c.confidence(), None);

        // --allow-no-tail: the TSD alone is not `high`.
        let options = ClusterOptions {
            require_tail: false,
            ..ClusterOptions::default()
        };
        let c = &cluster(untailed(), &options)[0];
        assert_eq!(c.tsd_length(), Some(12));
        assert_eq!(c.confidence(), Some(Confidence::Medium));
    }

    #[test]
    fn mode_prefers_the_smallest_on_ties() {
        assert_eq!(mode([5, 3, 5, 3, 7].into_iter()), Some(3));
        assert_eq!(mode([5, 5, 3].into_iter()), Some(5));
        assert_eq!(mode(std::iter::empty()), None);
    }
}
