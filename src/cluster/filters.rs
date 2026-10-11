//! Filters applied once the candidate is typed and its depth known: element
//! typing, probe-edge artefacts of capture data and the ratio of supporting
//! reads to the local depth.

use super::{Candidate, Filter};
use crate::reader::signals::SignalKind;
use crate::regions::targets::TargetSet;

/// Capture filter options.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FilterOptions {
    /// Minimum number of supporting reads per read of local depth.
    pub min_signal_ratio: f64,
    /// Maximum distance, in bp, between a soft-clip and a probe start or end
    /// for the clip to count as stopping at the probe edge.
    pub probe_edge_tolerance: u64,
    /// Minimum number of soft-clips at the junctions (as
    /// [`super::ClusterOptions::min_clips`]), here without the short ones
    /// when no mate is an MEI mate.
    pub min_clips: usize,
    /// Minimum number of junction soft-clips typed to an MEI consensus.
    pub min_typed_clips: usize,
    /// Minimum fraction of MEI mates among the checked mates.
    pub min_mei_mate_fraction: f64,
    /// Minimum quality score (0: no minimum).
    pub min_qual: u32,
}

impl Default for FilterOptions {
    fn default() -> Self {
        Self {
            min_signal_ratio: 0.05,
            probe_edge_tolerance: 2,
            min_clips: 2,
            min_typed_clips: 2,
            min_mei_mate_fraction: 0.3,
            min_qual: 0,
        }
    }
}

/// Flag a candidate once typed and its depth known.
///
/// Short soft-clips count only when at least one mate is an MEI mate: alone,
/// they are too often a piece of a nearby repeat (#32). Without one, they are
/// left out of the junction clips ([`Filter::NoJunction`], against
/// `min_clips`) and of the typed clips.
///
/// - [`Filter::Untyped`]: fewer junction soft-clips typed to an MEI consensus
///   than the minimum.
/// - [`Filter::NonMeiMates`]: some mates are checked, and MEI mates are a
///   smaller fraction of them than the minimum.
/// - [`Filter::ProbeEdge`]: it has soft-clips, every one stops within the
///   tolerance of a probe start or end, none has a poly(A/T) tail, none is
///   typed and there is no TSD (digestion or ligation artefact of the
///   capture). A real insertion next to a probe edge keeps its hallmarks.
/// - [`Filter::LowQual`]: quality score ([`crate::genotype::score`]) below
///   the minimum, also set as [`Candidate::qual`].
/// - [`Filter::LowRatio`]: supporting reads per read of local depth below
///   the threshold.
pub fn apply(
    candidate: &mut Candidate,
    contig: &str,
    targets: &TargetSet,
    options: &FilterOptions,
) {
    let typed_clips = candidate.typing.as_ref().map_or(0, |t| t.typed_clips);
    let short_counts = candidate.mei_mates > 0;
    let junction_clips = candidate.left_junction_clips + candidate.right_junction_clips;
    if !short_counts
        && junction_clips - candidate.short_junction_clips < options.min_clips
        && !candidate.filters.contains(&Filter::NoJunction)
    {
        // In filter order: only `low_support` comes before.
        let at = candidate
            .filters
            .partition_point(|f| *f < Filter::NoJunction);
        candidate.filters.insert(at, Filter::NoJunction);
    }
    let counted_typed_clips = if short_counts {
        typed_clips
    } else {
        typed_clips - candidate.typing.as_ref().map_or(0, |t| t.short_typed_clips)
    };
    // One typed clip is enough with a poly(A/T) tail and as many MEI mates
    // as typed clips required: on short reads, few clips are long enough to
    // be typed (#32).
    let tail = candidate.polya_clips + candidate.polyt_clips > 0;
    let typed_with_mates =
        counted_typed_clips >= 1 && tail && candidate.mei_mates >= options.min_typed_clips;
    if counted_typed_clips < options.min_typed_clips && !typed_with_mates {
        candidate.filters.push(Filter::Untyped);
    }
    #[allow(clippy::cast_precision_loss)]
    if candidate.checked_mates > 0
        && (candidate.mei_mates as f64)
            < options.min_mei_mate_fraction * candidate.checked_mates as f64
    {
        candidate.filters.push(Filter::NonMeiMates);
    }
    let hallmarks = candidate.polya_clips + candidate.polyt_clips > 0
        || candidate.tsd_length().is_some()
        || typed_clips > 0;
    let mut clips = candidate
        .signals
        .iter()
        .filter(|s| matches!(s.kind, SignalKind::SoftClip { .. }))
        .peekable();
    let at_probe_edges = clips.peek().is_some()
        && clips.all(|s| {
            targets
                .nearest_probe_edge(contig, s.position, options.probe_edge_tolerance)
                .is_some()
        });
    if at_probe_edges && !hallmarks {
        candidate.filters.push(Filter::ProbeEdge);
    }
    if candidate
        .signal_ratio()
        .is_some_and(|ratio| ratio < options.min_signal_ratio)
    {
        candidate.filters.push(Filter::LowRatio);
    }
    let qual = crate::genotype::score::quality(candidate);
    candidate.qual = Some(qual);
    if qual < options.min_qual {
        candidate.filters.push(Filter::LowQual);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cluster::tests::{clip, mate};
    use crate::cluster::{ClusterOptions, cluster};
    use crate::reader::signals::Side;
    use crate::regions::bed::Probe;
    use crate::regions::targets::Contig;

    fn targets() -> TargetSet {
        let probe = Probe {
            contig: "chr1".to_owned(),
            start: 1_000,
            end: 1_200,
            name: None,
        };
        let contig = Contig {
            name: "chr1".to_owned(),
            length: 10_000,
        };
        TargetSet::new(vec![probe], 100, &[contig]).unwrap()
    }

    fn candidate(signals: Vec<crate::reader::signals::Signal>, depth: u64) -> Candidate {
        let mut signals = signals;
        signals.sort();
        // Tails are not what these tests are about.
        let options = ClusterOptions {
            require_tail: false,
            ..ClusterOptions::default()
        };
        let mut c = cluster(signals, &options).remove(0);
        c.depth = Some(depth);
        c
    }

    #[test]
    fn probe_edge() {
        let options = untyped_ok();
        let mut c = candidate(
            (0..5)
                .map(|i| clip(1_000, Side::Left, &format!("r{i}"), None))
                .collect(),
            20,
        );
        apply(&mut c, "chr1", &targets(), &options);
        assert_eq!(c.filters, [Filter::ProbeEdge]);

        // The probe end, too.
        let mut c = candidate(
            (0..5)
                .map(|i| clip(1_200, Side::Right, &format!("r{i}"), None))
                .collect(),
            20,
        );
        apply(&mut c, "chr1", &targets(), &options);
        assert_eq!(c.filters, [Filter::ProbeEdge]);

        // Within the tolerance (2 bp by default), but not beyond.
        let mut c = candidate(
            (0..5)
                .map(|i| clip(1_002, Side::Left, &format!("r{i}"), None))
                .collect(),
            20,
        );
        apply(&mut c, "chr1", &targets(), &options);
        assert_eq!(c.filters, [Filter::ProbeEdge]);
        let exact = FilterOptions {
            probe_edge_tolerance: 0,
            ..options
        };
        apply_fresh(&mut c, &exact);
        assert!(c.is_pass());

        // One clip elsewhere is enough to keep it.
        let mut signals: Vec<_> = (0..5)
            .map(|i| clip(1_000, Side::Left, &format!("r{i}"), None))
            .collect();
        signals.push(clip(1_003, Side::Left, "x", None));
        let mut c = candidate(signals, 20);
        apply(&mut c, "chr1", &targets(), &options);
        assert!(c.is_pass());

        // Mate signals only: not a probe-edge artefact (but no junction).
        let mut c = candidate(
            (0..5)
                .map(|i| {
                    mate(
                        1_000,
                        Side::Left,
                        &format!("r{i}"),
                        SignalKind::MateUnmapped,
                    )
                })
                .collect(),
            20,
        );
        apply(&mut c, "chr1", &targets(), &options);
        assert_eq!(c.filters, [Filter::NoJunction]);
    }

    /// Default options, without the typing requirement (these tests are
    /// not about it).
    fn untyped_ok() -> FilterOptions {
        FilterOptions {
            min_typed_clips: 0,
            ..FilterOptions::default()
        }
    }

    #[test]
    fn non_mei_mates() {
        let signals: Vec<_> = (0..5)
            .map(|i| clip(1_100, Side::Left, &format!("r{i}"), Some(b'A')))
            .collect();
        let mut c = candidate(signals, 20);
        let options = untyped_ok();
        // Nothing checked: not flagged.
        apply(&mut c, "chr1", &targets(), &options);
        assert!(c.is_pass());
        // 3 MEI mates of 10: at the threshold; 2 of 10: below.
        c.checked_mates = 10;
        c.mei_mates = 3;
        apply_fresh(&mut c, &options);
        assert!(c.is_pass());
        c.mei_mates = 2;
        apply_fresh(&mut c, &options);
        assert_eq!(c.filters, [Filter::NonMeiMates]);
        let off = FilterOptions {
            min_mei_mate_fraction: 0.0,
            ..options
        };
        apply_fresh(&mut c, &off);
        assert!(c.is_pass());
    }

    #[test]
    fn untyped() {
        let signals: Vec<_> = (0..5)
            .map(|i| clip(1_100, Side::Left, &format!("r{i}"), Some(b'A')))
            .collect();
        let mut c = candidate(signals, 20);
        apply(&mut c, "chr1", &targets(), &FilterOptions::default());
        assert_eq!(c.filters, [Filter::Untyped]);
        c.filters.clear();
        c.typing = Some(alu_typing(2, 0));
        apply(&mut c, "chr1", &targets(), &FilterOptions::default());
        assert!(c.is_pass(), "{:?}", c.filters);
    }

    fn alu_typing(typed_clips: usize, short_typed_clips: usize) -> crate::typer::Typing {
        crate::typer::Typing {
            family: crate::typer::consensus::Family::Alu,
            subfamily: "AluY".to_owned(),
            reverse: false,
            element_start: 250,
            element_end: 281,
            identity: 1.0,
            typed_clips,
            short_typed_clips,
        }
    }

    fn short_clip(position: u64, side: Side, read: &str) -> crate::reader::signals::Signal {
        let mut signal = clip(position, side, read, None);
        if let SignalKind::SoftClip { short, .. } = &mut signal.kind {
            *short = true;
        }
        signal
    }

    #[test]
    fn short_clips_need_mei_mates() {
        // One long and two short junction clips, all typed.
        let mut c = candidate(
            vec![
                clip(1_100, Side::Left, "long", None),
                short_clip(1_100, Side::Left, "s1"),
                short_clip(1_100, Side::Left, "s2"),
            ],
            20,
        );
        assert_eq!(c.short_junction_clips, 2);
        c.typing = Some(alu_typing(3, 2));
        let options = FilterOptions::default();
        // Without an MEI mate, the short clips do not count.
        apply(&mut c, "chr1", &targets(), &options);
        assert_eq!(c.filters, [Filter::NoJunction, Filter::Untyped]);
        // With one, they do.
        c.checked_mates = 1;
        c.mei_mates = 1;
        apply_fresh(&mut c, &options);
        assert!(c.is_pass(), "{:?}", c.filters);
    }

    #[test]
    fn one_typed_clip_with_a_tail_and_mei_mates() {
        let signals: Vec<_> = (0..3)
            .map(|i| clip(1_100, Side::Left, &format!("r{i}"), Some(b'A')))
            .collect();
        let mut c = candidate(signals, 20);
        c.typing = Some(alu_typing(1, 0));
        let options = FilterOptions::default();
        apply(&mut c, "chr1", &targets(), &options);
        assert_eq!(c.filters, [Filter::Untyped]);
        // As many MEI mates as typed clips required.
        c.checked_mates = 2;
        c.mei_mates = 2;
        apply_fresh(&mut c, &options);
        assert!(c.is_pass(), "{:?}", c.filters);
        // Not without a tail.
        let mut c = candidate(
            (0..3)
                .map(|i| clip(1_100, Side::Left, &format!("r{i}"), None))
                .collect(),
            20,
        );
        c.typing = Some(alu_typing(1, 0));
        c.checked_mates = 2;
        c.mei_mates = 2;
        apply(&mut c, "chr1", &targets(), &options);
        assert_eq!(c.filters, [Filter::Untyped]);
    }

    /// Re-apply the filters from scratch.
    fn apply_fresh(c: &mut Candidate, options: &FilterOptions) {
        c.filters.clear();
        apply(c, "chr1", &targets(), options);
    }

    #[test]
    fn probe_edge_spares_insertion_hallmarks() {
        let options = untyped_ok();
        // A poly(A/T) tail at the edge: not an artefact.
        let mut signals: Vec<_> = (0..5)
            .map(|i| clip(1_000, Side::Left, &format!("r{i}"), None))
            .collect();
        signals.push(clip(1_000, Side::Left, "tail", Some(b'A')));
        let mut c = candidate(signals, 20);
        apply(&mut c, "chr1", &targets(), &options);
        assert!(c.is_pass(), "{:?}", c.filters);

        // A TSD between two junctions both within reach of the edge.
        let mut signals: Vec<_> = (0..3)
            .map(|i| clip(1_000, Side::Left, &format!("l{i}"), None))
            .collect();
        signals.extend((0..3).map(|i| clip(1_002, Side::Right, &format!("r{i}"), None)));
        let mut c = candidate(signals, 20);
        assert_eq!(c.tsd_length(), Some(2));
        apply(&mut c, "chr1", &targets(), &options);
        assert!(c.is_pass(), "{:?}", c.filters);
    }

    #[test]
    fn low_ratio() {
        let signals: Vec<_> = (0..5)
            .map(|i| clip(1_100, Side::Left, &format!("r{i}"), None))
            .collect();
        let options = untyped_ok();
        let mut c = candidate(signals.clone(), 100);
        apply(&mut c, "chr1", &targets(), &options);
        assert!(c.is_pass());
        let mut c = candidate(signals, 101);
        apply(&mut c, "chr1", &targets(), &options);
        assert_eq!(c.filters, [Filter::LowRatio]);
    }
}
