//! Filters specific to capture data: probe-edge artefacts and the ratio of
//! supporting reads to the local depth.

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
}

impl Default for FilterOptions {
    fn default() -> Self {
        Self {
            min_signal_ratio: 0.05,
            probe_edge_tolerance: 2,
        }
    }
}

/// Flag a candidate whose depth is known.
///
/// - [`Filter::ProbeEdge`]: it has soft-clips, every one stops within the
///   tolerance of a probe start or end, none has a poly(A/T) tail and there
///   is no TSD (digestion or ligation artefact of the capture). A real
///   insertion next to a probe edge keeps its hallmarks.
/// - [`Filter::LowRatio`]: supporting reads per read of local depth below
///   the threshold.
pub fn apply(
    candidate: &mut Candidate,
    contig: &str,
    targets: &TargetSet,
    options: &FilterOptions,
) {
    let hallmarks =
        candidate.polya_clips + candidate.polyt_clips > 0 || candidate.tsd_length().is_some();
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
        let options = FilterOptions::default();
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

    /// Re-apply the filters from scratch.
    fn apply_fresh(c: &mut Candidate, options: &FilterOptions) {
        c.filters.clear();
        apply(c, "chr1", &targets(), options);
    }

    #[test]
    fn probe_edge_spares_insertion_hallmarks() {
        let options = FilterOptions::default();
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
        let options = FilterOptions::default();
        let mut c = candidate(signals.clone(), 100);
        apply(&mut c, "chr1", &targets(), &options);
        assert!(c.is_pass());
        let mut c = candidate(signals, 101);
        apply(&mut c, "chr1", &targets(), &options);
        assert_eq!(c.filters, [Filter::LowRatio]);
    }
}
