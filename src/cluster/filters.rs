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
}

impl Default for FilterOptions {
    fn default() -> Self {
        Self {
            min_signal_ratio: 0.05,
        }
    }
}

/// Flag a candidate whose depth is known.
///
/// - [`Filter::ProbeEdge`]: it has soft-clips and every one stops exactly at
///   the start or end of a probe (digestion or ligation artefact of the
///   capture).
/// - [`Filter::LowRatio`]: supporting reads per read of local depth below
///   the threshold.
pub fn apply(
    candidate: &mut Candidate,
    contig: &str,
    targets: &TargetSet,
    options: &FilterOptions,
) {
    let mut clips = candidate
        .signals
        .iter()
        .filter(|s| matches!(s.kind, SignalKind::SoftClip { .. }))
        .peekable();
    let at_probe_edges = clips.peek().is_some()
        && clips.all(|s| targets.nearest_probe_edge(contig, s.position, 0) == Some(0));
    if at_probe_edges {
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
