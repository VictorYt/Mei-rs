//! Capture targets: the probes of the BED file, padded and merged into the
//! regions to scan, with an interval tree per contig for probe queries.
//!
//! All coordinates are 0-based half-open. A *boundary* is a position between
//! two bases: boundary `b` lies between bases `b - 1` and `b`, so a probe
//! `[start, end)` has its edges at boundaries `start` and `end`.

use std::collections::HashMap;

use coitrees::{COITree, GenericInterval, Interval, IntervalTree};
use thiserror::Error;

use super::Region;
use super::bed::Probe;

/// A contig of the alignment header.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Contig {
    /// Name.
    pub name: String,
    /// Length in bases.
    pub length: u64,
}

/// A region to scan: padded probes merged when they overlap.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ScanRegion {
    /// Contig name.
    pub contig: String,
    /// First base, 0-based.
    pub start: u64,
    /// Base after the last one, 0-based.
    pub end: u64,
}

impl ScanRegion {
    /// Number of bases.
    #[must_use]
    pub fn len(&self) -> u64 {
        self.end - self.start
    }

    /// Whether the region holds no base.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.start >= self.end
    }
}

/// Errors raised while building a [`TargetSet`].
#[derive(Debug, Error)]
pub enum TargetError {
    /// A probe is on a contig absent from the alignment header.
    #[error("target {contig}:{start}-{end} is on a contig absent from the alignment header")]
    UnknownContig {
        /// Contig name.
        contig: String,
        /// Probe start, 0-based.
        start: u64,
        /// Probe end, 0-based exclusive.
        end: u64,
    },
    /// A probe starts after the end of its contig.
    #[error("target {contig}:{start}-{end} starts after the end of {contig} ({length} bp)")]
    BeyondContig {
        /// Contig name.
        contig: String,
        /// Probe start, 0-based.
        start: u64,
        /// Probe end, 0-based exclusive.
        end: u64,
        /// Contig length.
        length: u64,
    },
    /// A contig is too long for the interval tree (32-bit coordinates).
    #[error(
        "contig {contig} is too long ({length} bp): at most {} bp are supported",
        i32::MAX
    )]
    ContigTooLong {
        /// Contig name.
        contig: String,
        /// Contig length.
        length: u64,
    },
}

/// The capture probes and the regions to scan.
pub struct TargetSet {
    padding: u32,
    /// Probes sorted by contig (header order), start and end.
    probes: Vec<Probe>,
    /// Regions to scan, sorted like the probes and disjoint.
    regions: Vec<ScanRegion>,
    /// Rank of each contig in the header.
    ranks: HashMap<String, usize>,
    /// Probes of each contig, as indices into `probes` (inclusive ends).
    trees: HashMap<String, COITree<usize, u32>>,
}

impl std::fmt::Debug for TargetSet {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("TargetSet")
            .field("padding", &self.padding)
            .field("probes", &self.probes.len())
            .field("regions", &self.regions.len())
            .finish_non_exhaustive()
    }
}

impl TargetSet {
    /// Pad the probes by `padding` bases on each side, clamped to the contig
    /// bounds, and merge the padded probes that overlap or touch.
    ///
    /// Probe ends beyond the contig end are clamped.
    ///
    /// # Errors
    ///
    /// Fails if a probe is on a contig absent from `contigs`, starts after the
    /// end of its contig, or if a contig is longer than `i32::MAX`.
    pub fn new(
        mut probes: Vec<Probe>,
        padding: u32,
        contigs: &[Contig],
    ) -> Result<Self, TargetError> {
        let ranks: HashMap<String, usize> = contigs
            .iter()
            .enumerate()
            .map(|(rank, contig)| (contig.name.clone(), rank))
            .collect();
        for probe in &mut probes {
            let Some(&rank) = ranks.get(&probe.contig) else {
                return Err(TargetError::UnknownContig {
                    contig: probe.contig.clone(),
                    start: probe.start,
                    end: probe.end,
                });
            };
            let contig = &contigs[rank];
            if contig.length > i32::MAX as u64 {
                return Err(TargetError::ContigTooLong {
                    contig: contig.name.clone(),
                    length: contig.length,
                });
            }
            if probe.start >= contig.length {
                return Err(TargetError::BeyondContig {
                    contig: probe.contig.clone(),
                    start: probe.start,
                    end: probe.end,
                    length: contig.length,
                });
            }
            probe.end = probe.end.min(contig.length);
        }
        probes.sort_by(|a, b| {
            (ranks[&a.contig], a.start, a.end).cmp(&(ranks[&b.contig], b.start, b.end))
        });

        let mut regions: Vec<ScanRegion> = Vec::new();
        for probe in &probes {
            let length = contigs[ranks[&probe.contig]].length;
            let start = probe.start.saturating_sub(u64::from(padding));
            let end = (probe.end + u64::from(padding)).min(length);
            match regions.last_mut() {
                Some(last) if last.contig == probe.contig && start <= last.end => {
                    last.end = last.end.max(end);
                }
                _ => regions.push(ScanRegion {
                    contig: probe.contig.clone(),
                    start,
                    end,
                }),
            }
        }

        let mut intervals: HashMap<&str, Vec<Interval<usize>>> = HashMap::new();
        for (i, probe) in probes.iter().enumerate() {
            intervals
                .entry(probe.contig.as_str())
                .or_default()
                .push(Interval::new(to_i32(probe.start), to_i32(probe.end - 1), i));
        }
        let trees = intervals
            .into_iter()
            .map(|(contig, intervals)| (contig.to_owned(), COITree::new(&intervals)))
            .collect();

        Ok(Self {
            padding,
            probes,
            regions,
            ranks,
            trees,
        })
    }

    /// Padding added on each side of the probes.
    #[must_use]
    pub fn padding(&self) -> u32 {
        self.padding
    }

    /// Probes, sorted by contig (header order), start and end.
    #[must_use]
    pub fn probes(&self) -> &[Probe] {
        &self.probes
    }

    /// Regions to scan, sorted and disjoint.
    #[must_use]
    pub fn regions(&self) -> &[ScanRegion] {
        &self.regions
    }

    /// Number of bases covered by at least one probe.
    #[must_use]
    pub fn probe_bases(&self) -> u64 {
        let mut total = 0;
        let mut current: Option<(&str, u64, u64)> = None;
        for probe in &self.probes {
            match &mut current {
                Some((contig, _, end)) if *contig == probe.contig && probe.start <= *end => {
                    *end = (*end).max(probe.end);
                }
                _ => {
                    if let Some((_, start, end)) = current {
                        total += end - start;
                    }
                    current = Some((&probe.contig, probe.start, probe.end));
                }
            }
        }
        total + current.map_or(0, |(_, start, end)| end - start)
    }

    /// Number of bases to scan.
    #[must_use]
    pub fn scanned_bases(&self) -> u64 {
        self.regions.iter().map(ScanRegion::len).sum()
    }

    /// Keep only the parts of the regions to scan that fall in `region`.
    /// The probes are kept, for the probe-edge queries.
    pub fn restrict(&mut self, region: &Region) {
        self.regions.retain_mut(|scan| {
            if scan.contig != region.contig {
                return false;
            }
            if let Some(span) = region.span {
                scan.start = scan.start.max(span.start);
                scan.end = scan.end.min(span.end);
            }
            !scan.is_empty()
        });
    }

    /// Whether base `pos` of `contig` lies in a region to scan.
    #[must_use]
    pub fn in_scan_region(&self, contig: &str, pos: u64) -> bool {
        let Some(&rank) = self.ranks.get(contig) else {
            return false;
        };
        let key = |r: &ScanRegion| (self.ranks[&r.contig], r.start);
        let i = self.regions.partition_point(|r| key(r) <= (rank, pos));
        i > 0 && {
            let r = &self.regions[i - 1];
            r.contig == contig && pos < r.end
        }
    }

    /// Whether base `pos` of `contig` lies in a probe.
    #[must_use]
    pub fn in_probe(&self, contig: &str, pos: u64) -> bool {
        self.trees.get(contig).is_some_and(|tree| {
            let pos = to_i32(pos);
            tree.query_count(pos, pos) > 0
        })
    }

    /// Distance from `boundary` to the nearest probe edge of `contig`, if one
    /// lies within `max_distance` bases.
    #[must_use]
    pub fn nearest_probe_edge(
        &self,
        contig: &str,
        boundary: u64,
        max_distance: u64,
    ) -> Option<u64> {
        let tree = self.trees.get(contig)?;
        // A probe with an edge within reach holds a base in this window.
        let first = to_i32(boundary.saturating_sub(max_distance + 1));
        let last = to_i32(boundary.saturating_add(max_distance));
        let mut nearest: Option<u64> = None;
        tree.query(first, last, |node| {
            let probe = &self.probes[*GenericInterval::<usize>::metadata(node)];
            for edge in [probe.start, probe.end] {
                let distance = edge.abs_diff(boundary);
                if distance <= max_distance && nearest.is_none_or(|d| distance < d) {
                    nearest = Some(distance);
                }
            }
        });
        nearest
    }
}

/// Convert a position to the 32-bit coordinates of the interval tree,
/// saturating (contig lengths are checked in [`TargetSet::new`]).
fn to_i32(pos: u64) -> i32 {
    i32::try_from(pos).unwrap_or(i32::MAX)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn probe(contig: &str, start: u64, end: u64) -> Probe {
        Probe {
            contig: contig.to_owned(),
            start,
            end,
            name: None,
        }
    }

    fn contigs() -> Vec<Contig> {
        vec![
            Contig {
                name: "chr1".to_owned(),
                length: 10_000,
            },
            Contig {
                name: "chr2".to_owned(),
                length: 5_000,
            },
        ]
    }

    fn region(contig: &str, start: u64, end: u64) -> ScanRegion {
        ScanRegion {
            contig: contig.to_owned(),
            start,
            end,
        }
    }

    #[test]
    fn pads_clamps_sorts_and_merges() {
        let probes = vec![
            probe("chr2", 4_900, 5_200), // end beyond the contig: clamped
            probe("chr1", 1_000, 1_100),
            probe("chr1", 50, 100),      // padding clamped at 0
            probe("chr1", 1_300, 1_400), // padded, touches the previous one
            probe("chr1", 1_500, 1_600), // padded, overlaps the previous one
            probe("chr1", 3_000, 3_100),
        ];
        let set = TargetSet::new(probes, 100, &contigs()).unwrap();
        assert_eq!(
            set.regions(),
            [
                region("chr1", 0, 200),
                region("chr1", 900, 1_700),
                region("chr1", 2_900, 3_200),
                region("chr2", 4_800, 5_000),
            ]
        );
        let starts: Vec<_> = set
            .probes()
            .iter()
            .map(|p| (p.contig.as_str(), p.start))
            .collect();
        assert_eq!(
            starts,
            [
                ("chr1", 50),
                ("chr1", 1_000),
                ("chr1", 1_300),
                ("chr1", 1_500),
                ("chr1", 3_000),
                ("chr2", 4_900)
            ]
        );
        assert_eq!(set.probes()[5].end, 5_000);
        assert_eq!(set.probe_bases(), 50 + 100 + 100 + 100 + 100 + 100);
        assert_eq!(set.scanned_bases(), 200 + 800 + 300 + 200);
    }

    #[test]
    fn overlapping_probes_are_counted_once() {
        let probes = vec![probe("chr1", 100, 200), probe("chr1", 150, 300)];
        let set = TargetSet::new(probes, 0, &contigs()).unwrap();
        assert_eq!(set.probe_bases(), 200);
        assert_eq!(set.regions(), [region("chr1", 100, 300)]);
    }

    #[test]
    fn invalid_probes_are_rejected() {
        let err = TargetSet::new(vec![probe("chr3", 1, 2)], 0, &contigs()).unwrap_err();
        assert!(matches!(err, TargetError::UnknownContig { .. }));
        let err = TargetSet::new(vec![probe("chr2", 5_000, 5_100)], 0, &contigs()).unwrap_err();
        assert!(matches!(
            err,
            TargetError::BeyondContig { length: 5_000, .. }
        ));
        let long = [Contig {
            name: "big".to_owned(),
            length: 1 << 31,
        }];
        let err = TargetSet::new(vec![probe("big", 1, 2)], 0, &long).unwrap_err();
        assert!(matches!(err, TargetError::ContigTooLong { .. }));
    }

    #[test]
    fn position_queries() {
        let probes = vec![probe("chr1", 1_000, 1_100), probe("chr2", 10, 20)];
        let set = TargetSet::new(probes, 50, &contigs()).unwrap();

        assert!(set.in_probe("chr1", 1_000));
        assert!(set.in_probe("chr1", 1_099));
        assert!(!set.in_probe("chr1", 1_100));
        assert!(!set.in_probe("chr1", 999));
        assert!(!set.in_probe("chr3", 1_000));

        assert!(set.in_scan_region("chr1", 950));
        assert!(set.in_scan_region("chr1", 1_149));
        assert!(!set.in_scan_region("chr1", 949));
        assert!(!set.in_scan_region("chr1", 1_150));
        assert!(set.in_scan_region("chr2", 0));
        assert!(!set.in_scan_region("chr2", 70));
        assert!(!set.in_scan_region("chr3", 0));
    }

    #[test]
    fn nearest_probe_edge() {
        let probes = vec![probe("chr1", 1_000, 1_100), probe("chr1", 1_120, 1_200)];
        let set = TargetSet::new(probes, 0, &contigs()).unwrap();
        assert_eq!(set.nearest_probe_edge("chr1", 1_000, 5), Some(0));
        assert_eq!(set.nearest_probe_edge("chr1", 997, 5), Some(3));
        assert_eq!(set.nearest_probe_edge("chr1", 1_100, 5), Some(0));
        assert_eq!(set.nearest_probe_edge("chr1", 1_104, 5), Some(4));
        assert_eq!(set.nearest_probe_edge("chr1", 1_110, 10), Some(10));
        assert_eq!(set.nearest_probe_edge("chr1", 1_050, 5), None);
        assert_eq!(set.nearest_probe_edge("chr1", 994, 5), None);
        assert_eq!(set.nearest_probe_edge("chr2", 1_000, 5), None);
    }

    #[test]
    fn restrict_to_a_region() {
        let probes = vec![
            probe("chr1", 1_000, 1_100),
            probe("chr1", 3_000, 3_100),
            probe("chr2", 10, 20),
        ];
        let mut set = TargetSet::new(probes, 0, &contigs()).unwrap();
        set.restrict(&"chr1:1051-3050".parse().unwrap());
        assert_eq!(
            set.regions(),
            [region("chr1", 1_050, 1_100), region("chr1", 3_000, 3_050)]
        );
        assert_eq!(set.probes().len(), 3);

        set.restrict(&"chr2".parse().unwrap());
        assert_eq!(set.regions(), []);
    }
}
