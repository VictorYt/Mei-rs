//! Quality score (QUAL) of a candidate: the phred-scaled probability that
//! it is not a mobile element insertion, from its evidence, with a naive
//! Bayes model.
//!
//! Each piece of evidence multiplies the prior odds of a true insertion by
//! its likelihood ratio when present, or when absent. The prior and the
//! ratios come from the `PASS` candidates of the 7 HGSVC3 exomes (#26): 91
//! true and 46 false, with one read added to each count. Leaving one sample
//! out, the score ranks a true candidate above a false one 88% of the time
//! (AUC 0.876; 0.893 on all samples), and a score of 5 or more keeps 73 of
//! the 91 true candidates and 10 of the 46 false ones. The evidence is not
//! independent, so the score ranks the candidates better than it estimates
//! the probability.

use crate::cluster::Candidate;

/// log10 prior odds of a true insertion among the `PASS` candidates.
pub const PRIOR: f64 = 0.292;

/// A piece of evidence and its log10 likelihood ratios (true against false
/// insertions) when present and when absent.
struct Evidence {
    present: fn(&Candidate) -> bool,
    log_ratios: (f64, f64),
}

/// Minimum identity of the typed clips.
pub const MIN_IDENTITY: f64 = 0.98;

/// Minimum fraction of distinct fragment starts among the supporting reads
/// (PCR duplicates).
pub const MIN_DISTINCT_STARTS: f64 = 0.7;

/// Minimum fraction of MEI mates among the checked mates.
pub const MIN_MEI_MATES: f64 = 0.6;

/// Minimum number of supporting reads per read of local depth.
pub const MIN_SIGNAL_RATIO: f64 = 0.2;

/// The evidence (true / false `PASS` candidates of HGSVC3 with it).
const EVIDENCE: [Evidence; 7] = [
    // Signals on both sides (37/91, 23/46).
    Evidence {
        present: |c| c.left_junction.is_some() && c.right_junction.is_some(),
        log_ratios: (-0.088, 0.073),
    },
    // A poly(A/T) tail (43/91, 3/46).
    Evidence {
        present: |c| c.polya_clips + c.polyt_clips > 0,
        log_ratios: (0.754, -0.240),
    },
    // A TSD (12/91, 5/46).
    Evidence {
        present: |c| c.tsd_length().is_some(),
        log_ratios: (0.049, -0.007),
    },
    // Typed clips of high identity (80/91, 22/46).
    Evidence {
        present: |c| {
            c.typing
                .as_ref()
                .is_some_and(|t| t.identity >= MIN_IDENTITY)
        },
        log_ratios: (0.260, -0.606),
    },
    // Distinct fragment starts (49/91, 18/46).
    Evidence {
        present: |c| {
            #[allow(clippy::cast_precision_loss)]
            let fraction = c.distinct_starts as f64 / c.support.max(1) as f64;
            fraction >= MIN_DISTINCT_STARTS
        },
        log_ratios: (0.133, -0.116),
    },
    // Mostly MEI mates (56/91, 8/46).
    Evidence {
        present: |c| {
            #[allow(clippy::cast_precision_loss)]
            let fraction = c.mei_mates as f64 / c.checked_mates.max(1) as f64;
            c.checked_mates > 0 && fraction >= MIN_MEI_MATES
        },
        log_ratios: (0.514, -0.322),
    },
    // Many supporting reads for the local depth (88/91, 34/46).
    Evidence {
        present: |c| c.signal_ratio().is_some_and(|r| r >= MIN_SIGNAL_RATIO),
        log_ratios: (0.118, -0.799),
    },
];

/// Quality score of a candidate, at most [`super::MAX_QUALITY`].
#[must_use]
pub fn quality(candidate: &Candidate) -> u32 {
    let log_odds = PRIOR
        + EVIDENCE
            .iter()
            .map(|e| {
                if (e.present)(candidate) {
                    e.log_ratios.0
                } else {
                    e.log_ratios.1
                }
            })
            .sum::<f64>();
    super::phred(1.0 / (1.0 + 10f64.powf(log_odds)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cluster::tests::clip;
    use crate::cluster::{ClusterOptions, cluster};
    use crate::reader::signals::Side;

    #[test]
    fn evidence_raises_the_score() {
        let tailed: Vec<_> = (0..4)
            .map(|i| clip(1_000, Side::Left, &format!("r{i}"), Some(b'A')))
            .collect();
        let bare: Vec<_> = (0..4)
            .map(|i| clip(1_000, Side::Left, &format!("r{i}"), None))
            .collect();
        let options = ClusterOptions {
            require_tail: false,
            ..ClusterOptions::default()
        };
        let mut with_tail = cluster(tailed, &options).remove(0);
        let without = cluster(bare, &options).remove(0);
        assert!(quality(&with_tail) > quality(&without));
        // MEI mates raise it further.
        let before = quality(&with_tail);
        with_tail.checked_mates = 10;
        with_tail.mei_mates = 8;
        assert!(quality(&with_tail) > before);
    }
}
