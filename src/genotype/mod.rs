//! Phase 5: genotype of the candidates from the reads of each allele at
//! their junctions (see [`crate::reader::alleles`]), and their quality
//! score ([`score`]).
//!
//! At each junction, the alternative reads among the reads counted follow a
//! binomial law of parameter [`ERROR`] for `0/0`, [`HET_FRACTION`] for `0/1`
//! and `1 - ERROR` for `1/1`. With two junctions, the log-likelihoods are
//! averaged: a reference read usually spans both junctions (they are a TSD
//! apart) and would count twice. The genotype is the most likely one, and
//! its quality (GQ) the phred-scaled probability that it is wrong given the
//! reads alone (the three genotypes, equal priors), at most [`MAX_QUALITY`].
//!
//! A candidate that passed the filters is an insertion in the sample: its
//! genotype is `0/1` or `1/1` (never `0/0`), with a prior that favours the
//! heterozygous genotype ([`HET_PRIOR`]). `0/0` stays possible for the
//! candidates that did not pass (`--keep-filtered`). The prior only chooses
//! the call: GQ measures the reads, so a call the reads do not support (no
//! alternative read at a mate-only insertion) has a GQ of 0 (#117).

pub mod score;

/// Fraction of alternative reads at a site without the insertion (reads
/// clipped there for another reason), and of reference reads at a
/// homozygous one.
pub const ERROR: f64 = 0.02;

/// Expected fraction of alternative reads at a heterozygous insertion:
/// below one half, as alternative reads are under-counted (capture probes
/// of the reference, clips of `overhang` bases of good quality). The median
/// fraction of the true heterozygous `PASS` insertions of the 7 HGSVC3
/// exomes (#26) is 0.28; with 0.3, 76 of their 91 genotypes are right
/// (68 with 0.5).
pub const HET_FRACTION: f64 = 0.3;

/// Prior probability of a heterozygous genotype among the carriers of an
/// insertion: 77% of the HGSVC3 genotypes of the seven exomes with an
/// insertion are heterozygous (4,471 + 4,367 haplotype pairs `0|1` / `1|0`
/// against 2,654 `1|1`). Used for the candidates that passed.
pub const HET_PRIOR: f64 = 0.77;

/// Smallest fraction of alternative reads among the reads at the site
/// (local depth) for a homozygous call when no reference read is counted.
///
/// Reference reads are counted only when they span the site by 15 bases on
/// each side; next to a probe edge, the reads of the reference allele end
/// at the site and none is counted, so a heterozygous insertion looks
/// homozygous. Its alternative reads are then a smaller part of the reads
/// at the site: among the true calls without reference read of the 18
/// tuning exomes, median 0.18 for the 27 heterozygous ones and 0.38 for the
/// 63 homozygous ones (#91). The value is the best of 0.15 to 0.3 there
/// (235 of 263 genotypes right, 229 before).
pub const HOM_MIN_ALT_FRACTION: f64 = 0.2;

/// Highest genotype quality.
pub const MAX_QUALITY: u32 = 99;

/// Reads of each allele at one junction.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct AlleleCounts {
    /// Reads holding the insertion.
    pub alternative: u32,
    /// Reads holding the site without the insertion.
    pub reference: u32,
}

/// A genotype call.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Call {
    /// No insertion.
    HomRef,
    /// Heterozygous insertion.
    Het,
    /// Homozygous insertion.
    HomAlt,
}

impl Call {
    /// VCF notation: `0/0`, `0/1` or `1/1`.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::HomRef => "0/0",
            Self::Het => "0/1",
            Self::HomAlt => "1/1",
        }
    }
}

/// Genotype of a candidate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Genotype {
    /// Most likely genotype.
    pub call: Call,
    /// Genotype quality (phred).
    pub quality: u32,
    /// Alternative reads, summed over the junctions.
    pub alternative: u32,
    /// Reference reads, summed over the junctions (a read spanning the site
    /// counts at each junction).
    pub reference: u32,
    /// log10 likelihoods of `0/0`, `0/1` and `1/1`.
    pub log_likelihoods: [f64; 3],
}

/// Genotype from the reads of each allele at the junctions (the reference
/// reads of the site at each junction); `None` without any read.
///
/// `present` is true for a candidate that passed the filters: the genotype
/// is then `0/1` or `1/1`, with [`HET_PRIOR`]; a homozygous call without
/// reference read is heterozygous when the alternative reads are less than
/// [`HOM_MIN_ALT_FRACTION`] of `depth`, the reads at the site. Otherwise the
/// three genotypes have equal priors.
#[must_use]
pub fn genotype(junctions: &[AlleleCounts], present: bool, depth: Option<u64>) -> Option<Genotype> {
    let alternative: u32 = junctions.iter().map(|c| c.alternative).sum();
    let reference: u32 = junctions.iter().map(|c| c.reference).sum();
    if alternative + reference == 0 {
        return None;
    }
    #[allow(clippy::cast_precision_loss)]
    let n = junctions.len() as f64;
    let mut log_likelihoods = [0.0; 3];
    for (ll, p) in log_likelihoods
        .iter_mut()
        .zip([ERROR, HET_FRACTION, 1.0 - ERROR])
    {
        *ll = junctions
            .iter()
            .map(|c| {
                f64::from(c.alternative) * p.log10() + f64::from(c.reference) * (1.0 - p).log10()
            })
            .sum::<f64>()
            / n;
    }
    let calls = [Call::HomRef, Call::Het, Call::HomAlt];
    // log10 prior of each genotype; `0/0` is excluded for a passed candidate.
    let priors = if present {
        [
            f64::NEG_INFINITY,
            HET_PRIOR.log10(),
            (1.0 - HET_PRIOR).log10(),
        ]
    } else {
        [0.0; 3]
    };
    let posteriors: Vec<f64> = log_likelihoods
        .iter()
        .zip(priors)
        .map(|(ll, prior)| ll + prior)
        .collect();
    let mut best = (0..3)
        .max_by(|&a, &b| posteriors[a].total_cmp(&posteriors[b]))
        .unwrap_or(0);
    #[allow(clippy::cast_precision_loss)]
    if present
        && best == 2
        && reference == 0
        && depth.is_some_and(|d| d > 0 && f64::from(alternative) < HOM_MIN_ALT_FRACTION * d as f64)
    {
        best = 1;
    }
    // Probability of the other genotypes from the reads alone, `0/0` included.
    let total: f64 = log_likelihoods
        .iter()
        .map(|ll| 10f64.powf(ll - log_likelihoods[best]))
        .sum();
    let wrong = 1.0 - 1.0 / total;
    Some(Genotype {
        call: calls[best],
        quality: phred(wrong).min(MAX_QUALITY),
        alternative,
        reference,
        log_likelihoods,
    })
}

/// Phred scale of a probability, rounded (`MAX_QUALITY` for 0).
#[must_use]
pub fn phred(probability: f64) -> u32 {
    if probability <= 0.0 {
        return MAX_QUALITY;
    }
    let q = (-10.0 * probability.log10()).round();
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    if q >= f64::from(MAX_QUALITY) {
        MAX_QUALITY
    } else {
        q.max(0.0) as u32
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn counts(alternative: u32, reference: u32) -> AlleleCounts {
        AlleleCounts {
            alternative,
            reference,
        }
    }

    #[test]
    fn calls() {
        let g = genotype(&[counts(10, 12), counts(11, 12)], false, None).unwrap();
        assert_eq!((g.call, g.alternative, g.reference), (Call::Het, 21, 24));
        assert_eq!(g.quality, MAX_QUALITY);
        assert_eq!(
            genotype(&[counts(15, 0), counts(14, 1)], false, None)
                .unwrap()
                .call,
            Call::HomAlt
        );
        assert_eq!(
            genotype(&[counts(0, 30)], false, None).unwrap().call,
            Call::HomRef
        );
        assert_eq!(genotype(&[counts(0, 0)], false, None), None);
    }

    #[test]
    fn a_passed_candidate_is_never_homozygous_reference() {
        // Few alternative reads among many reference reads: `0/0` without the
        // filters, a heterozygous insertion for a candidate that passed.
        assert_eq!(
            genotype(&[counts(1, 56)], false, None).unwrap().call,
            Call::HomRef
        );
        let g = genotype(&[counts(1, 56)], true, None).unwrap();
        assert_eq!(g.call, Call::Het);
        assert!(g.quality <= MAX_QUALITY);
        // Without reference reads, the homozygous genotype is still called.
        assert_eq!(
            genotype(&[counts(15, 0)], true, None).unwrap().call,
            Call::HomAlt
        );
    }

    #[test]
    fn the_quality_of_a_passed_candidate_measures_the_reads() {
        // No alternative read, or one of 57: the call stays `0/1` but the
        // reads do not support it.
        for reads in [counts(0, 30), counts(1, 56)] {
            let g = genotype(&[reads], true, None).unwrap();
            assert_eq!((g.call, g.quality), (Call::Het, 0), "{reads:?}");
        }
        // Alternative and reference reads at a heterozygous fraction do.
        let g = genotype(&[counts(10, 20)], true, None).unwrap();
        assert_eq!((g.call, g.quality), (Call::Het, 88));
        // Huge counts do not overflow into NaN.
        let g = genotype(&[counts(0, 100_000)], true, None).unwrap();
        assert_eq!((g.call, g.quality), (Call::Het, 0));
    }

    #[test]
    fn missing_reference_reads_next_to_a_probe_edge() {
        // 9 alternative reads, none of the reference, 68 reads at the site:
        // the reference reads end at the site, the insertion is heterozygous.
        let g = genotype(&[counts(9, 0)], true, Some(68)).unwrap();
        assert_eq!(g.call, Call::Het);
        // The reads do not tell it: the quality stays low.
        assert!(g.quality < 10, "{g:?}");
        // Most reads at the site are alternative: homozygous.
        assert_eq!(
            genotype(&[counts(9, 0)], true, Some(28)).unwrap().call,
            Call::HomAlt
        );
        // Without depth, or for a filtered candidate, the reads decide.
        assert_eq!(
            genotype(&[counts(9, 0)], true, None).unwrap().call,
            Call::HomAlt
        );
        assert_eq!(
            genotype(&[counts(9, 0)], false, Some(68)).unwrap().call,
            Call::HomAlt
        );
    }

    #[test]
    fn the_prior_favours_heterozygous_insertions() {
        // 4 alternative reads and 1 reference read: `1/1` from the reads
        // alone, `0/1` once the prior of a passed candidate is applied.
        assert_eq!(
            genotype(&[counts(4, 1)], false, None).unwrap().call,
            Call::HomAlt
        );
        assert_eq!(
            genotype(&[counts(4, 1)], true, None).unwrap().call,
            Call::Het
        );
    }

    #[test]
    fn likelihoods_are_averaged_over_the_junctions() {
        // A read spanning the site counts at both junctions: the same counts
        // at two junctions give the genotype and quality of one.
        let one = genotype(&[counts(3, 4)], false, None).unwrap();
        let two = genotype(&[counts(3, 4), counts(3, 4)], false, None).unwrap();
        assert!(one.quality < MAX_QUALITY);
        assert_eq!((two.call, two.quality), (one.call, one.quality));
    }

    #[test]
    fn homozygous_at_the_minimum_alternative_fraction() {
        // Exactly `HOM_MIN_ALT_FRACTION` of the reads at the site are
        // alternative: still homozygous.
        assert_eq!(
            genotype(&[counts(2, 0)], true, Some(10)).unwrap().call,
            Call::HomAlt
        );
        assert_eq!(
            genotype(&[counts(2, 0)], true, Some(11)).unwrap().call,
            Call::Het
        );
    }

    #[test]
    fn few_reads_give_a_low_quality() {
        let g = genotype(&[counts(1, 0)], false, None).unwrap();
        assert_eq!(g.call, Call::HomAlt);
        assert!(g.quality < 10, "{g:?}");
    }

    #[test]
    fn phred_scale() {
        assert_eq!(phred(0.1), 10);
        assert_eq!(phred(0.001), 30);
        assert_eq!(phred(0.0), MAX_QUALITY);
        assert_eq!(phred(1.0), 0);
    }
}
