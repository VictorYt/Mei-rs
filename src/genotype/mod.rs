//! Phase 5: genotype of the candidates from the reads of each allele at
//! their junctions (see [`crate::reader::alleles`]), and their quality
//! score ([`score`]).
//!
//! At each junction, the alternative reads among the reads counted follow a
//! binomial law of parameter [`ERROR`] for `0/0`, [`HET_FRACTION`] for `0/1`
//! and `1 - ERROR` for `1/1`. With two junctions, the log-likelihoods are
//! averaged: a reference read usually spans both junctions (they are a TSD
//! apart) and would count twice. The genotype is the most likely one, and
//! its quality (GQ) the phred-scaled probability that it is wrong, with
//! equal priors, at most [`MAX_QUALITY`].
//!
//! A candidate that passed the filters is an insertion in the sample: its
//! genotype is `0/1` or `1/1` (never `0/0`), with a prior that favours the
//! heterozygous genotype ([`HET_PRIOR`]). `0/0` stays possible for the
//! candidates that did not pass (`--keep-filtered`).

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
/// is then `0/1` or `1/1`, with [`HET_PRIOR`]. Otherwise the three genotypes
/// have equal priors.
#[must_use]
pub fn genotype(junctions: &[AlleleCounts], present: bool) -> Option<Genotype> {
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
    let best = (0..3)
        .max_by(|&a, &b| posteriors[a].total_cmp(&posteriors[b]))
        .unwrap_or(0);
    // Posterior of the other genotypes.
    let total: f64 = posteriors
        .iter()
        .map(|ll| 10f64.powf(ll - posteriors[best]))
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
        let g = genotype(&[counts(10, 12), counts(11, 12)], false).unwrap();
        assert_eq!((g.call, g.alternative, g.reference), (Call::Het, 21, 24));
        assert_eq!(g.quality, MAX_QUALITY);
        assert_eq!(
            genotype(&[counts(15, 0), counts(14, 1)], false)
                .unwrap()
                .call,
            Call::HomAlt
        );
        assert_eq!(
            genotype(&[counts(0, 30)], false).unwrap().call,
            Call::HomRef
        );
        assert_eq!(genotype(&[counts(0, 0)], false), None);
    }

    #[test]
    fn a_passed_candidate_is_never_homozygous_reference() {
        // Few alternative reads among many reference reads: `0/0` without the
        // filters, a heterozygous insertion for a candidate that passed.
        assert_eq!(
            genotype(&[counts(1, 56)], false).unwrap().call,
            Call::HomRef
        );
        let g = genotype(&[counts(1, 56)], true).unwrap();
        assert_eq!(g.call, Call::Het);
        assert!(g.quality <= MAX_QUALITY);
        // Without reference reads, the homozygous genotype is still called.
        assert_eq!(genotype(&[counts(15, 0)], true).unwrap().call, Call::HomAlt);
    }

    #[test]
    fn few_reads_give_a_low_quality() {
        let g = genotype(&[counts(1, 0)], false).unwrap();
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
