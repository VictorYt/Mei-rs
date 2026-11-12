//! Ploidy of the sex chromosomes (#127).
//!
//! A male has one X and one Y: outside their pseudoautosomal regions (PAR),
//! an insertion on them is haploid, `0` or `1`. The sex is given
//! (`--sex female|male`) or inferred from the reads: the reads per targeted
//! base of the X against those of the autosomes, about 1 for a female and
//! 0.5 for a male ([`infer`]).
//!
//! The PAR are known for `GRCh37` and `GRCh38` only, recognised by the length
//! of the X (or the Y): on another assembly the whole X and Y are taken as
//! haploid in a male. Contigs are named `X`, `Y`, `1`..`22`, with or without
//! `chr`. Coordinates are 0-based.

/// Sex of the sample.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Sex {
    /// Two X.
    Female,
    /// One X and one Y.
    Male,
}

impl Sex {
    /// `female` or `male`.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::Female => "female",
            Self::Male => "male",
        }
    }
}

/// Kind of chromosome of a contig.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Chromosome {
    /// `1` to `22`.
    Autosome,
    /// `X`.
    X,
    /// `Y`.
    Y,
}

/// Kind of chromosome of a contig name, if it is one of the 24.
#[must_use]
pub fn chromosome(contig: &str) -> Option<Chromosome> {
    let name = contig.strip_prefix("chr").unwrap_or(contig);
    match name {
        "X" => Some(Chromosome::X),
        "Y" => Some(Chromosome::Y),
        _ => name
            .parse::<u8>()
            .ok()
            .filter(|n| (1..=22).contains(n) && !name.starts_with('0'))
            .map(|_| Chromosome::Autosome),
    }
}

/// Pseudoautosomal regions of an assembly: contig length, then the regions
/// (0-based, half-open).
type ParTable = [(u64, [(u64, u64); 2]); 4];

/// PAR1 and PAR2 of the X and the Y of `GRCh38` and `GRCh37` (Genome Reference
/// Consortium), by contig length.
const PAR: ParTable = [
    // `GRCh38` X, Y.
    (
        156_040_895,
        [(10_000, 2_781_479), (155_701_382, 156_030_895)],
    ),
    (57_227_415, [(10_000, 2_781_479), (56_887_902, 57_217_415)]),
    // `GRCh37` X, Y.
    (
        155_270_560,
        [(60_000, 2_699_520), (154_931_043, 155_260_560)],
    ),
    (59_373_566, [(10_000, 2_649_520), (59_034_049, 59_363_566)]),
];

/// Whether the PAR of a sex chromosome of this length are known.
#[must_use]
pub fn known_par(length: u64) -> bool {
    PAR.iter().any(|(l, _)| *l == length)
}

/// Whether position `pos` of a contig (name, length) is haploid in a sample
/// of sex `sex`: on the X or the Y of a male, outside the PAR.
#[must_use]
pub fn is_haploid(sex: Option<Sex>, contig: &str, length: u64, pos: u64) -> bool {
    if sex != Some(Sex::Male) || chromosome(contig).is_none_or(|c| c == Chromosome::Autosome) {
        return false;
    }
    let in_par = PAR
        .iter()
        .filter(|(l, _)| *l == length)
        .flat_map(|(_, regions)| regions)
        .any(|(start, end)| (*start..*end).contains(&pos));
    !in_par
}

/// Fewest reads on the targets of the autosomes, and on those of the X, to
/// infer the sex.
pub const MIN_SEX_READS: u64 = 200;

/// Reads per targeted base of the X against the autosomes below which the
/// sample is male. One X does not give half the reads of the autosomes on
/// every capture: on 36 1000 Genomes exomes, the 22 males had 0.61 to 0.84
/// (4 of them above 0.8) and the 14 females 1.12 to 1.63.
pub const MALE_MAX_X_RATIO: f64 = 0.95;

/// Sex from the reads and bases of the targets of the autosomes and of the
/// X: male when the X has fewer than [`MALE_MAX_X_RATIO`] of the reads per
/// base of the autosomes, female otherwise; `None` with fewer than
/// [`MIN_SEX_READS`] reads on either. Returns the ratio as well.
#[must_use]
pub fn infer(autosomes: (u64, u64), x: (u64, u64)) -> (Option<Sex>, Option<f64>) {
    let ((a_reads, a_bases), (x_reads, x_bases)) = (autosomes, x);
    if a_reads < MIN_SEX_READS || x_reads < MIN_SEX_READS || a_bases == 0 || x_bases == 0 {
        return (None, None);
    }
    #[allow(clippy::cast_precision_loss)]
    let ratio = (x_reads as f64 / x_bases as f64) / (a_reads as f64 / a_bases as f64);
    let sex = if ratio < MALE_MAX_X_RATIO {
        Sex::Male
    } else {
        Sex::Female
    };
    (Some(sex), Some(ratio))
}

#[cfg(test)]
mod tests {
    use super::*;

    const GRCH38_X: u64 = 156_040_895;
    const GRCH37_Y: u64 = 59_373_566;

    #[test]
    fn chromosomes() {
        assert_eq!(chromosome("chrX"), Some(Chromosome::X));
        assert_eq!(chromosome("Y"), Some(Chromosome::Y));
        assert_eq!(chromosome("chr22"), Some(Chromosome::Autosome));
        assert_eq!(chromosome("1"), Some(Chromosome::Autosome));
        for other in [
            "chr23",
            "chr0",
            "chr01",
            "chrM",
            "chr1_KI270706v1_random",
            "MT",
        ] {
            assert_eq!(chromosome(other), None, "{other}");
        }
    }

    #[test]
    fn haploid_outside_the_par_of_a_male() {
        // X of GRCh38: PAR1 ends at 2,781,479.
        assert!(!is_haploid(Some(Sex::Male), "chrX", GRCH38_X, 2_781_478));
        assert!(is_haploid(Some(Sex::Male), "chrX", GRCH38_X, 2_781_479));
        assert!(!is_haploid(Some(Sex::Male), "chrX", GRCH38_X, 155_701_382));
        // Y of GRCh37: PAR1 from 10,000.
        assert!(is_haploid(Some(Sex::Male), "Y", GRCH37_Y, 9_999));
        assert!(!is_haploid(Some(Sex::Male), "Y", GRCH37_Y, 10_000));
        // A female, an autosome, an unknown sex: diploid.
        assert!(!is_haploid(Some(Sex::Female), "chrX", GRCH38_X, 50_000_000));
        assert!(!is_haploid(
            Some(Sex::Male),
            "chr7",
            159_345_973,
            50_000_000
        ));
        assert!(!is_haploid(None, "chrX", GRCH38_X, 50_000_000));
        // Another assembly: no PAR known, the whole X is haploid.
        assert!(!known_par(20_000));
        assert!(is_haploid(Some(Sex::Male), "chrX", 20_000, 15_000));
    }

    #[test]
    fn sex_from_the_depth_of_the_x() {
        // Same targeted bases; half the reads per base on the X: male.
        assert_eq!(infer((10_000, 1_000), (500, 100)).0, Some(Sex::Male));
        assert_eq!(infer((10_000, 1_000), (1_000, 100)).0, Some(Sex::Female));
        // The highest ratio of a male exome seen.
        assert_eq!(infer((10_000, 1_000), (840, 100)).0, Some(Sex::Male));
        let (_, ratio) = infer((10_000, 1_000), (500, 100));
        assert!((ratio.unwrap() - 0.5).abs() < 1e-9);
        // Too few reads: unknown.
        assert_eq!(infer((10_000, 1_000), (150, 100)), (None, None));
        assert_eq!(infer((0, 0), (500, 100)), (None, None));
    }
}
