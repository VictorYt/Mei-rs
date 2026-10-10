//! Loading and consistency checks of the `scan` inputs: alignment file,
//! capture BED, reference and `--region`.
//!
//! Incompatible inputs fail early with a message that says how to fix them;
//! recoverable problems (targets on contigs absent from the alignment, stale
//! index) are logged as warnings.

use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};

use thiserror::Error;
use tracing::warn;

use crate::cli::ScanArgs;
use crate::reader::alignment::{AlignmentError, AlignmentFile};
use crate::regions::Region;
use crate::regions::bed::{self, BedError, Probe};
use crate::regions::targets::{Contig, TargetError, TargetSet};

/// Number of examples listed in messages.
const EXAMPLES: usize = 5;

/// Errors raised while loading the inputs.
#[derive(Debug, Error)]
pub enum InputError {
    /// The alignment file cannot be used.
    #[error(transparent)]
    Alignment(#[from] AlignmentError),
    /// The BED file cannot be read.
    #[error(transparent)]
    Bed(#[from] BedError),
    /// The targets cannot be built.
    #[error(transparent)]
    Targets(#[from] TargetError),
    /// The BED and the alignment name their contigs differently.
    #[error(
        "the contigs of {bed} are named differently from those of {alignment} (e.g. {}); rename the contigs of one of the files",
        .examples.iter().map(|(b, a)| format!("`{b}` in the BED, `{a}` in the alignment")).collect::<Vec<_>>().join(", ")
    )]
    ContigNaming {
        /// Path of the BED file.
        bed: PathBuf,
        /// Path of the alignment file.
        alignment: PathBuf,
        /// BED name and alignment name of a few contigs.
        examples: Vec<(String, String)>,
    },
    /// No contig of the BED is in the alignment header.
    #[error("no contig of {bed} ({}) is in the header of {alignment}", .contigs.join(", "))]
    NoKnownContig {
        /// Path of the BED file.
        bed: PathBuf,
        /// Path of the alignment file.
        alignment: PathBuf,
        /// A few BED contigs.
        contigs: Vec<String>,
    },
    /// The reference has no `.fai` index, or it cannot be read.
    #[error("cannot read the index {path} of the reference (run `samtools faidx`): {source}")]
    ReferenceIndex {
        /// Path of the `.fai` file.
        path: PathBuf,
        /// Underlying I/O error.
        source: std::io::Error,
    },
    /// A target contig is missing from the reference.
    #[error(
        "contig {contig} of {alignment} is missing from the reference {reference}; is it the same genome build?"
    )]
    ReferenceContig {
        /// Contig name.
        contig: String,
        /// Path of the alignment file.
        alignment: PathBuf,
        /// Path of the reference.
        reference: PathBuf,
    },
    /// A target contig has a different length in the reference.
    #[error(
        "contig {contig} has {reference_length} bp in the reference {reference} but {header_length} bp in the header of {alignment}; is it the same genome build?"
    )]
    ReferenceLength {
        /// Contig name.
        contig: String,
        /// Path of the alignment file.
        alignment: PathBuf,
        /// Path of the reference.
        reference: PathBuf,
        /// Length in the alignment header.
        header_length: u64,
        /// Length in the reference.
        reference_length: u64,
    },
    /// The contig of `--region` is not in the alignment header.
    #[error("contig {contig} of --region is not in the header of {alignment}{}",
        .suggestion.as_ref().map(|s| format!(" (did you mean {s}?)")).unwrap_or_default())]
    RegionContig {
        /// Contig name.
        contig: String,
        /// Path of the alignment file.
        alignment: PathBuf,
        /// Contig of the header with the other naming style, if any.
        suggestion: Option<String>,
    },
    /// `--region` holds no target.
    #[error("--region {region} holds no target")]
    EmptyRegion {
        /// The region.
        region: Region,
    },
}

/// The checked inputs of `scan`.
#[derive(Debug)]
pub struct Inputs {
    /// Alignment file.
    pub alignment: AlignmentFile,
    /// Capture targets and regions to scan.
    pub targets: TargetSet,
    /// Sample name.
    pub sample: String,
}

impl Inputs {
    /// Open and check the inputs of `scan`.
    ///
    /// # Errors
    ///
    /// Fails if a file cannot be read or if the files are inconsistent with
    /// each other (see [`InputError`]).
    pub fn load(args: &ScanArgs) -> Result<Self, InputError> {
        let alignment = AlignmentFile::open(&args.bam, args.reference.as_deref())?;
        warn_if_stale_index(&alignment);
        let contigs = alignment.contigs();

        let probes = bed::read(&args.bed)?;
        let probes = check_probes(probes, &contigs, &args.bed, alignment.path())?;
        if let Some(reference) = &args.reference {
            check_reference(reference, &probes, &contigs, alignment.path())?;
        }

        let mut targets = TargetSet::new(probes, args.padding, &contigs)?;
        if let Some(region) = &args.region {
            if !contigs.iter().any(|c| c.name == region.contig) {
                return Err(InputError::RegionContig {
                    contig: region.contig.clone(),
                    alignment: alignment.path().to_path_buf(),
                    suggestion: other_style(&region.contig, &contigs),
                });
            }
            targets.restrict(region);
            if targets.regions().is_empty() {
                return Err(InputError::EmptyRegion {
                    region: region.clone(),
                });
            }
        }

        let sample = sample_name(&alignment);
        Ok(Self {
            alignment,
            targets,
            sample,
        })
    }
}

/// Keep the probes the alignment can be queried on.
///
/// A naming mismatch (`chr1` vs `1`) or a BED without any known contig is an
/// error; probes on other unknown contigs, or starting after the end of their
/// contig, are dropped with a warning.
fn check_probes(
    probes: Vec<Probe>,
    contigs: &[Contig],
    bed: &Path,
    alignment: &Path,
) -> Result<Vec<Probe>, InputError> {
    let lengths: HashMap<&str, u64> = contigs
        .iter()
        .map(|c| (c.name.as_str(), c.length))
        .collect();

    let mut missing: Vec<&str> = Vec::new();
    for probe in &probes {
        if !lengths.contains_key(probe.contig.as_str()) && !missing.contains(&probe.contig.as_str())
        {
            missing.push(&probe.contig);
        }
    }
    let renamed: Vec<(String, String)> = missing
        .iter()
        .filter_map(|&name| other_style(name, contigs).map(|other| (name.to_owned(), other)))
        .take(EXAMPLES)
        .collect();
    if !renamed.is_empty() {
        return Err(InputError::ContigNaming {
            bed: bed.to_path_buf(),
            alignment: alignment.to_path_buf(),
            examples: renamed,
        });
    }
    if missing.len()
        == probes
            .iter()
            .map(|p| &p.contig)
            .collect::<HashSet<_>>()
            .len()
    {
        return Err(InputError::NoKnownContig {
            bed: bed.to_path_buf(),
            alignment: alignment.to_path_buf(),
            contigs: missing
                .iter()
                .take(EXAMPLES)
                .map(|&c| c.to_owned())
                .collect(),
        });
    }

    let missing = examples(missing.iter().map(|&c| c.to_owned()));
    let total = probes.len();
    let (kept, unknown): (Vec<Probe>, Vec<Probe>) = probes
        .into_iter()
        .partition(|p| lengths.contains_key(p.contig.as_str()));
    if !unknown.is_empty() {
        warn!(
            "{} of {total} targets are on contigs absent from {} and are skipped: {}",
            unknown.len(),
            alignment.display(),
            missing
        );
    }
    let (kept, beyond): (Vec<Probe>, Vec<Probe>) = kept
        .into_iter()
        .partition(|p| p.start < lengths[p.contig.as_str()]);
    if !beyond.is_empty() {
        warn!(
            "{} targets start after the end of their contig and are skipped: {}",
            beyond.len(),
            examples(
                beyond
                    .iter()
                    .map(|p| format!("{}:{}-{}", p.contig, p.start, p.end))
            )
        );
    }
    let clamped = kept
        .iter()
        .filter(|p| p.end > lengths[p.contig.as_str()])
        .count();
    if clamped > 0 {
        warn!("{clamped} targets end after the end of their contig and are clamped");
    }
    Ok(kept)
}

/// Check that the reference has the target contigs, with the lengths of the
/// alignment header.
fn check_reference(
    reference: &Path,
    probes: &[Probe],
    contigs: &[Contig],
    alignment: &Path,
) -> Result<(), InputError> {
    let mut fai_path = reference.as_os_str().to_owned();
    fai_path.push(".fai");
    let fai_path = PathBuf::from(fai_path);
    let index =
        noodles::fasta::fai::fs::read(&fai_path).map_err(|source| InputError::ReferenceIndex {
            path: fai_path.clone(),
            source,
        })?;
    let reference_lengths: HashMap<String, u64> = index
        .as_ref()
        .iter()
        .map(|record| (record.name().to_string(), record.length()))
        .collect();

    let used: HashSet<&str> = probes.iter().map(|p| p.contig.as_str()).collect();
    for contig in contigs.iter().filter(|c| used.contains(c.name.as_str())) {
        match reference_lengths.get(&contig.name) {
            None => {
                return Err(InputError::ReferenceContig {
                    contig: contig.name.clone(),
                    alignment: alignment.to_path_buf(),
                    reference: reference.to_path_buf(),
                });
            }
            Some(&length) if length != contig.length => {
                return Err(InputError::ReferenceLength {
                    contig: contig.name.clone(),
                    alignment: alignment.to_path_buf(),
                    reference: reference.to_path_buf(),
                    header_length: contig.length,
                    reference_length: length,
                });
            }
            Some(_) => {}
        }
    }
    Ok(())
}

/// Name of `name` in the other naming style (`chr1` ↔ `1`, `chrM` ↔ `MT`), if
/// a contig of the header has it.
fn other_style(name: &str, contigs: &[Contig]) -> Option<String> {
    let other = match name {
        "chrM" => "MT".to_owned(),
        "MT" => "chrM".to_owned(),
        _ => match name.strip_prefix("chr") {
            Some(stripped) => stripped.to_owned(),
            None => format!("chr{name}"),
        },
    };
    contigs.iter().any(|c| c.name == other).then_some(other)
}

/// Warn if the index is older than the alignment file, as samtools does: the
/// index may describe a previous version of the file.
fn warn_if_stale_index(alignment: &AlignmentFile) {
    let modified = |path: &Path| fs::metadata(path).and_then(|m| m.modified()).ok();
    if let (Some(data), Some(index)) =
        (modified(alignment.path()), modified(alignment.index_path()))
        && index < data
    {
        warn!(
            "the index {} is older than {}; re-run `samtools index` if the file has changed",
            alignment.index_path().display(),
            alignment.path().display()
        );
    }
}

/// Sample name: the `SM` of the read groups, or else the file name.
fn sample_name(alignment: &AlignmentFile) -> String {
    let samples = alignment.samples();
    match samples.as_slice() {
        [sample] => sample.clone(),
        [] => {
            let name = alignment
                .path()
                .file_stem()
                .map_or_else(|| "sample".to_owned(), |s| s.to_string_lossy().into_owned());
            warn!(
                "no sample name (SM) in the read groups of {}; using `{name}`",
                alignment.path().display()
            );
            name
        }
        [first, ..] => {
            warn!(
                "{} holds several samples ({}); using `{first}`",
                alignment.path().display(),
                samples.join(", ")
            );
            first.clone()
        }
    }
}

/// The first few items, comma-separated, with an ellipsis if there are more.
fn examples(items: impl Iterator<Item = String>) -> String {
    let mut items = items.peekable();
    let mut shown: Vec<String> = items.by_ref().take(EXAMPLES).collect();
    if items.peek().is_some() {
        shown.push("...".to_owned());
    }
    shown.join(", ")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{Cli, Command};
    use clap::Parser;

    fn data(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data")
            .join(name)
    }

    /// A temporary directory, removed when dropped.
    struct TempDir(PathBuf);

    impl TempDir {
        fn new(name: &str) -> Self {
            let dir = std::env::temp_dir().join(format!("mei-rs-{name}-{}", std::process::id()));
            fs::create_dir_all(&dir).unwrap();
            Self(dir)
        }

        fn file(&self, name: &str, content: &str) -> PathBuf {
            let path = self.0.join(name);
            fs::write(&path, content).unwrap();
            path
        }
    }

    impl Drop for TempDir {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    fn args(bed: &Path) -> ScanArgs {
        let bam = data("positive.bam");
        let cli = Cli::try_parse_from([
            "mei-rs".as_ref(),
            "scan".as_ref(),
            "--bam".as_ref(),
            bam.as_os_str(),
            "--bed".as_ref(),
            bed.as_os_str(),
        ])
        .unwrap();
        let Command::Scan(args) = cli.command;
        args
    }

    #[test]
    fn loads_the_fixtures() {
        for (bam, reference) in [
            ("positive.bam", None),
            ("positive.bam", Some(data("reference.fa"))),
            ("positive.cram", Some(data("reference.fa"))),
        ] {
            let mut args = args(&data("targets.bed"));
            args.bam = data(bam);
            args.reference = reference;
            let inputs = Inputs::load(&args).unwrap();
            assert_eq!(inputs.sample, "positive");
            assert_eq!(inputs.targets.probes().len(), 3);
            assert_eq!(inputs.targets.probe_bases(), 300 + 200 + 300);
            assert_eq!(inputs.targets.scanned_bases(), 900 + 800 + 900);
        }
    }

    #[test]
    fn contig_naming_mismatch() {
        let dir = TempDir::new("naming");
        let bed = dir.file("t.bed", "1\t9850\t10150\nchrM\t1\t10\n");
        let err = Inputs::load(&args(&bed)).unwrap_err();
        assert!(matches!(err, InputError::ContigNaming { .. }));
        assert!(
            err.to_string()
                .contains("(e.g. `1` in the BED, `chr1` in the alignment)"),
            "{err}"
        );
    }

    #[test]
    fn no_known_contig() {
        let dir = TempDir::new("unknown");
        let bed = dir.file("t.bed", "chrZ\t1\t10\n");
        let err = Inputs::load(&args(&bed)).unwrap_err();
        assert!(matches!(err, InputError::NoKnownContig { .. }));
        assert!(err.to_string().starts_with("no contig of "), "{err}");
    }

    #[test]
    fn unknown_and_out_of_range_targets_are_skipped() {
        let dir = TempDir::new("skipped");
        let bed = dir.file(
            "t.bed",
            "chr1\t9850\t10150\nchrUn_x\t1\t10\nchr1\t25000\t25100\nchr2\t9900\t10100\n",
        );
        let inputs = Inputs::load(&args(&bed)).unwrap();
        let probes: Vec<_> = inputs
            .targets
            .probes()
            .iter()
            .map(|p| (p.contig.as_str(), p.start, p.end))
            .collect();
        assert_eq!(probes, [("chr1", 9850, 10150), ("chr2", 9900, 10000)]);
    }

    #[test]
    fn reference_mismatch() {
        let dir = TempDir::new("reference");
        let mut args = args(&data("targets.bed"));

        let fasta = dir.file("other.fa", ">chr1\nACGT\n>chr2\nACGT\n");
        dir.file("other.fa.fai", "chr1\t4\t6\t4\t5\nchr2\t4\t17\t4\t5\n");
        args.reference = Some(fasta);
        let err = Inputs::load(&args).unwrap_err();
        assert!(
            matches!(
                err,
                InputError::ReferenceLength {
                    header_length: 20_000,
                    reference_length: 4,
                    ..
                }
            ),
            "{err}"
        );

        let fasta = dir.file("partial.fa", ">chr2\nACGT\n");
        dir.file("partial.fa.fai", "chr2\t4\t6\t4\t5\n");
        args.reference = Some(fasta);
        let err = Inputs::load(&args).unwrap_err();
        assert!(matches!(err, InputError::ReferenceContig { ref contig, .. } if contig == "chr1"));

        args.reference = Some(dir.file("noindex.fa", ">chr1\nACGT\n"));
        let err = Inputs::load(&args).unwrap_err();
        assert!(matches!(err, InputError::ReferenceIndex { .. }));
    }

    #[test]
    fn region() {
        let mut args = args(&data("targets.bed"));
        args.region = Some("chr1:10001-20000".parse().unwrap());
        let inputs = Inputs::load(&args).unwrap();
        assert_eq!(inputs.targets.regions().len(), 1);
        assert_eq!(inputs.targets.scanned_bases(), 10_450 - 10_000);

        args.region = Some("1:1-100".parse().unwrap());
        let err = Inputs::load(&args).unwrap_err();
        assert_eq!(
            err.to_string(),
            format!(
                "contig 1 of --region is not in the header of {} (did you mean chr1?)",
                data("positive.bam").display()
            )
        );

        args.region = Some("chr1:1-100".parse().unwrap());
        let err = Inputs::load(&args).unwrap_err();
        assert_eq!(err.to_string(), "--region chr1:1-100 holds no target");
    }

    #[test]
    fn other_naming_style() {
        let contigs: Vec<Contig> = ["chr1", "MT", "2"]
            .iter()
            .map(|name| Contig {
                name: (*name).to_owned(),
                length: 1,
            })
            .collect();
        assert_eq!(other_style("1", &contigs), Some("chr1".to_owned()));
        assert_eq!(other_style("chrM", &contigs), Some("MT".to_owned()));
        assert_eq!(other_style("chr2", &contigs), Some("2".to_owned()));
        assert_eq!(other_style("chr3", &contigs), None);
    }

    #[test]
    fn examples_are_truncated() {
        let items = (1..=7).map(|i| i.to_string());
        assert_eq!(examples(items), "1, 2, 3, 4, 5, ...");
        assert_eq!(examples(["a".to_owned()].into_iter()), "a");
    }
}
