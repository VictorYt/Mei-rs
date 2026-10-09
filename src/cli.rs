//! Command-line interface.

use std::num::NonZeroUsize;
use std::path::PathBuf;

use clap::{ArgAction, Args, Parser, Subcommand, ValueEnum};

use crate::regions::Region;

/// Detect mobile element insertions (Alu, LINE-1, SVA) in targeted capture
/// and exome sequencing data.
#[derive(Debug, Parser)]
#[command(name = "mei-rs", version, propagate_version = true)]
pub struct Cli {
    /// Number of worker threads [default: all available cores]
    #[arg(short, long, global = true, value_name = "N")]
    pub threads: Option<NonZeroUsize>,

    /// Log more details (-v: debug, -vv: trace)
    #[arg(short, long, global = true, action = ArgAction::Count, conflicts_with = "quiet")]
    pub verbose: u8,

    /// Only log errors
    #[arg(short, long, global = true)]
    pub quiet: bool,

    /// Subcommand to run.
    #[command(subcommand)]
    pub command: Command,
}

/// Subcommands.
#[derive(Debug, Subcommand)]
pub enum Command {
    /// Find candidate insertion sites in the capture targets
    Scan(ScanArgs),
}

/// Arguments of `mei-rs scan`.
#[derive(Debug, Args)]
pub struct ScanArgs {
    /// Alignments to scan: BAM or CRAM, coordinate-sorted and indexed
    #[arg(long, value_name = "FILE")]
    pub bam: PathBuf,

    /// Capture targets (BED)
    #[arg(long, value_name = "FILE")]
    pub bed: PathBuf,

    /// Reference genome (FASTA, indexed); required for CRAM
    #[arg(long, value_name = "FILE")]
    pub reference: Option<PathBuf>,

    /// Bases added on each side of every target
    #[arg(long, value_name = "BP", default_value_t = 300)]
    pub padding: u32,

    /// Only scan this region (chr, or chr:start-end, 1-based inclusive)
    #[arg(long, value_name = "REGION")]
    pub region: Option<Region>,

    /// Output file [default: standard output]
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,

    /// Output format
    #[arg(long, value_enum, default_value_t = OutputFormat::Tsv)]
    pub format: OutputFormat,
}

/// Format of the candidate list.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum OutputFormat {
    /// Tab-separated values, one candidate per line, 1-based coordinates
    Tsv,
    /// BED, 0-based half-open coordinates
    Bed,
    /// JSON, with the detailed signals of every candidate
    Json,
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;
    use clap::error::ErrorKind;

    fn parse(args: &[&str]) -> Result<Cli, clap::Error> {
        Cli::try_parse_from(std::iter::once("mei-rs").chain(args.iter().copied()))
    }

    fn scan(args: &[&str]) -> Result<ScanArgs, clap::Error> {
        let base = ["scan", "--bam", "s.bam", "--bed", "t.bed"];
        let all: Vec<&str> = base.iter().chain(args).copied().collect();
        parse(&all).map(|cli| match cli.command {
            Command::Scan(args) => args,
        })
    }

    #[test]
    fn definition_is_valid() {
        Cli::command().debug_assert();
    }

    #[test]
    fn scan_defaults() {
        let args = scan(&[]).unwrap();
        assert_eq!(args.bam, PathBuf::from("s.bam"));
        assert_eq!(args.bed, PathBuf::from("t.bed"));
        assert_eq!(args.reference, None);
        assert_eq!(args.padding, 300);
        assert_eq!(args.region, None);
        assert_eq!(args.output, None);
        assert_eq!(args.format, OutputFormat::Tsv);
    }

    #[test]
    fn scan_all_options() {
        let args = scan(&[
            "--reference",
            "g.fa",
            "--padding",
            "50",
            "--region",
            "chr1:10-20",
            "-o",
            "out.bed",
            "--format",
            "bed",
        ])
        .unwrap();
        assert_eq!(args.reference, Some(PathBuf::from("g.fa")));
        assert_eq!(args.padding, 50);
        assert_eq!(args.region, Some("chr1:10-20".parse().unwrap()));
        assert_eq!(args.output, Some(PathBuf::from("out.bed")));
        assert_eq!(args.format, OutputFormat::Bed);
    }

    #[test]
    fn scan_requires_bam_and_bed() {
        let err = parse(&["scan", "--bed", "t.bed"]).unwrap_err();
        assert_eq!(err.kind(), ErrorKind::MissingRequiredArgument);
        let err = parse(&["scan", "--bam", "s.bam"]).unwrap_err();
        assert_eq!(err.kind(), ErrorKind::MissingRequiredArgument);
    }

    #[test]
    fn invalid_values_are_rejected() {
        assert_eq!(
            scan(&["--format", "vcf"]).unwrap_err().kind(),
            ErrorKind::InvalidValue
        );
        assert_eq!(
            scan(&["--padding", "-1"]).unwrap_err().kind(),
            ErrorKind::UnknownArgument
        );
        assert_eq!(
            scan(&["--padding", "abc"]).unwrap_err().kind(),
            ErrorKind::ValueValidation
        );
        assert_eq!(
            scan(&["--region", "chr1:20-10"]).unwrap_err().kind(),
            ErrorKind::ValueValidation
        );
    }

    #[test]
    fn global_options() {
        let cli = parse(&["-vv", "scan", "--bam", "s.bam", "--bed", "t.bed", "-t", "4"]).unwrap();
        assert_eq!(cli.verbose, 2);
        assert!(!cli.quiet);
        assert_eq!(cli.threads, NonZeroUsize::new(4));

        let cli = parse(&["scan", "--bam", "s.bam", "--bed", "t.bed", "-q"]).unwrap();
        assert!(cli.quiet);
        assert_eq!(cli.threads, None);
    }

    #[test]
    fn verbose_conflicts_with_quiet() {
        let err = parse(&["-v", "-q", "scan", "--bam", "s.bam", "--bed", "t.bed"]).unwrap_err();
        assert_eq!(err.kind(), ErrorKind::ArgumentConflict);
    }

    #[test]
    fn zero_threads_is_rejected() {
        let err = parse(&["-t", "0", "scan", "--bam", "s.bam", "--bed", "t.bed"]).unwrap_err();
        assert_eq!(err.kind(), ErrorKind::ValueValidation);
    }
}
