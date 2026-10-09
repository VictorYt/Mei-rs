//! Command-line interface.

use std::num::NonZeroUsize;
use std::path::PathBuf;

use clap::{ArgAction, Args, Parser, Subcommand, ValueEnum};

use crate::reader::extract::ExtractOptions;
use crate::reader::filter::ReadFilter;
use crate::reader::signals::ClipOptions;
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
    #[arg(long, value_name = "FILE", help_heading = "Input")]
    pub bam: PathBuf,

    /// Capture targets (BED)
    #[arg(long, value_name = "FILE", help_heading = "Input")]
    pub bed: PathBuf,

    /// Reference genome (FASTA, indexed); required for CRAM
    #[arg(long, value_name = "FILE", help_heading = "Input")]
    pub reference: Option<PathBuf>,

    /// Bases added on each side of every target
    #[arg(long, value_name = "BP", default_value_t = 300, help_heading = "Input")]
    pub padding: u32,

    /// Only scan this region (chr, or chr:start-end, 1-based inclusive)
    #[arg(long, value_name = "REGION", help_heading = "Input")]
    pub region: Option<Region>,

    /// Minimum mapping quality of an anchored read
    #[arg(
        long,
        value_name = "Q",
        default_value_t = 20,
        help_heading = "Read filtering"
    )]
    pub min_mapq: u8,

    /// Keep reads flagged as duplicates (files without duplicate marking)
    #[arg(long, help_heading = "Read filtering")]
    pub include_duplicates: bool,

    /// Minimum length of a soft-clip
    #[arg(long, value_name = "BP", default_value_t = 20,
        value_parser = clap::value_parser!(u32).range(1..), help_heading = "Signal extraction")]
    pub min_clip_len: u32,

    /// Minimum median base quality of a soft-clip
    #[arg(
        long,
        value_name = "Q",
        default_value_t = 20,
        help_heading = "Signal extraction"
    )]
    pub min_clip_quality: u8,

    /// Minimum length of a poly(A/T) tail at the junction end of a soft-clip
    #[arg(long, value_name = "BP", default_value_t = 10,
        value_parser = clap::value_parser!(u32).range(1..), help_heading = "Signal extraction")]
    pub polya_min_len: u32,

    /// Minimum fraction of A (or T) in a poly(A/T) tail
    #[arg(long, value_name = "FRACTION", default_value_t = 0.8,
        value_parser = parse_fraction, help_heading = "Signal extraction")]
    pub polya_min_frac: f64,

    /// Output file [default: standard output]
    #[arg(short, long, value_name = "FILE", help_heading = "Output")]
    pub output: Option<PathBuf>,

    /// Output format
    #[arg(long, value_enum, default_value_t = OutputFormat::Tsv, help_heading = "Output")]
    pub format: OutputFormat,
}

impl ScanArgs {
    /// Signal extraction options.
    #[must_use]
    pub fn extract_options(&self) -> ExtractOptions {
        ExtractOptions {
            filter: ReadFilter {
                min_mapq: self.min_mapq,
                include_duplicates: self.include_duplicates,
            },
            clip: ClipOptions {
                min_length: self.min_clip_len,
                min_quality: self.min_clip_quality,
                polya_min_length: self.polya_min_len,
                polya_min_fraction: self.polya_min_frac,
            },
        }
    }
}

/// Parse a fraction in (0, 1].
fn parse_fraction(s: &str) -> Result<f64, String> {
    let value: f64 = s.parse().map_err(|_| format!("`{s}` is not a number"))?;
    if value > 0.0 && value <= 1.0 {
        Ok(value)
    } else {
        Err(format!("{value} is not in (0, 1]"))
    }
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
        assert_eq!(args.extract_options(), ExtractOptions::default());
    }

    #[test]
    fn extraction_options() {
        let args = scan(&[
            "--min-mapq",
            "30",
            "--include-duplicates",
            "--min-clip-len",
            "15",
            "--min-clip-quality",
            "10",
            "--polya-min-len",
            "8",
            "--polya-min-frac",
            "0.9",
        ])
        .unwrap();
        let options = args.extract_options();
        assert_eq!(options.filter.min_mapq, 30);
        assert!(options.filter.include_duplicates);
        assert_eq!(options.clip.min_length, 15);
        assert_eq!(options.clip.min_quality, 10);
        assert_eq!(options.clip.polya_min_length, 8);
        assert!((options.clip.polya_min_fraction - 0.9).abs() < f64::EPSILON);
        for (option, value) in [
            ("--polya-min-frac", "0"),
            ("--polya-min-frac", "1.5"),
            ("--polya-min-frac", "x"),
            ("--min-clip-len", "0"),
            ("--polya-min-len", "0"),
            ("--min-mapq", "256"),
        ] {
            assert_eq!(
                scan(&[option, value]).unwrap_err().kind(),
                ErrorKind::ValueValidation,
                "{option} {value}"
            );
        }
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
