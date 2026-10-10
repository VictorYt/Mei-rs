//! The `scan` subcommand: discovery of candidate insertion sites.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use tracing::{info, warn};

use crate::cli::ScanArgs;
use crate::cluster::{Candidate, Confidence, Filter, cluster, filters};
use crate::inputs::Inputs;
use crate::reader::depth::local_depths;
use crate::reader::extract::{Extraction, extract};
use crate::writer::{self, Report};

/// Run `mei-rs scan`.
///
/// # Errors
///
/// Fails if the inputs are unusable or inconsistent, or if the output cannot
/// be written.
pub fn run(args: &ScanArgs) -> Result<()> {
    let inputs = Inputs::load(args)?;
    let targets = &inputs.targets;
    info!(
        "sample {}: {} {}",
        inputs.sample,
        inputs.alignment.format(),
        inputs.alignment.path().display()
    );
    info!(
        "{} targets ({} bp), {} regions to scan ({} bp with {} bp of padding)",
        targets.probes().len(),
        targets.probe_bases(),
        targets.regions().len(),
        targets.scanned_bases(),
        targets.padding()
    );

    let extract_options = args.extract_options();
    let extraction = extract(&inputs.alignment, targets, &extract_options)?;
    log_extraction(&extraction);

    let contigs: Vec<String> = inputs
        .alignment
        .contigs()
        .into_iter()
        .map(|c| c.name)
        .collect();
    let mut candidates = cluster(extraction.signals, &args.cluster_options());
    let clusters = candidates.len();
    candidates.retain(|c| !c.is_low_support());

    let boundaries: Vec<(String, u64)> = candidates
        .iter()
        .map(|c| (contigs[c.contig].clone(), c.breakpoint))
        .collect();
    let depths = local_depths(&inputs.alignment, &boundaries, &extract_options.filter)?;
    let filter_options = args.filter_options();
    for (candidate, depth) in candidates.iter_mut().zip(depths) {
        candidate.depth = Some(depth);
        filters::apply(
            candidate,
            &contigs[candidate.contig],
            targets,
            &filter_options,
        );
    }
    let pass = candidates.iter().filter(|c| c.is_pass()).count();
    let level = |confidence| {
        candidates
            .iter()
            .filter(|c| c.confidence() == Some(confidence))
            .count()
    };
    info!(
        "{clusters} clusters, {} with enough support: {pass} PASS ({} high, {} medium, {} low confidence), {} filtered{}",
        candidates.len(),
        level(Confidence::High),
        level(Confidence::Medium),
        level(Confidence::Low),
        candidates.len() - pass,
        if args.keep_filtered { " (kept)" } else { "" }
    );
    for filter in [
        Filter::NoJunction,
        Filter::NoTail,
        Filter::ProbeEdge,
        Filter::LowRatio,
    ] {
        let n = candidates
            .iter()
            .filter(|c| c.filters.contains(&filter))
            .count();
        info!("  {}: {n}", filter.name());
    }
    if !args.keep_filtered {
        candidates.retain(Candidate::is_pass);
    }

    let command = command_line();
    let report = Report {
        sample: &inputs.sample,
        command: &command,
        contigs: &contigs,
        candidates: &candidates,
    };
    write_output(&report, args)
}

/// Write the report to the output file, or to standard output.
fn write_output(report: &Report<'_>, args: &ScanArgs) -> Result<()> {
    let open = |path: &Path| {
        File::create(path).with_context(|| format!("cannot create {}", path.display()))
    };
    let mut out: Box<dyn Write> = match &args.output {
        Some(path) => Box::new(BufWriter::new(open(path)?)),
        None => Box::new(BufWriter::new(io::stdout().lock())),
    };
    writer::write(report, args.format, &mut out)
        .and_then(|()| out.flush())
        .context("cannot write the candidates")?;
    if let Some(path) = &args.output {
        info!(
            "{} candidates written to {}",
            report.candidates.len(),
            path.display()
        );
    }
    Ok(())
}

/// The command line, with the program name only (not its path).
fn command_line() -> String {
    let mut args = std::env::args();
    let program = args
        .next()
        .map(|p| {
            Path::new(&p)
                .file_name()
                .map_or(p.clone(), |name| name.to_string_lossy().into_owned())
        })
        .unwrap_or_default();
    std::iter::once(program)
        .chain(args)
        .collect::<Vec<_>>()
        .join(" ")
}

/// Log the read counts and the signals found.
fn log_extraction(extraction: &Extraction) {
    let stats = &extraction.stats;
    info!(
        "{} reads: {} kept, {} duplicates, {} below the MAPQ threshold, {} secondary or supplementary, {} QC fail, {} unmapped",
        stats.total,
        stats.kept,
        stats.duplicate,
        stats.low_mapq,
        stats.secondary + stats.supplementary,
        stats.qc_fail,
        stats.unmapped
    );
    if let Some(size) = extraction.insert_size {
        info!(
            "insert size: median {} bp, MAD {} bp ({} proper pairs); large inserts above {} bp",
            size.median, size.mad, size.pairs, size.threshold
        );
    } else {
        warn!("no proper pair: large inserts are not detected");
    }
    info!(
        "{} signals: {} soft-clips ({} with a poly(A/T) tail, {} more set aside for their base quality, {} of which are masked read ends), {} unmapped mates, {} low-MAPQ mates, {} mates on another contig, {} large inserts",
        extraction.signals.len(),
        extraction.count("soft_clip"),
        extraction.tailed_clips(),
        extraction.low_quality_clips,
        extraction.masked_clips,
        extraction.count("mate_unmapped"),
        extraction.count("mate_low_mapq"),
        extraction.count("mate_elsewhere"),
        extraction.count("large_insert")
    );
}
