//! The `scan` subcommand: discovery of candidate insertion sites.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use noodles::core::{Position, Region};
use noodles::fasta;
use rayon::prelude::*;
use tracing::{info, warn};

use crate::cli::ScanArgs;
use crate::cluster::{Candidate, Confidence, Filter, cluster, filters, tracts};
use crate::inputs::Inputs;
use crate::reader::extract::{Extraction, extract};
use crate::reader::{self, alignment::AlignmentFile, depth::local_depths, mates::MateRequest};
use crate::regions::targets::Contig;
use crate::typer::{
    self,
    consensus::{ConsensusSet, Family},
    mates,
};
use crate::writer::{self, Report};

/// Run `mei-rs scan`.
///
/// # Errors
///
/// Fails if the inputs are unusable or inconsistent, or if the output cannot
/// be written.
pub fn run(args: &ScanArgs) -> Result<()> {
    let inputs = Inputs::load(args)?;
    let consensus = ConsensusSet::load(args.consensus.as_deref())?;
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

    let header_contigs = inputs.alignment.contigs();
    let contigs: Vec<String> = header_contigs.iter().map(|c| c.name.clone()).collect();
    let mut candidates = cluster(extraction.signals, &args.cluster_options());
    let clusters = candidates.len();
    candidates.retain(|c| !c.is_low_support());

    let boundaries: Vec<(String, u64)> = candidates
        .iter()
        .map(|c| (contigs[c.contig].clone(), c.breakpoint))
        .collect();
    let depths = local_depths(&inputs.alignment, &boundaries, &extract_options.filter)?;
    if let Some(reference) = &args.reference {
        flag_tandem_repeats(&mut candidates, &header_contigs, reference)?;
    } else {
        warn!(
            "no --reference: the tandem_repeat filter is off (length variants of reference poly(A) tracts and tandem repeats can pass as insertions)"
        );
    }
    type_candidates(&mut candidates, &consensus, args.consensus.as_deref());
    check_mates(
        &mut candidates,
        &inputs.alignment,
        &header_contigs,
        &consensus,
        args.reference.as_deref(),
    )?;
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
        Filter::Untyped,
        Filter::NonMeiMates,
        Filter::TandemRepeat,
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

/// Type the candidates with junction soft-clips ([`typer`]), in parallel.
fn type_candidates(candidates: &mut [Candidate], consensus: &ConsensusSet, path: Option<&Path>) {
    candidates
        .par_iter_mut()
        .filter(|c| c.left_junction.is_some() || c.right_junction.is_some())
        .for_each(|c| c.typing = typer::type_candidate(consensus, c));
    let typed = |family| {
        candidates
            .iter()
            .filter(|c| c.typing.as_ref().is_some_and(|t| t.family == family))
            .count()
    };
    info!(
        "{} consensus sequences ({}); candidates typed: {} ALU, {} LINE1, {} SVA",
        consensus.sequences().len(),
        path.map_or_else(
            || "built-in Dfam set".to_owned(),
            |p| p.display().to_string()
        ),
        typed(Family::Alu),
        typed(Family::Line1),
        typed(Family::Sva)
    );
}

/// Check a sample of the mates of each candidate against the consensus
/// sequences ([`mates`]); mapped mates need the reference.
fn check_mates(
    candidates: &mut [Candidate],
    alignment: &AlignmentFile,
    contigs: &[Contig],
    consensus: &ConsensusSet,
    reference: Option<&Path>,
) -> Result<()> {
    candidates.iter_mut().for_each(mates::sample_mates);
    let names: Vec<String> = contigs.iter().map(|c| c.name.clone()).collect();
    let requests = mates::unmapped_requests(candidates);
    let groups: Vec<Vec<MateRequest>> = requests
        .iter()
        .map(|(_, group)| group.iter().map(|&(_, r)| r).collect())
        .collect();
    let sequences = reader::mates::fetch_unmapped_mates(alignment, &names, &groups)?;
    for ((c, group), found) in requests.into_iter().zip(sequences) {
        for ((k, _), sequence) in group.into_iter().zip(found) {
            candidates[c].mate_checks[k].sequence = sequence;
        }
    }
    if let Some(reference) = reference {
        mates::fetch_mapped_mates(candidates, contigs, reference)
            .with_context(|| format!("cannot read the mate loci from {}", reference.display()))?;
    }
    candidates
        .par_iter_mut()
        .filter(|c| mates::checked(c))
        .for_each(|c| mates::type_mates(consensus, c));
    info!(
        "mates checked against the consensus sequences{}: {} of {} MEI mates",
        if reference.is_some() {
            ""
        } else {
            " (unmapped ones only: no --reference)"
        },
        candidates.iter().map(|c| c.mei_mates).sum::<usize>(),
        candidates.iter().map(|c| c.checked_mates).sum::<usize>()
    );
    Ok(())
}

/// Flag the candidates whose soft-clips are length variants of reference
/// tracts ([`tracts`]), reading the reference around each candidate with
/// junction soft-clips.
fn flag_tandem_repeats(
    candidates: &mut [Candidate],
    contigs: &[Contig],
    reference: &Path,
) -> Result<()> {
    let mut reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(reference)
        .with_context(|| format!("cannot open the reference {}", reference.display()))?;
    let mut sequence = Vec::new();
    for candidate in candidates
        .iter_mut()
        .filter(|c| c.left_junction.is_some() || c.right_junction.is_some())
    {
        let contig = &contigs[candidate.contig];
        // 0-based half-open, read as a 1-based inclusive region.
        let (start, end) = tracts::window(candidate, contig.length);
        let (Some(first), Some(last)) = (
            usize::try_from(start + 1)
                .ok()
                .and_then(|p| Position::try_from(p).ok()),
            usize::try_from(end)
                .ok()
                .and_then(|p| Position::try_from(p).ok()),
        ) else {
            continue;
        };
        let region = Region::new(contig.name.as_str(), first..=last);
        let record = reader
            .query(&region)
            .with_context(|| format!("cannot read {region} from {}", reference.display()))?;
        sequence.clear();
        sequence.extend_from_slice(record.sequence().as_ref());
        sequence.make_ascii_uppercase();
        tracts::apply(
            candidate,
            &tracts::ReferenceWindow {
                start,
                sequence: &sequence,
            },
        );
    }
    Ok(())
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
