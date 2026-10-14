//! The `scan` subcommand: discovery of candidate insertion sites.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use noodles::core::{Position, Region};
use noodles::{bgzf, csi, fasta, tabix, vcf};
use rayon::prelude::*;
use tracing::{info, warn};

use crate::assembly::{self, Assembly};
use crate::cli::{OutputFormat, ScanArgs};
use crate::cluster::{Candidate, Confidence, Filter, cluster, filters, tracts};
use crate::genotype::{self, AlleleCounts};
use crate::inputs::Inputs;
use crate::reader::alleles::Site;
use crate::reader::breakpoints::BreakpointRequest;
use crate::reader::extract::ExtractOptions;
use crate::reader::extract::{Extraction, extract};
use crate::reader::filter::ReadFilter;
use crate::reader::junctions::JunctionRequest;
use crate::reader::signals::{ClipOptions, Signal, SignalKind};
use crate::reader::{self, alignment::AlignmentFile, depth::local_depths, mates::MateRequest};
use crate::regions::targets::Contig;
use crate::regions::targets::ScanRegion;
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
    let mut extraction = extract(&inputs.alignment, targets, &extract_options)?;
    let short_clips = keep_short_clips(&mut extraction.signals, &consensus);
    log_extraction(&extraction, short_clips, &extract_options.clip);

    let header_contigs = inputs.alignment.contigs();
    let contigs: Vec<String> = header_contigs.iter().map(|c| c.name.clone()).collect();
    let mut candidates = cluster(extraction.signals, &args.cluster_options());
    let clusters = candidates.len();
    candidates.retain(|c| !c.is_low_support());

    set_depths(
        &mut candidates,
        &inputs.alignment,
        &contigs,
        extract_options.filter,
    )?;
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
    locate_mate_only(
        &mut candidates,
        &inputs.alignment,
        &contigs,
        extract_options.filter,
    )?;
    if let Some(reference) = &args.reference {
        mates::flag_reference_copies(&mut candidates, &header_contigs, reference, &consensus)
            .with_context(|| format!("cannot read {}", reference.display()))?;
    }
    assemble_candidates(
        &mut candidates,
        &inputs.alignment,
        &header_contigs,
        args.reference.as_deref(),
        &extract_options,
    )?;
    let filter_options = args.filter_options();
    for candidate in &mut candidates {
        filters::apply(
            candidate,
            &contigs[candidate.contig],
            targets,
            &filter_options,
        );
    }
    log_candidates(&candidates, clusters, args.keep_filtered);
    if !args.keep_filtered {
        candidates.retain(Candidate::is_pass);
    }
    genotype_candidates(
        &mut candidates,
        &inputs.alignment,
        &contigs,
        &extract_options,
    )?;

    report(args, &inputs.sample, &header_contigs, &candidates)
}

/// Write the candidates in the output format.
fn report(
    args: &ScanArgs,
    sample: &str,
    header_contigs: &[Contig],
    candidates: &[Candidate],
) -> Result<()> {
    let contigs: Vec<String> = header_contigs.iter().map(|c| c.name.clone()).collect();
    let command = command_line();
    let contig_lengths: Vec<u64> = header_contigs.iter().map(|c| c.length).collect();
    let reference_bases = match (&args.reference, args.format) {
        (Some(reference), OutputFormat::Vcf) => {
            Some(reference_bases(candidates, header_contigs, reference)?)
        }
        _ => None,
    };
    let reference_name = args
        .reference
        .as_deref()
        .and_then(Path::file_name)
        .map(|n| n.to_string_lossy().into_owned());
    let report = Report {
        sample,
        command: &command,
        contigs: &contigs,
        contig_lengths: &contig_lengths,
        reference: reference_name.as_deref(),
        reference_bases: reference_bases.as_deref(),
        candidates,
    };
    write_output(&report, args)
}

/// Drop the short soft-clips (below `--min-clip-len`) without a poly(A/T)
/// tail or MEI sequence ([`typer::is_mei_clip`]), in parallel. Returns the
/// number of short clips and of those kept.
fn keep_short_clips(signals: &mut Vec<Signal>, consensus: &ConsensusSet) -> (usize, usize) {
    let short = |s: &Signal| matches!(s.kind, SignalKind::SoftClip { short: true, .. });
    let keep: Vec<bool> = signals
        .par_iter()
        .map(|s| !short(s) || typer::is_mei_clip(consensus, &s.kind))
        .collect();
    let total = signals.iter().filter(|s| short(s)).count();
    let mut keep = keep.into_iter();
    signals.retain(|_| keep.next().unwrap_or(true));
    let kept = signals.iter().filter(|s| short(s)).count();
    (total, kept)
}

/// Log the PASS candidates by confidence and the filtered ones by filter.
fn log_candidates(candidates: &[Candidate], clusters: usize, keep_filtered: bool) {
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
        if keep_filtered { " (kept)" } else { "" }
    );
    for filter in [
        Filter::NoJunction,
        Filter::NoTail,
        Filter::Untyped,
        Filter::NonMeiMates,
        Filter::TandemRepeat,
        Filter::ProbeEdge,
        Filter::LowRatio,
        Filter::LowQual,
    ] {
        let n = candidates
            .iter()
            .filter(|c| c.filters.contains(&filter))
            .count();
        info!("  {}: {n}", filter.name());
    }
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

/// Local depth of each candidate at its breakpoint.
fn set_depths(
    candidates: &mut [Candidate],
    alignment: &AlignmentFile,
    contigs: &[String],
    filter: ReadFilter,
) -> Result<()> {
    let boundaries: Vec<(String, u64)> = candidates
        .iter()
        .map(|c| (contigs[c.contig].clone(), c.breakpoint))
        .collect();
    for (candidate, depth) in
        candidates
            .iter_mut()
            .zip(local_depths(alignment, &boundaries, &filter)?)
    {
        candidate.depth = Some(depth);
    }
    Ok(())
}

/// Look for the breakpoint of the candidates placed by their mates among
/// the soft-clips of their reads ([`reader::breakpoints`]), for those with
/// MEI mates (see [`mates::MIN_TYPED_MATES`]): the others cannot pass. The
/// depth of the candidates moved is measured again at their breakpoint.
fn locate_mate_only(
    candidates: &mut [Candidate],
    alignment: &AlignmentFile,
    contigs: &[String],
    filter: ReadFilter,
) -> Result<()> {
    let mut selected: Vec<(usize, BreakpointRequest)> = candidates
        .iter()
        .enumerate()
        .filter(|(_, c)| c.mei_mates >= mates::MIN_TYPED_MATES)
        .filter_map(|(i, c)| {
            let (first, last) = c.mate_interval?;
            Some((
                i,
                BreakpointRequest {
                    contig: contigs[c.contig].clone(),
                    first,
                    last,
                    estimate: c.breakpoint,
                },
            ))
        })
        .collect();
    // Sorted by contig and first boundary, as the batches need.
    selected.sort_by_key(|(i, request)| (candidates[*i].contig, request.first));
    let requests: Vec<BreakpointRequest> = selected.iter().map(|(_, r)| r.clone()).collect();
    let found = reader::breakpoints::find_breakpoints(alignment, &requests, &filter)?;
    let mut moved: Vec<usize> = Vec::new();
    for ((i, _), breakpoint) in selected.iter().zip(found) {
        if let Some(breakpoint) = breakpoint {
            let c = &mut candidates[*i];
            c.breakpoint = breakpoint.position;
            c.breakpoint_clips = breakpoint.clips;
            c.mate_interval = None;
            moved.push(*i);
        }
    }
    let boundaries: Vec<(String, u64)> = moved
        .iter()
        .map(|&i| {
            (
                contigs[candidates[i].contig].clone(),
                candidates[i].breakpoint,
            )
        })
        .collect();
    for (&i, depth) in moved
        .iter()
        .zip(local_depths(alignment, &boundaries, &filter)?)
    {
        candidates[i].depth = Some(depth);
    }
    info!(
        "candidates placed by their mates, with MEI mates: {} looked for among the soft-clips of their reads, {} found",
        selected.len(),
        moved.len()
    );
    Ok(())
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
        "mates checked against the consensus sequences{}: {} of {} MEI mates; {} candidates typed by their mates",
        if reference.is_some() {
            ""
        } else {
            " (unmapped ones only: no --reference)"
        },
        candidates.iter().map(|c| c.mei_mates).sum::<usize>(),
        candidates.iter().map(|c| c.checked_mates).sum::<usize>(),
        candidates
            .iter()
            .filter(|c| c.typing.as_ref().is_some_and(|t| t.by_mates))
            .count()
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
    for candidate in candidates
        .iter_mut()
        .filter(|c| c.left_junction.is_some() || c.right_junction.is_some())
    {
        let contig = &contigs[candidate.contig];
        let (start, end) = tracts::window(candidate, contig.length);
        let Some(sequence) = reference_window(&mut reader, reference, &contig.name, start, end)?
        else {
            continue;
        };
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

/// Genotype the candidates with junctions ([`genotype`]) from the reads of
/// each allele at their junctions: alternative reads are clipped there by
/// at least `--min-short-clip-len` bases of median quality at least
/// `--min-clip-quality`, reference reads span the whole site by as many
/// bases on each side.
fn genotype_candidates(
    candidates: &mut [Candidate],
    alignment: &AlignmentFile,
    contigs: &[String],
    options: &ExtractOptions,
) -> Result<()> {
    let sites: Vec<Site> = candidates
        .iter()
        .map(|c| Site {
            contig: contigs[c.contig].clone(),
            left: c.left_junction,
            right: c.right_junction,
        })
        .collect();
    let counts = reader::alleles::allele_counts(
        alignment,
        &sites,
        options.filter,
        u64::from(options.clip.min_extracted_length()),
        options.clip.min_quality,
    )?;
    for ((c, site), counts) in candidates.iter_mut().zip(&sites).zip(counts) {
        // One count per junction, against the reference reads of the site.
        let junctions: Vec<AlleleCounts> = [(site.left, counts.left), (site.right, counts.right)]
            .into_iter()
            .filter(|(junction, _)| junction.is_some())
            .map(|(_, alternative)| AlleleCounts {
                alternative,
                reference: counts.reference,
            })
            .collect();
        c.genotype = genotype::genotype(&junctions);
    }
    Ok(())
}

/// Bases added on each side of a candidate to the window read again for its
/// assembly: unmapped mates are placed at their anchor, up to a read length
/// before the signals.
const ASSEMBLY_MARGIN: u64 = 300;

/// Assemble the junctions of the candidates that can pass (typed by their
/// clips or mates, with junction soft-clips: [`mates::may_pass`]), in parallel
/// ([`assembly`]); it needs the reference.
fn assemble_candidates(
    candidates: &mut [Candidate],
    alignment: &AlignmentFile,
    contigs: &[Contig],
    reference: Option<&Path>,
    options: &ExtractOptions,
) -> Result<()> {
    let Some(reference) = reference else {
        warn!("no --reference: the junctions are not assembled (TSD estimated from the clips)");
        return Ok(());
    };
    let mut reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(reference)
        .with_context(|| format!("cannot open the reference {}", reference.display()))?;
    // Candidates are sorted by contig and start, as the batches need.
    let mut selected: Vec<(usize, u64, Vec<u8>)> = Vec::new();
    let mut requests: Vec<JunctionRequest> = Vec::new();
    for (i, c) in candidates
        .iter()
        .enumerate()
        .filter(|(_, c)| mates::may_pass(c))
    {
        let contig = &contigs[c.contig];
        let (start, end) = tracts::window(c, contig.length);
        let Some(sequence) = reference_window(&mut reader, reference, &contig.name, start, end)?
        else {
            continue;
        };
        let junctions = [c.left_junction, c.right_junction];
        let first = junctions.iter().flatten().fold(c.start, |a, &j| a.min(j));
        let last = junctions.iter().flatten().fold(c.end, |a, &j| a.max(j));
        requests.push(JunctionRequest {
            window: ScanRegion {
                contig: contig.name.clone(),
                start: first.saturating_sub(ASSEMBLY_MARGIN),
                end: (last + ASSEMBLY_MARGIN).min(contig.length),
            },
            left: c.left_junction,
            right: c.right_junction,
        });
        selected.push((i, start, sequence));
    }
    let reads = reader::junctions::fetch_junction_reads(alignment, &requests, &options.filter)?;
    let assemblies: Vec<Assembly> = selected
        .par_iter()
        .zip(reads.par_iter())
        .map(|((i, start, sequence), reads)| {
            let c = &candidates[*i];
            assembly::assemble(
                reads,
                &tracts::ReferenceWindow {
                    start: *start,
                    sequence,
                },
                c.left_junction,
                c.right_junction,
                &options.clip,
            )
        })
        .collect();
    let mut counts = [0usize; 3];
    for ((i, _, _), assembly) in selected.iter().zip(assemblies) {
        counts[usize::from(assembly.left.is_some()) + usize::from(assembly.right.is_some())] += 1;
        candidates[*i].set_assembly(assembly);
    }
    info!(
        "junctions assembled for {} candidates: both {}, one {}, none {}",
        selected.len(),
        counts[2],
        counts[1],
        counts[0]
    );
    Ok(())
}

/// Reference bases of a 0-based half-open interval, upper case; `None` if
/// the interval is empty.
fn reference_window(
    reader: &mut fasta::io::IndexedReader<fasta::io::BufReader<std::fs::File>>,
    path: &Path,
    contig: &str,
    start: u64,
    end: u64,
) -> Result<Option<Vec<u8>>> {
    // Read as a 1-based inclusive region.
    let (Some(first), Some(last)) = (
        usize::try_from(start + 1)
            .ok()
            .and_then(|p| Position::try_from(p).ok()),
        usize::try_from(end)
            .ok()
            .and_then(|p| Position::try_from(p).ok()),
    ) else {
        return Ok(None);
    };
    if first > last {
        return Ok(None);
    }
    let region = Region::new(contig, first..=last);
    let record = reader
        .query(&region)
        .with_context(|| format!("cannot read {region} from {}", path.display()))?;
    let mut sequence = record.sequence().as_ref().to_vec();
    sequence.make_ascii_uppercase();
    Ok(Some(sequence))
}

/// Write the report to the output file, or to standard output.
fn write_output(report: &Report<'_>, args: &ScanArgs) -> Result<()> {
    let open = |path: &Path| {
        File::create(path).with_context(|| format!("cannot create {}", path.display()))
    };
    let compressed = args.format == OutputFormat::Vcf
        && args
            .output
            .as_deref()
            .is_some_and(|p| p.extension().is_some_and(|e| e == "gz"));
    let mut out: Box<dyn Write> = match &args.output {
        Some(path) if compressed => Box::new(bgzf::io::Writer::new(open(path)?)),
        Some(path) => Box::new(BufWriter::new(open(path)?)),
        None => Box::new(BufWriter::new(io::stdout().lock())),
    };
    writer::write(report, args.format, &mut out)
        .and_then(|()| out.flush())
        .context("cannot write the candidates")?;
    // The BGZF writer writes its end-of-file block when dropped.
    drop(out);
    if let Some(path) = &args.output {
        info!(
            "{} candidates written to {}",
            report.candidates.len(),
            path.display()
        );
        if compressed {
            index_vcf(path)?;
        }
    }
    Ok(())
}

/// Index a bgzipped VCF: tabix (`.tbi`), or CSI (`.csi`) for contigs longer
/// than 2^29 bp.
fn index_vcf(path: &Path) -> Result<()> {
    let index = vcf::fs::index(path).with_context(|| format!("cannot index {}", path.display()))?;
    let with = |extension: &str| {
        let mut name = path.as_os_str().to_owned();
        name.push(extension);
        std::path::PathBuf::from(name)
    };
    let written = match index {
        vcf::Index::Tabix(index) => {
            let dst = with(".tbi");
            tabix::fs::write(&dst, &index).map(|()| dst)
        }
        vcf::Index::Csi(index) => {
            let dst = with(".csi");
            csi::fs::write(&dst, &index).map(|()| dst)
        }
    }
    .with_context(|| format!("cannot write the index of {}", path.display()))?;
    info!("index written to {}", written.display());
    Ok(())
}

/// Reference base at the position of each candidate (the base right before
/// the insertion, VCF `REF`).
fn reference_bases(
    candidates: &[Candidate],
    contigs: &[Contig],
    reference: &Path,
) -> Result<Vec<u8>> {
    let mut reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(reference)
        .with_context(|| format!("cannot open the reference {}", reference.display()))?;
    candidates
        .iter()
        .map(|c| {
            let position = c.breakpoint.max(1);
            let base = reference_window(
                &mut reader,
                reference,
                &contigs[c.contig].name,
                position - 1,
                position,
            )?;
            Ok(base.and_then(|b| b.first().copied()).unwrap_or(b'N'))
        })
        .collect()
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

/// Log the read counts, the signals found and the short soft-clips kept.
fn log_extraction(extraction: &Extraction, (short, kept): (usize, usize), clip: &ClipOptions) {
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
    if clip.min_short_length < clip.min_length {
        info!(
            "{short} short soft-clips ({} to {} bp): {kept} kept with a poly(A/T) tail or MEI sequence",
            clip.min_short_length,
            clip.min_length - 1
        );
    }
}
