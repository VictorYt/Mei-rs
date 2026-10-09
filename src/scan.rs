//! The `scan` subcommand: discovery of candidate insertion sites.

use anyhow::{Result, bail};
use tracing::{info, warn};

use crate::cli::ScanArgs;
use crate::inputs::Inputs;
use crate::reader::extract::{Extraction, extract};

/// Run `mei-rs scan`.
///
/// # Errors
///
/// Fails if the inputs are unusable or inconsistent. Clustering is not
/// implemented yet: once the signals are extracted, it always fails.
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

    let extraction = extract(&inputs.alignment, targets, &args.extract_options())?;
    log_extraction(&extraction);
    bail!("clustering is not implemented yet")
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
        "{} signals: {} soft-clips ({} with a poly(A/T) tail, {} more set aside for their base quality), {} unmapped mates, {} low-MAPQ mates, {} mates on another contig, {} large inserts",
        extraction.signals.len(),
        extraction.count("soft_clip"),
        extraction.tailed_clips(),
        extraction.low_quality_clips,
        extraction.count("mate_unmapped"),
        extraction.count("mate_low_mapq"),
        extraction.count("mate_elsewhere"),
        extraction.count("large_insert")
    );
}
