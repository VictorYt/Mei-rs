//! The `scan` subcommand: discovery of candidate insertion sites.

use anyhow::{Result, bail};
use tracing::info;

use crate::cli::ScanArgs;
use crate::inputs::Inputs;

/// Run `mei-rs scan`.
///
/// # Errors
///
/// Fails if the inputs are unusable or inconsistent. Signal extraction is not
/// implemented yet: once the inputs are checked, it always fails.
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
    bail!("signal extraction is not implemented yet")
}
