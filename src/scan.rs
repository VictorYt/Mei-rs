//! The `scan` subcommand: discovery of candidate insertion sites.

use anyhow::{Result, bail};
use tracing::info;

use crate::cli::ScanArgs;

/// Run `mei-rs scan`.
///
/// # Errors
///
/// Always fails for now: the pipeline is built in the next pull requests of
/// the v0.1.0 milestone.
pub fn run(args: &ScanArgs) -> Result<()> {
    info!(
        bam = %args.bam.display(),
        bed = %args.bed.display(),
        padding = args.padding,
        "scan"
    );
    bail!("the scan pipeline is not implemented yet")
}
