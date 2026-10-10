//! Mei-rs entry point.

use anyhow::{Context, Result};
use clap::Parser;
use tracing_subscriber::EnvFilter;

use mei_rs::cli::{Cli, Command};
use mei_rs::scan;

fn main() -> Result<()> {
    let cli = Cli::parse();
    init_logging(cli.verbose, cli.quiet);

    if let Some(threads) = cli.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads.get())
            .build_global()
            .context("cannot start the thread pool")?;
    }

    match cli.command {
        Command::Scan(args) => scan::run(&args),
    }
}

/// Log to standard error, which keeps standard output for the results.
/// `RUST_LOG` takes precedence over `-v` and `-q`.
fn init_logging(verbose: u8, quiet: bool) {
    let level = match (quiet, verbose) {
        (true, _) => "error",
        (false, 0) => "info",
        (false, 1) => "debug",
        (false, _) => "trace",
    };
    let filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new(level));
    tracing_subscriber::fmt()
        .with_env_filter(filter)
        .with_writer(std::io::stderr)
        .with_target(false)
        .init();
}
