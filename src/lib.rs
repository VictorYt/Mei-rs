//! Mei-rs: detection of mobile element insertions (Alu, LINE-1, SVA) in
//! targeted capture and exome sequencing data.
//!
//! The pipeline runs in phases: targeted signal extraction ([`reader`]),
//! clustering and capture filters ([`cluster`]), element typing
//! ([`typer`]), then candidate export ([`writer`]). Local assembly and
//! genotyping come in later releases.

pub mod assembly;
pub mod cli;
pub mod cluster;
pub mod inputs;
pub mod reader;
pub mod regions;
pub mod scan;
pub mod typer;
pub mod writer;
