//! Mei-rs: detection of mobile element insertions (Alu, LINE-1, SVA) in
//! targeted capture and exome sequencing data.
//!
//! The pipeline runs in phases: targeted signal extraction ([`reader`]),
//! clustering and capture filters ([`cluster`]), then candidate export
//! ([`writer`]). Element typing, local assembly and
//! genotyping come in later releases.

pub mod cli;
pub mod cluster;
pub mod inputs;
pub mod reader;
pub mod regions;
pub mod scan;
pub mod writer;
