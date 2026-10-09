//! Phase 1: targeted reading of BAM/CRAM files and extraction of the insertion
//! signals (soft-clipped reads, discordant pairs, unmapped mates).

pub mod alignment;
pub mod extract;
pub mod filter;
pub mod signals;
