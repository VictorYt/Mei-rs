#![no_main]
//! The consensus FASTA parser on any text.

use libfuzzer_sys::fuzz_target;
use mei_rs::typer::consensus::ConsensusSet;

fuzz_target!(|data: &[u8]| {
    let set = ConsensusSet::parse(&String::from_utf8_lossy(data));
    assert!(set.sequences().iter().all(|c| !c.sequence.is_empty()));
});
