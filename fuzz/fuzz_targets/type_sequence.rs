#![no_main]
//! Typing of any sequence (any bytes) against the built-in consensus set.

use std::sync::OnceLock;

use libfuzzer_sys::fuzz_target;
use mei_rs::typer::consensus::ConsensusSet;
use mei_rs::typer::{self, MIN_IDENTITY};

fuzz_target!(|data: &[u8]| {
    static SET: OnceLock<ConsensusSet> = OnceLock::new();
    let set = SET.get_or_init(|| ConsensusSet::load(None).unwrap());
    // Keep the sequence at the size of a read (longer inputs only slow the run).
    let sequence = &data[..data.len().min(400)];
    for hit in typer::type_sequence(set, sequence) {
        assert!(hit.identity >= MIN_IDENTITY && hit.identity <= 1.0);
        assert!(hit.element_start < hit.element_end);
    }
});
