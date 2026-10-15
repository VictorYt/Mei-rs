#![no_main]
//! The BED parser on any bytes: an error or probes, never a panic.

use std::io::Cursor;
use std::path::Path;

use libfuzzer_sys::fuzz_target;
use mei_rs::regions::bed;

fuzz_target!(|data: &[u8]| {
    if let Ok(probes) = bed::parse(Cursor::new(data), Path::new("fuzz.bed")) {
        assert!(!probes.is_empty());
        assert!(probes.iter().all(|p| p.start < p.end));
    }
});
