#![no_main]
//! Local assembly of the junctions on a reference window and reads built
//! from the input: reads are random bases or pieces of the window with
//! inserted bases, so that the walks leave the flanks.

mod common;

use libfuzzer_sys::fuzz_target;
use mei_rs::assembly::{MAX_TSD, MIN_INSERT, assemble};
use mei_rs::cluster::tracts::ReferenceWindow;
use mei_rs::reader::junctions::JunctionReads;
use mei_rs::reader::signals::ClipOptions;

fuzz_target!(|data: &[u8]| {
    let set = common::consensus();
    let mut input = common::Input(data);
    let start = u64::from(input.u16());
    let length = 60 + usize::from(input.u8()) * 2;
    let window: Vec<u8> = (0..length)
        .map(|_| b"ACGT"[usize::from(input.u8()) % 4])
        .collect();
    let left = (input.u8() & 1 != 0).then(|| start + u64::from(input.u16()) % (length as u64 + 1));
    let right = (input.u8() & 1 != 0).then(|| start + u64::from(input.u16()) % (length as u64 + 1));
    let mut reads = JunctionReads::default();
    for _ in 0..usize::from(input.u8() % 40) {
        if input.is_empty() {
            break;
        }
        let byte = input.u8();
        let read = if byte & 1 != 0 {
            // A piece of the window followed (or preceded) by inserted bases.
            let at = usize::from(input.u16()) % window.len();
            let flank = &window[at.saturating_sub(60)..at];
            let insert = input.sequence(set, 80);
            if byte & 2 != 0 {
                [flank, &insert].concat()
            } else {
                let after = &window[at..(at + 60).min(window.len())];
                [&insert[..], after].concat()
            }
        } else {
            input.sequence(set, 120)
        };
        match (byte >> 2) % 3 {
            0 => reads.left.push(read),
            1 => reads.right.push(read),
            _ => reads.inside.push(read),
        }
    }
    let reference = ReferenceWindow {
        start,
        sequence: &window,
    };
    let end = start + window.len() as u64;
    let assembly = assemble(&reads, &reference, left, right, &ClipOptions::default());
    for contig in [&assembly.left, &assembly.right].into_iter().flatten() {
        assert!((start..=end).contains(&contig.junction));
        assert!(contig.flank <= contig.contig.len());
        assert!(contig.contig.len() - contig.flank >= MIN_INSERT);
    }
    if let Some(tsd) = &assembly.tsd {
        assert!(assembly.both());
        assert!(!tsd.is_empty() && tsd.len() as u64 <= MAX_TSD);
    }
    if let Some(tail) = &assembly.tail {
        let longest = [&assembly.left, &assembly.right]
            .into_iter()
            .flatten()
            .map(|c| c.contig.len() - c.flank)
            .max()
            .unwrap_or(0);
        assert!(tail.length as usize <= longest);
    }
});
