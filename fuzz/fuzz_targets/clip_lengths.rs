#![no_main]
//! Soft-clip lengths of a CIGAR made of the input bytes (operation, length).

use libfuzzer_sys::fuzz_target;
use mei_rs::reader::signals::clip_lengths;
use noodles::sam::alignment::record::cigar::op::Kind;

const KINDS: [Kind; 9] = [
    Kind::Match,
    Kind::Insertion,
    Kind::Deletion,
    Kind::Skip,
    Kind::SoftClip,
    Kind::HardClip,
    Kind::Pad,
    Kind::SequenceMatch,
    Kind::SequenceMismatch,
];

fuzz_target!(|data: &[u8]| {
    let ops: Vec<(Kind, usize)> = data
        .chunks_exact(3)
        .map(|c| (KINDS[usize::from(c[0]) % KINDS.len()], usize::from(u16::from_le_bytes([c[1], c[2]]))))
        .collect();
    let (left, right) = clip_lengths(&ops);
    let soft: usize = ops.iter().filter(|op| op.0 == Kind::SoftClip).map(|op| op.1).sum();
    assert!(left + right <= soft);
});
