//! False typing rate of short clips on real non-MEI clips (#133).
//!
//! ```bash
//! mei-rs scan --bam sample.cram --bed targets.bed --reference genome.fa \
//!     --profile discovery --keep-filtered --format json -o sample.json
//! cargo run --release --example clip_typing -- sample.json genome.fa
//! ```
//!
//! Short clips without a poly(A/T) tail are dropped before the JSON output
//! unless they are typed, so the clips measured are long ones, cut down to
//! the bases next to their junction (the part a shorter clip would hold).
//! Two kinds of clips (at least [`MIN_LENGTH`] bases, without a tail) are
//! known not to be insertions:
//! - local: found on either strand, with at most [`MAX_MISMATCHES`] of its
//!   bases different, in the reference within [`WINDOW`] bp of its junction
//!   (a local rearrangement or an alignment artefact), the matched reference
//!   bases not themselves typed to an MEI (a piece of a reference copy);
//! - adapter: holding the start of the Illumina adapter (read-through of a
//!   short fragment).
//!
//! Each clip is counted once (contig, junction, sequence); the pieces of
//! each length are typed with the default consensus set. Random sequences of
//! the same lengths are typed for comparison.

use std::collections::HashSet;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use noodles::core::{Position, Region};
use noodles::fasta;

use mei_rs::typer::consensus::{ConsensusSet, reverse_complement};
use mei_rs::typer::type_sequence;

/// Shortest whole clip used (bases).
const MIN_LENGTH: usize = 40;
/// Distance to the junction (bp) within which the clip is looked for.
const WINDOW: u64 = 1_000;
/// Largest fraction of the clip bases that may differ from the reference.
const MAX_MISMATCHES: f64 = 0.1;
/// Start of the Illumina adapter, and its reverse complement.
const ADAPTER: [&[u8]; 2] = [b"AGATCGGAAGAGC", b"GCTCTTCCGATCT"];
/// Lengths of the pieces typed.
const LENGTHS: [usize; 8] = [15, 17, 19, 20, 25, 30, 35, 40];
/// Random sequences typed per length.
const RANDOM: usize = 20_000;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let [_, json, reference] = &args[..] else {
        return Err("usage: clip_typing SAMPLE.json GENOME.fa".into());
    };
    let set = ConsensusSet::load(None)?;
    let report: serde_json::Value = serde_json::from_reader(BufReader::new(File::open(json)?))?;
    let mut fasta =
        fasta::io::indexed_reader::Builder::default().build_from_path(Path::new(reference))?;

    let mut seen = HashSet::new();
    let (mut clips, mut repeat) = (0, 0);
    // Clips and typed pieces per length: local, adapter.
    let mut measured = [0usize; 2];
    let mut typed = [[0usize; LENGTHS.len()]; 2];
    for candidate in report["candidates"].as_array().ok_or("no candidates")? {
        let contig = candidate["contig"].as_str().ok_or("no contig")?;
        for signal in candidate["signals"].as_array().ok_or("no signals")? {
            let (Some(sequence), None) = (signal["clipped_sequence"].as_str(), signal.get("tail"))
            else {
                continue;
            };
            let junction = signal["position"].as_u64().ok_or("no position")?;
            let sequence = sequence.as_bytes().to_vec();
            if sequence.len() < MIN_LENGTH
                || !seen.insert((contig.to_owned(), junction, sequence.clone()))
            {
                continue;
            }
            clips += 1;
            let kind = if ADAPTER
                .iter()
                .any(|a| sequence.windows(a.len()).any(|w| w == *a))
            {
                1
            } else {
                let start = junction.saturating_sub(WINDOW);
                // Near a contig end the window cannot be read: the clip is left out.
                let Ok(window) = read(&mut fasta, contig, start, junction + WINDOW) else {
                    continue;
                };
                let Some(at) = find(&window, &sequence) else {
                    continue;
                };
                let matched = &window[at..at + sequence.len()];
                if !type_sequence(&set, matched).is_empty() {
                    repeat += 1;
                    continue;
                }
                0
            };
            measured[kind] += 1;
            // The inserted bases next to the junction: the end of a left clip
            // (the insertion before the read), the start of a right one.
            let left = signal["side"].as_str() == Some("left");
            for (k, &n) in LENGTHS.iter().enumerate() {
                let piece = if left {
                    &sequence[sequence.len() - n..]
                } else {
                    &sequence[..n]
                };
                typed[kind][k] += usize::from(!type_sequence(&set, piece).is_empty());
            }
        }
    }
    println!(
        "clips without tail of {MIN_LENGTH}+ bases: {clips}; local: {}; adapter: {}; \
         local next to a reference MEI copy (left out): {repeat}",
        measured[0], measured[1]
    );
    println!(
        "length\tlocal_typed\tlocal\tadapter_typed\tadapter\trandom_typed\trandom\trandom_rate"
    );
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    for (k, &n) in LENGTHS.iter().enumerate() {
        let mut random_typed = 0;
        for _ in 0..RANDOM {
            let piece: Vec<u8> = (0..n)
                .map(|_| {
                    // xorshift64: a fixed sequence, no dependency.
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    b"ACGT"[(state >> 62) as usize]
                })
                .collect();
            random_typed += usize::from(!type_sequence(&set, &piece).is_empty());
        }
        #[allow(clippy::cast_precision_loss)]
        let rate = random_typed as f64 / RANDOM as f64;
        println!(
            "{n}\t{}\t{}\t{}\t{}\t{random_typed}\t{RANDOM}\t{rate:.4}",
            typed[0][k], measured[0], typed[1][k], measured[1]
        );
    }
    Ok(())
}

/// Reference bases of `[start, end)` (0-based), upper case.
fn read(
    fasta: &mut fasta::io::IndexedReader<fasta::io::BufReader<File>>,
    contig: &str,
    start: u64,
    end: u64,
) -> Result<Vec<u8>, Box<dyn std::error::Error>> {
    let first = Position::try_from(usize::try_from(start + 1)?)?;
    let last = Position::try_from(usize::try_from(end)?)?;
    let record = fasta.query(&Region::new(contig, first..=last))?;
    Ok(record.sequence().as_ref().to_ascii_uppercase())
}

/// Start of the best match of `needle` or of its reverse complement in
/// `haystack`, without gaps, with at most [`MAX_MISMATCHES`] of its bases
/// different.
fn find(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss
    )]
    let allowed = (needle.len() as f64 * MAX_MISMATCHES) as usize;
    let reverse = reverse_complement(needle);
    let mismatches = |w: &[u8], n: &[u8]| w.iter().zip(n).filter(|(a, b)| a != b).count();
    haystack
        .windows(needle.len())
        .enumerate()
        .map(|(at, w)| (mismatches(w, needle).min(mismatches(w, &reverse)), at))
        .filter(|&(m, _)| m <= allowed)
        .min()
        .map(|(_, at)| at)
}
