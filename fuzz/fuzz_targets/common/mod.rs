//! Inputs built from the fuzz bytes, shared by the targets that run the
//! pipeline after clustering (`mates`, `writers`): signals around one site of
//! a small contig, whose clips and mates are random bases or pieces of the
//! default consensus sequences (so that typing is reached), then the
//! clustering, the typing of clips and mates, the filters, a genotype and
//! the quality score, as `scan` chains them.

#![allow(dead_code)]

use std::num::NonZeroU32;
use std::sync::LazyLock;

use mei_rs::cluster::filters::{self, FilterOptions};
use mei_rs::cluster::{Candidate, ClusterOptions, cluster};
use mei_rs::genotype::{self, AlleleCounts, score};
use mei_rs::reader::signals::{ClipOptions, Locus, Side, Signal, SignalKind, poly_tail};
use mei_rs::regions::bed::Probe;
use mei_rs::regions::targets::{Contig, TargetSet};
use mei_rs::typer::consensus::{ConsensusSet, reverse_complement};
use mei_rs::typer::{self, mates};

/// Length of the single contig.
pub const CONTIG_LENGTH: u64 = 10_000;

/// Reads the fuzz bytes in order; zeros once they are used up.
pub struct Input<'a>(pub &'a [u8]);

impl Input<'_> {
    pub fn u8(&mut self) -> u8 {
        let Some((&first, rest)) = self.0.split_first() else {
            return 0;
        };
        self.0 = rest;
        first
    }

    pub fn u16(&mut self) -> u16 {
        u16::from_le_bytes([self.u8(), self.u8()])
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Up to `max` bases: random (with N), or a piece of a consensus
    /// sequence, either strand, after an optional poly(A) or poly(T) run.
    pub fn sequence(&mut self, set: &ConsensusSet, max: usize) -> Vec<u8> {
        let mode = self.u8();
        let length = usize::from(self.u8()) % max.max(1) + 1;
        let mut sequence = Vec::new();
        let run = usize::from(mode >> 3);
        if mode & 2 != 0 {
            sequence.resize(run, if mode & 4 != 0 { b'T' } else { b'A' });
        }
        if mode & 1 != 0 && !set.sequences().is_empty() {
            let consensus =
                &set.sequences()[usize::from(self.u8()) % set.sequences().len()].sequence;
            let start = usize::from(self.u16()) % consensus.len();
            let piece = &consensus[start..(start + length).min(consensus.len())];
            if mode & 4 != 0 {
                sequence.extend(reverse_complement(piece));
            } else {
                sequence.extend_from_slice(piece);
            }
        } else {
            sequence.extend((0..length).map(|_| b"ACGTN"[usize::from(self.u8()) % 5]));
        }
        sequence.truncate(max);
        sequence
    }
}

/// The default consensus sequences, parsed once.
pub fn consensus() -> &'static ConsensusSet {
    static SET: LazyLock<ConsensusSet> =
        LazyLock::new(|| ConsensusSet::load(None).expect("the default consensus sequences"));
    &SET
}

/// Signals of up to 48 reads within 400 bp of position 5,000.
pub fn signals(input: &mut Input<'_>, set: &ConsensusSet) -> Vec<Signal> {
    let clip = ClipOptions::default();
    let mut signals = Vec::new();
    for read in 0..usize::from(input.u8() % 48) {
        if input.is_empty() {
            break;
        }
        let kind_byte = input.u8();
        let position = 4_800 + u64::from(input.u16() % 400);
        let side = if kind_byte & 1 == 0 {
            Side::Left
        } else {
            Side::Right
        };
        let start = position.saturating_sub(u64::from(input.u8()));
        let end = (position + u64::from(input.u8())).max(start + 1);
        let locus = Locus {
            contig: 0,
            span: NonZeroU32::new(u32::from(input.u8())),
            start: u64::from(input.u16()) % CONTIG_LENGTH,
        };
        let kind = match (kind_byte >> 1) % 6 {
            0 | 1 => {
                let sequence = input.sequence(set, 120);
                let tail = poly_tail(&sequence, side == Side::Right, &clip);
                SignalKind::SoftClip {
                    short: clip.is_short(sequence.len()),
                    sequence,
                    median_quality: input.u8() % 42,
                    tail,
                    supplementary: kind_byte & 0x80 != 0,
                }
            }
            2 => SignalKind::MateUnmapped,
            3 => SignalKind::MateLowMapq {
                mapq: input.u8() % 60,
                mate: (kind_byte & 0x40 != 0).then_some(locus),
            },
            4 => SignalKind::MateElsewhere { mate: locus },
            _ => SignalKind::LargeInsert {
                length: u64::from(input.u16()),
                mate: (kind_byte & 0x40 != 0).then_some(locus),
            },
        };
        signals.push(Signal {
            contig: 0,
            position,
            side,
            read: format!("r{read}"),
            first_segment: kind_byte & 0x20 != 0,
            kind,
            reverse: side == Side::Left,
            start,
            end,
            mapq: Some(60),
        });
    }
    signals.sort();
    signals
}

/// Candidates of the signals, typed, with their mates checked (sequences
/// from the input), filtered, genotyped and scored.
pub fn candidates(input: &mut Input<'_>, set: &ConsensusSet) -> Vec<Candidate> {
    let signals = signals(input, set);
    let mut candidates = cluster(signals, &ClusterOptions::default());
    let contigs = [Contig {
        name: "chr1".to_owned(),
        length: CONTIG_LENGTH,
    }];
    let probe = Probe {
        contig: "chr1".to_owned(),
        start: 4_900,
        end: 5_100,
        name: None,
    };
    let targets = TargetSet::new(vec![probe], 100, &contigs).expect("one valid probe");
    for candidate in &mut candidates {
        candidate.depth = Some(u64::from(input.u8()));
        candidate.typing = typer::type_candidate(set, candidate);
        mates::sample_mates(candidate);
        for check in &mut candidate.mate_checks {
            let byte = input.u8();
            check.local = byte & 1 != 0;
            check.sequence = (byte & 2 != 0).then(|| input.sequence(set, 150));
        }
        if mates::checked(candidate) {
            mates::type_mates(set, candidate);
        }
        candidate.reference_copy = input.u8() & 1 != 0;
        mates::discount_tail_mates(candidate);
        filters::apply(candidate, "chr1", &targets, &FilterOptions::default());
        let junctions = [
            AlleleCounts {
                alternative: u32::from(input.u8()),
                reference: u32::from(input.u8()),
            },
            AlleleCounts {
                alternative: u32::from(input.u8()),
                reference: u32::from(input.u8()),
            },
        ];
        candidate.genotype = genotype::genotype(
            &junctions,
            candidate.passes_before_quality(),
            candidate.depth,
        );
        candidate.qual = Some(score::quality(candidate));
    }
    candidates
}
