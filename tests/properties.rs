//! Property tests of the parsers and of the sequence routines: no panic on
//! any input, and the invariants of their results. The same functions are
//! fuzzed by `fuzz/` (`cargo fuzz`), with coverage-guided inputs.

use std::fmt::Write;
use std::io::Cursor;
use std::path::Path;

use noodles::sam::alignment::record::cigar::op::Kind;
use proptest::prelude::*;

use mei_rs::reader::signals::{ClipOptions, clip_lengths, poly_tail};
use mei_rs::regions::bed;
use mei_rs::typer::align;
use mei_rs::typer::consensus::{ConsensusSet, kmers, reverse_complement};
use mei_rs::typer::{self, MIN_IDENTITY};

fn dna(max: usize) -> impl Strategy<Value = Vec<u8>> {
    prop::collection::vec(prop::sample::select(b"ACGT".to_vec()), 0..max)
}

fn dna_with_n(max: usize) -> impl Strategy<Value = Vec<u8>> {
    prop::collection::vec(prop::sample::select(b"ACGTNacgtn".to_vec()), 0..max)
}

fn kind() -> impl Strategy<Value = Kind> {
    prop::sample::select(vec![
        Kind::Match,
        Kind::Insertion,
        Kind::Deletion,
        Kind::Skip,
        Kind::SoftClip,
        Kind::HardClip,
        Kind::Pad,
        Kind::SequenceMatch,
        Kind::SequenceMismatch,
    ])
}

proptest! {
    /// Any bytes: an error or well-formed probes, never a panic.
    #[test]
    fn bed_parse_any_bytes(bytes in prop::collection::vec(any::<u8>(), 0..400)) {
        if let Ok(probes) = bed::parse(Cursor::new(bytes), Path::new("fuzz.bed")) {
            prop_assert!(!probes.is_empty());
            for probe in probes {
                prop_assert!(probe.start < probe.end);
                prop_assert!(!probe.contig.is_empty());
            }
        }
    }

    /// Text close to BED (tabs, numbers, comments): valid lines come back
    /// unchanged and in order.
    #[test]
    fn bed_parse_valid_lines(
        lines in prop::collection::vec(("chr[0-9XY]{1,2}", 0u64..1_000_000, 1u64..10_000), 1..30),
    ) {
        let mut text = String::new();
        for (contig, start, length) in &lines {
            writeln!(text, "{contig}\t{start}\t{}", start + length).unwrap();
        }
        let probes = bed::parse(Cursor::new(format!("# comment\n{text}")), Path::new("a.bed")).unwrap();
        prop_assert_eq!(probes.len(), lines.len());
        for (probe, (contig, start, length)) in probes.iter().zip(&lines) {
            prop_assert_eq!(&probe.contig, contig);
            prop_assert_eq!(probe.start, *start);
            prop_assert_eq!(probe.end, start + length);
        }
    }

    /// Soft-clip lengths come from the first and last non-hard-clip operations.
    #[test]
    fn clip_lengths_follow_the_ends(ops in prop::collection::vec((kind(), 1usize..500), 0..10)) {
        let (left, right) = clip_lengths(&ops);
        let aligned: Vec<_> = ops.iter().filter(|op| op.0 != Kind::HardClip).collect();
        if aligned.len() <= 1 {
            prop_assert_eq!((left, right), (0, 0));
        } else {
            let expect = |op: &(Kind, usize)| if op.0 == Kind::SoftClip { op.1 } else { 0 };
            prop_assert_eq!(left, expect(aligned[0]));
            prop_assert_eq!(right, expect(aligned[aligned.len() - 1]));
        }
        let soft: usize = ops.iter().filter(|op| op.0 == Kind::SoftClip).map(|op| op.1).sum();
        prop_assert!(left + right <= soft);
    }

    /// A tail is made of A or T, at least `polya_min_length` long, within the
    /// clip, and holds the required fraction of its base; it ends on that base.
    #[test]
    #[allow(clippy::cast_precision_loss)]
    fn poly_tail_invariants(
        clip in dna_with_n(80),
        junction_at_end in any::<bool>(),
        min_length in 1u32..20,
        fraction in 0.5f64..1.0,
    ) {
        let options = ClipOptions {
            polya_min_length: min_length,
            polya_min_fraction: fraction,
            ..ClipOptions::default()
        };
        if let Some(tail) = poly_tail(&clip, junction_at_end, &options) {
            prop_assert!(tail.base == b'A' || tail.base == b'T');
            prop_assert!(tail.length >= min_length);
            prop_assert!(tail.length as usize <= clip.len());
            let from_junction: Vec<u8> = if junction_at_end {
                clip.iter().rev().copied().collect()
            } else {
                clip.clone()
            };
            let span = &from_junction[..tail.length as usize];
            let matching = span.iter().filter(|b| b.eq_ignore_ascii_case(&tail.base)).count();
            prop_assert!(matching as f64 >= fraction * f64::from(tail.length) - 1e-9);
            prop_assert!(span[span.len() - 1].eq_ignore_ascii_case(&tail.base));
        }
    }

    /// A pure run is found whole, on either side.
    #[test]
    fn pure_runs_are_found(length in 10usize..60, base in prop::sample::select(vec![b'A', b'T']), junction_at_end in any::<bool>()) {
        let run = vec![base; length];
        let tail = poly_tail(&run, junction_at_end, &ClipOptions::default()).unwrap();
        prop_assert_eq!(tail.base, base);
        prop_assert_eq!(tail.length as usize, length);
    }

    #[test]
    fn reverse_complement_is_an_involution(sequence in dna(200)) {
        prop_assert_eq!(reverse_complement(&reverse_complement(&sequence)), sequence);
    }

    /// One k-mer per window, at increasing positions.
    #[test]
    fn kmers_are_ordered(sequence in dna_with_n(120)) {
        let positions: Vec<usize> = kmers(&sequence).map(|(p, _)| p).collect();
        prop_assert!(positions.windows(2).all(|w| w[0] < w[1]));
        prop_assert!(positions.iter().all(|&p| p + 11 <= sequence.len()));
    }

    /// Any text: a set of consensus sequences without a panic.
    #[test]
    fn consensus_parse_any_text(text in "\\PC{0,300}") {
        let _ = ConsensusSet::parse(&text);
    }

    /// Alignment coordinates stay within the sequences; identity is a fraction.
    #[test]
    fn local_alignment_bounds(query in dna_with_n(60), target in dna_with_n(80)) {
        if let Some(a) = align::local(&query, &target) {
            prop_assert!(a.query_start < a.query_end && a.query_end <= query.len());
            prop_assert!(a.target_start < a.target_end && a.target_end <= target.len());
            prop_assert!(a.matches <= a.columns);
            prop_assert!(a.score > 0);
            prop_assert!((0.0..=1.0).contains(&a.identity()));
        }
    }

    /// A piece of the target aligns to it exactly.
    #[test]
    fn local_alignment_finds_a_substring(target in dna(120), a in 0usize..120, b in 0usize..120) {
        let (lo, hi) = (a.min(b).min(target.len()), a.max(b).min(target.len()));
        prop_assume!(hi - lo >= 8);
        let query = &target[lo..hi];
        let alignment = align::local(query, &target).unwrap();
        prop_assert!(alignment.matches >= query.len());
        prop_assert!(alignment.score >= align::MATCH * i32::try_from(query.len()).unwrap());
    }

    /// Random sequences are typed without a panic, and every hit respects the
    /// thresholds.
    #[test]
    fn typing_random_sequences(sequence in dna(200)) {
        let set = ConsensusSet::load(None).unwrap();
        let hits = typer::type_sequence(&set, &sequence);
        prop_assert!(hits.len() <= typer::ALIGNED_STRANDS);
        for hit in hits {
            prop_assert!(hit.identity >= MIN_IDENTITY && hit.identity <= 1.0);
            prop_assert!(hit.element_start < hit.element_end);
            prop_assert!(hit.element_end <= set.sequences()[hit.consensus].sequence.len());
        }
    }

    /// A piece of a consensus, on either strand, types to its family.
    #[test]
    fn typing_a_consensus_piece(index in 0usize..1000, start in 0usize..1000, length in 40usize..120, reverse in any::<bool>()) {
        let set = ConsensusSet::load(None).unwrap();
        let consensus = &set.sequences()[index % set.sequences().len()];
        let lo = start % consensus.sequence.len();
        let hi = (lo + length).min(consensus.sequence.len());
        prop_assume!(hi - lo >= 40);
        let mut piece = consensus.sequence[lo..hi].to_vec();
        // A piece that is mostly poly(A) (the tail of the consensus) is not
        // typed by design: a tail carries no information on the element.
        let most_common = b"ACGT".iter().map(|b| piece.iter().filter(|p| *p == b).count()).max().unwrap();
        prop_assume!(most_common * 10 < piece.len() * 6);
        if reverse {
            piece = reverse_complement(&piece);
        }
        let hits = typer::type_sequence(&set, &piece);
        prop_assert!(
            hits.iter().any(|h| set.sequences()[h.consensus].family == consensus.family),
            "{} {lo}..{hi} not typed",
            consensus.name
        );
    }
}
