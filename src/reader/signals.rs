//! Insertion signals carried by one anchored read.
//!
//! - **Soft-clips** (split reads): the clipped part of the read is the start
//!   of the inserted sequence; the junction gives the breakpoint to the base.
//! - **Mate signals** (discordant pairs): the read is anchored near the
//!   insertion and its mate lies inside the inserted element, so the mate is
//!   unmapped, multi-mapped (low MAPQ) or mapped to another copy of the
//!   element elsewhere in the genome.
//!
//! Every signal says on which [`Side`] of the anchored read the inserted
//! sequence lies, and gives a *boundary* position (between two bases,
//! 0-based: boundary `b` lies between bases `b - 1` and `b`).

use std::io;

use noodles::sam::Header;
use noodles::sam::alignment::Record;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;

use super::filter::mapq;

/// Side of the anchored read where the inserted sequence lies.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Side {
    /// Before the read (left soft-clip, or reverse read whose mate is inside
    /// the insertion).
    Left,
    /// After the read (right soft-clip, or forward read whose mate is inside
    /// the insertion).
    Right,
}

/// A poly(A) or poly(T) tail at the junction end of a soft-clip.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PolyTail {
    /// `b'A'` or `b'T'`.
    pub base: u8,
    /// Length of the tail, from the junction.
    pub length: u32,
}

/// What a signal is made of.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SignalKind {
    /// Soft-clipped read.
    SoftClip {
        /// Clipped sequence, in reference orientation.
        sequence: Vec<u8>,
        /// Median base quality of the clipped part.
        median_quality: u8,
        /// Poly(A/T) tail at the junction end of the clip.
        tail: Option<PolyTail>,
        /// The read has supplementary alignments (`SA` tag).
        supplementary: bool,
    },
    /// The mate is unmapped.
    MateUnmapped,
    /// The mate maps with a low mapping quality (`MQ` tag).
    MateLowMapq {
        /// Mapping quality of the mate.
        mapq: u8,
    },
    /// The mate maps to another contig.
    MateElsewhere {
        /// Contig of the mate (index in the header).
        contig: usize,
    },
    /// The mate maps to the same contig, unusually far away.
    LargeInsert {
        /// Absolute template length.
        length: u64,
    },
}

impl SignalKind {
    /// Short name, for logs and outputs.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Self::SoftClip { .. } => "soft_clip",
            Self::MateUnmapped => "mate_unmapped",
            Self::MateLowMapq { .. } => "mate_low_mapq",
            Self::MateElsewhere { .. } => "mate_elsewhere",
            Self::LargeInsert { .. } => "large_insert",
        }
    }

    /// Whether this is a mate (discordant pair) signal.
    #[must_use]
    pub fn is_mate_signal(&self) -> bool {
        !matches!(self, Self::SoftClip { .. })
    }
}

/// One insertion signal, anchored by one read.
///
/// Signals sort by contig, position, side, then read.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Signal {
    /// Contig (index in the header).
    pub contig: usize,
    /// Boundary nearest to the insertion: the junction of a soft-clip, or the
    /// end of the read facing the insertion for a mate signal.
    pub position: u64,
    /// Side of the read where the inserted sequence lies.
    pub side: Side,
    /// Read name.
    pub read: String,
    /// First segment of its template (read 1).
    pub first_segment: bool,
    /// Kind of signal.
    pub kind: SignalKind,
    /// The read maps to the reverse strand.
    pub reverse: bool,
    /// Alignment start of the read, 0-based.
    pub start: u64,
    /// Alignment end of the read, 0-based exclusive.
    pub end: u64,
    /// Mapping quality of the read (`None` if unavailable).
    pub mapq: Option<u8>,
}

/// Soft-clip and poly(A/T) detection options.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ClipOptions {
    /// Minimum length of a soft-clip.
    pub min_length: u32,
    /// Minimum median base quality of a soft-clip.
    pub min_quality: u8,
    /// Minimum length of a poly(A/T) tail.
    pub polya_min_length: u32,
    /// Minimum fraction of A (or T) in a poly(A/T) tail.
    pub polya_min_fraction: f64,
}

impl Default for ClipOptions {
    fn default() -> Self {
        Self {
            min_length: 20,
            min_quality: 20,
            polya_min_length: 10,
            polya_min_fraction: 0.8,
        }
    }
}

/// Signals and side information found in one read.
#[derive(Debug, Default)]
pub struct ReadSignals {
    /// Signals of the read (soft-clips and mate signals).
    pub signals: Vec<Signal>,
    /// Same-contig pair that is not a proper pair: a large insert if its
    /// template length exceeds the threshold of the sample.
    pub large_insert: Option<Signal>,
    /// Template length of a proper pair (read 1 only), to estimate the insert
    /// size distribution.
    pub insert_size: Option<u64>,
    /// Soft-clips set aside for their low base quality.
    pub low_quality_clips: u32,
}

/// Find the signals of a read that passed the filter.
///
/// # Errors
///
/// Fails if a field of the record cannot be decoded.
pub fn read_signals(
    record: &dyn Record,
    header: &Header,
    clip: &ClipOptions,
    min_mapq: u8,
) -> io::Result<ReadSignals> {
    let mut out = ReadSignals::default();
    let (Some(contig), Some(start), Some(end)) = (
        record.reference_sequence_id(header).transpose()?,
        record.alignment_start().transpose()?,
        record.alignment_end().transpose()?,
    ) else {
        return Ok(out);
    };
    let flags = record.flags()?;
    let anchor = Signal {
        contig,
        position: 0,
        side: Side::Left,
        read: record.name().map_or_else(String::new, ToString::to_string),
        first_segment: flags.is_first_segment(),
        kind: SignalKind::MateUnmapped,
        reverse: flags.is_reverse_complemented(),
        start: start.get() as u64 - 1,
        end: end.get() as u64,
        mapq: mapq(record)?,
    };

    clip_signals(record, &anchor, clip, &mut out)?;

    if flags.is_segmented() {
        // The inserted sequence lies on the side the read points to.
        let (position, side) = if anchor.reverse {
            (anchor.start, Side::Left)
        } else {
            (anchor.end, Side::Right)
        };
        let mate = |kind| Signal {
            position,
            side,
            kind,
            ..anchor.clone()
        };
        let proper = flags.is_properly_segmented();
        if flags.is_mate_unmapped() {
            out.signals.push(mate(SignalKind::MateUnmapped));
        } else {
            let mate_contig = record.mate_reference_sequence_id(header).transpose()?;
            let mate_mapq = mate_mapq(record)?;
            if let Some(mapq) = mate_mapq.filter(|&q| q < min_mapq && !proper) {
                out.signals.push(mate(SignalKind::MateLowMapq { mapq }));
            } else if let Some(other) = mate_contig.filter(|&c| c != contig) {
                out.signals
                    .push(mate(SignalKind::MateElsewhere { contig: other }));
            } else {
                let length = u64::from(record.template_length()?.unsigned_abs());
                if !proper {
                    out.large_insert = Some(mate(SignalKind::LargeInsert { length }));
                } else if anchor.first_segment && length > 0 {
                    out.insert_size = Some(length);
                }
            }
        }
    }
    Ok(out)
}

/// Add the soft-clip signals of a read (one per clipped end).
fn clip_signals(
    record: &dyn Record,
    anchor: &Signal,
    options: &ClipOptions,
    out: &mut ReadSignals,
) -> io::Result<()> {
    let ops = record
        .cigar()
        .iter()
        .map(|op| op.map(|op| (op.kind(), op.len())))
        .collect::<io::Result<Vec<_>>>()?;
    let (left, right) = clip_lengths(&ops);
    let min = options.min_length as usize;
    if left < min && right < min {
        return Ok(());
    }

    let sequence: Vec<u8> = record.sequence().iter().collect();
    let qualities = record
        .quality_scores()
        .iter()
        .collect::<io::Result<Vec<u8>>>()?;
    let supplementary = record.data().get(&Tag::OTHER_ALIGNMENTS).is_some();
    let n = sequence.len();

    for (length, side, position) in [
        (left, Side::Left, anchor.start),
        (right, Side::Right, anchor.end),
    ] {
        if length < min || length > n {
            continue;
        }
        let range = match side {
            Side::Left => 0..length,
            Side::Right => n - length..n,
        };
        let median_quality = qualities.get(range.clone()).map_or(u8::MAX, median);
        if median_quality < options.min_quality {
            out.low_quality_clips += 1;
            continue;
        }
        let clipped = &sequence[range];
        // The tail is at the junction end of the clip.
        let tail = poly_tail(clipped, side == Side::Left, options);
        out.signals.push(Signal {
            position,
            side,
            kind: SignalKind::SoftClip {
                sequence: clipped.to_vec(),
                median_quality,
                tail,
                supplementary,
            },
            ..anchor.clone()
        });
    }
    Ok(())
}

/// Soft-clipped lengths at the start and end of a CIGAR (hard clips skipped).
fn clip_lengths(ops: &[(Kind, usize)]) -> (usize, usize) {
    let soft = |op: Option<&(Kind, usize)>| match op {
        Some(&(Kind::SoftClip, len)) => len,
        _ => 0,
    };
    let not_hard = |op: &&(Kind, usize)| op.0 != Kind::HardClip;
    let left = soft(ops.iter().find(not_hard));
    let right = soft(ops.iter().rev().find(not_hard));
    if ops.iter().filter(|op| op.0 != Kind::HardClip).count() == 1 {
        // Fully clipped: not an alignment.
        return (0, 0);
    }
    (left, right)
}

/// Longest poly(A) or poly(T) tail starting at the junction end of a clip.
///
/// `junction_at_end` is true when the junction is at the end of `clip` (left
/// soft-clip). The tail extends while the fraction of A (or T) stays at least
/// `polya_min_fraction` and ends on a matching base, so isolated sequencing
/// errors are tolerated.
fn poly_tail(clip: &[u8], junction_at_end: bool, options: &ClipOptions) -> Option<PolyTail> {
    let mut best: Option<PolyTail> = None;
    for base in *b"AT" {
        let bases: Box<dyn Iterator<Item = &u8>> = if junction_at_end {
            Box::new(clip.iter().rev())
        } else {
            Box::new(clip.iter())
        };
        let mut matches = 0u32;
        let mut length = 0u32;
        for (i, b) in (1u32..).zip(bases) {
            if b.eq_ignore_ascii_case(&base) {
                matches += 1;
                if f64::from(matches) >= options.polya_min_fraction * f64::from(i) {
                    length = i;
                }
            }
        }
        if length >= options.polya_min_length && best.is_none_or(|t| length > t.length) {
            best = Some(PolyTail { base, length });
        }
    }
    best
}

/// Median of the qualities (lower median), 255 if empty or unavailable.
fn median(qualities: &[u8]) -> u8 {
    if qualities.is_empty() {
        return u8::MAX;
    }
    let mut sorted = qualities.to_vec();
    sorted.sort_unstable();
    sorted[(sorted.len() - 1) / 2]
}

/// Mapping quality of the mate (`MQ` tag), if present.
fn mate_mapq(record: &dyn Record) -> io::Result<Option<u8>> {
    record
        .data()
        .get(&Tag::MATE_MAPPING_QUALITY)
        .transpose()
        .map(|value| {
            value
                .and_then(|v| v.as_int())
                .and_then(|q| u8::try_from(q).ok())
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::core::Position;
    use noodles::sam::alignment::RecordBuf;
    use noodles::sam::alignment::record::cigar::Op;
    use noodles::sam::alignment::record::{Flags, MappingQuality};
    use noodles::sam::alignment::record_buf::data::field::Value;
    use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};
    use noodles::sam::header::record::value::{Map, map::ReferenceSequence};
    use std::num::NonZero;

    fn header() -> Header {
        let contig = || Map::<ReferenceSequence>::new(NonZero::new(100_000).unwrap());
        Header::builder()
            .add_reference_sequence("chr1", contig())
            .add_reference_sequence("chr2", contig())
            .build()
    }

    fn cigar(spec: &[(Kind, usize)]) -> Cigar {
        spec.iter().map(|&(kind, len)| Op::new(kind, len)).collect()
    }

    /// A read at 1-based position 1001, with the given CIGAR and sequence.
    fn read(
        flags: Flags,
        ops: &[(Kind, usize)],
        sequence: &str,
    ) -> noodles::sam::alignment::record_buf::Builder {
        RecordBuf::builder()
            .set_name("r1")
            .set_flags(flags)
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(1001).unwrap())
            .set_mapping_quality(MappingQuality::new(60).unwrap())
            .set_cigar(cigar(ops))
            .set_sequence(Sequence::from(sequence.as_bytes().to_vec()))
            .set_quality_scores(QualityScores::from(vec![30; sequence.len()]))
    }

    fn signals(record: &RecordBuf) -> ReadSignals {
        read_signals(record, &header(), &ClipOptions::default(), 20).unwrap()
    }

    #[test]
    fn clip_lengths_skip_hard_clips() {
        use Kind::{HardClip, Match, SoftClip};
        assert_eq!(clip_lengths(&[(SoftClip, 30), (Match, 70)]), (30, 0));
        assert_eq!(clip_lengths(&[(Match, 70), (SoftClip, 30)]), (0, 30));
        assert_eq!(
            clip_lengths(&[
                (HardClip, 5),
                (SoftClip, 25),
                (Match, 50),
                (SoftClip, 20),
                (HardClip, 3)
            ]),
            (25, 20)
        );
        assert_eq!(clip_lengths(&[(Match, 100)]), (0, 0));
        assert_eq!(clip_lengths(&[(SoftClip, 100)]), (0, 0));
    }

    #[test]
    fn poly_tails() {
        let options = ClipOptions::default();
        let tail = |clip: &str, at_end| poly_tail(clip.as_bytes(), at_end, &options);
        assert_eq!(
            tail("AAAAAAAAAAAACGTGCC", false),
            Some(PolyTail {
                base: b'A',
                length: 12
            })
        );
        assert_eq!(
            tail("GCAGCGTTTTTTTTTTT", true),
            Some(PolyTail {
                base: b'T',
                length: 11
            })
        );
        // One sequencing error in the tail, and one at the junction.
        assert_eq!(
            tail("AAAAGAAAAAAAGTCCGT", false),
            Some(PolyTail {
                base: b'A',
                length: 12
            })
        );
        assert_eq!(
            tail("CAAAAAAAAAAA", false),
            Some(PolyTail {
                base: b'A',
                length: 12
            })
        );
        // Too short, or not at the junction.
        assert_eq!(tail("AAAAAAAAACGT", false), None);
        assert_eq!(tail("CGTCGTAAAAAAAAAAAA", false), None);
        // An A-rich sequence is not a tail.
        assert_eq!(tail("AACAGATACAAGTAACAGTA", false), None);
    }

    #[test]
    fn median_quality() {
        assert_eq!(median(&[30, 10, 20]), 20);
        assert_eq!(median(&[30, 10, 20, 40]), 20);
        assert_eq!(median(&[]), 255);
    }

    #[test]
    fn soft_clips_on_both_sides() {
        let seq = format!("{}{}{}", "A".repeat(20), "C".repeat(60), "GT".repeat(10));
        let record = read(
            Flags::empty(),
            &[
                (Kind::SoftClip, 20),
                (Kind::Match, 60),
                (Kind::SoftClip, 20),
            ],
            &seq,
        )
        .build();
        let out = signals(&record);
        let clips: Vec<_> = out
            .signals
            .iter()
            .map(|s| (s.position, s.side, s.start, s.end))
            .collect();
        assert_eq!(
            clips,
            [
                (1000, Side::Left, 1000, 1060),
                (1060, Side::Right, 1000, 1060)
            ]
        );
        let SignalKind::SoftClip {
            sequence,
            tail,
            supplementary,
            ..
        } = &out.signals[0].kind
        else {
            panic!("not a soft-clip");
        };
        assert_eq!(sequence, "A".repeat(20).as_bytes());
        assert_eq!(
            *tail,
            Some(PolyTail {
                base: b'A',
                length: 20
            })
        );
        assert!(!supplementary);
        let SignalKind::SoftClip { tail, .. } = &out.signals[1].kind else {
            panic!("not a soft-clip");
        };
        assert_eq!(*tail, None);
    }

    #[test]
    fn short_and_low_quality_clips() {
        let seq = "C".repeat(100);
        let short = read(
            Flags::empty(),
            &[(Kind::SoftClip, 19), (Kind::Match, 81)],
            &seq,
        )
        .build();
        assert_eq!(signals(&short).signals, []);

        let mut qualities = vec![30; 100];
        qualities[..20].fill(10);
        let low = read(
            Flags::empty(),
            &[(Kind::SoftClip, 20), (Kind::Match, 80)],
            &seq,
        )
        .set_quality_scores(QualityScores::from(qualities))
        .build();
        let out = signals(&low);
        assert_eq!(out.signals, []);
        assert_eq!(out.low_quality_clips, 1);
    }

    #[test]
    fn supplementary_tag() {
        let seq = "C".repeat(100);
        let record = read(
            Flags::empty(),
            &[(Kind::SoftClip, 30), (Kind::Match, 70)],
            &seq,
        )
        .set_data(
            [(
                Tag::OTHER_ALIGNMENTS,
                Value::from("chr2,500,+,30M70S,60,0;"),
            )]
            .into_iter()
            .collect(),
        )
        .build();
        let out = signals(&record);
        assert!(matches!(
            out.signals[0].kind,
            SignalKind::SoftClip {
                supplementary: true,
                ..
            }
        ));
    }

    fn pair(flags: Flags) -> noodles::sam::alignment::record_buf::Builder {
        read(
            Flags::SEGMENTED | flags,
            &[(Kind::Match, 100)],
            &"C".repeat(100),
        )
    }

    fn mate_kinds(record: &RecordBuf) -> Vec<(u64, Side, SignalKind)> {
        signals(record)
            .signals
            .into_iter()
            .map(|s| (s.position, s.side, s.kind))
            .collect()
    }

    #[test]
    fn mate_unmapped() {
        let forward = pair(Flags::MATE_UNMAPPED).build();
        assert_eq!(
            mate_kinds(&forward),
            [(1100, Side::Right, SignalKind::MateUnmapped)]
        );
        let reverse = pair(Flags::MATE_UNMAPPED | Flags::REVERSE_COMPLEMENTED).build();
        assert_eq!(
            mate_kinds(&reverse),
            [(1000, Side::Left, SignalKind::MateUnmapped)]
        );
    }

    #[test]
    fn mate_low_mapq_and_elsewhere() {
        let mq = |q: u8| {
            [(Tag::MATE_MAPPING_QUALITY, Value::from(q))]
                .into_iter()
                .collect()
        };
        let low = pair(Flags::empty())
            .set_mate_reference_sequence_id(1)
            .set_mate_alignment_start(Position::try_from(500).unwrap())
            .set_data(mq(0))
            .build();
        assert_eq!(
            mate_kinds(&low),
            [(1100, Side::Right, SignalKind::MateLowMapq { mapq: 0 })]
        );
        let elsewhere = pair(Flags::empty())
            .set_mate_reference_sequence_id(1)
            .set_mate_alignment_start(Position::try_from(500).unwrap())
            .set_data(mq(60))
            .build();
        assert_eq!(
            mate_kinds(&elsewhere),
            [(1100, Side::Right, SignalKind::MateElsewhere { contig: 1 })]
        );
        // A proper pair with a low-MAPQ mate (a repeat next to the read) is
        // not a signal.
        let proper = pair(Flags::PROPERLY_SEGMENTED | Flags::FIRST_SEGMENT)
            .set_mate_reference_sequence_id(0)
            .set_mate_alignment_start(Position::try_from(1200).unwrap())
            .set_template_length(300)
            .set_data(mq(0))
            .build();
        let out = signals(&proper);
        assert_eq!(out.signals, []);
        assert_eq!(out.insert_size, Some(300));
    }

    #[test]
    fn same_contig_pairs() {
        let far = pair(Flags::REVERSE_COMPLEMENTED)
            .set_mate_reference_sequence_id(0)
            .set_mate_alignment_start(Position::try_from(60_000).unwrap())
            .set_template_length(-59_100)
            .build();
        let out = signals(&far);
        assert_eq!(out.signals, []);
        let large = out.large_insert.unwrap();
        assert_eq!(
            (large.position, large.side, large.kind),
            (1000, Side::Left, SignalKind::LargeInsert { length: 59_100 })
        );
        assert_eq!(out.insert_size, None);

        // Read 2 of a proper pair: not sampled for the insert size.
        let second = pair(Flags::PROPERLY_SEGMENTED | Flags::LAST_SEGMENT)
            .set_mate_reference_sequence_id(0)
            .set_mate_alignment_start(Position::try_from(800).unwrap())
            .set_template_length(-300)
            .build();
        assert_eq!(signals(&second).insert_size, None);
    }

    #[test]
    fn unpaired_reads_have_no_mate_signal() {
        let record = read(Flags::empty(), &[(Kind::Match, 100)], &"C".repeat(100)).build();
        let out = signals(&record);
        assert_eq!(out.signals, []);
        assert!(out.large_insert.is_none());
    }
}
