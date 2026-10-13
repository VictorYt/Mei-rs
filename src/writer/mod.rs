//! Export of the candidates as TSV, BED, JSON or VCF ([`vcf`]).
//!
//! Positions in TSV, JSON and VCF are 1-based: a position `p` means the
//! insertion lies right after base `p`. The BED output is 0-based
//! half-open.

pub mod vcf;

use std::io::{self, Write};

use serde::Serialize;

use crate::assembly::Assembly;
use crate::cli::OutputFormat;
use crate::cluster::{Candidate, Confidence};
use crate::genotype::Genotype;
use crate::reader::signals::{Side, Signal, SignalKind};
use crate::typer::Typing;

/// What the outputs describe.
#[derive(Clone, Copy, Debug)]
pub struct Report<'a> {
    /// Sample name.
    pub sample: &'a str,
    /// Command line, for the record.
    pub command: &'a str,
    /// Contig names, by header index.
    pub contigs: &'a [String],
    /// Contig lengths, by header index.
    pub contig_lengths: &'a [u64],
    /// Name of the reference file, if any.
    pub reference: Option<&'a str>,
    /// Reference base at the position of each candidate (VCF `REF`), with
    /// the reference.
    pub reference_bases: Option<&'a [u8]>,
    /// Candidates, sorted by contig and position.
    pub candidates: &'a [Candidate],
}

/// Write the report in the given format.
///
/// # Errors
///
/// Fails if writing fails.
pub fn write(report: &Report<'_>, format: OutputFormat, out: &mut dyn Write) -> io::Result<()> {
    match format {
        OutputFormat::Tsv => write_tsv(report, out),
        OutputFormat::Bed => write_bed(report, out),
        OutputFormat::Json => write_json(report, out),
        OutputFormat::Vcf => vcf::write_vcf(report, out),
    }
}

/// Columns of the TSV output.
pub const TSV_COLUMNS: [&str; 46] = [
    "id",
    "contig",
    "position",
    "left_junction",
    "right_junction",
    "tsd_length",
    "strand",
    "filter",
    "confidence",
    "support",
    "left_clips",
    "right_clips",
    "left_junction_clips",
    "right_junction_clips",
    "short_junction_clips",
    "unmapped_mates",
    "discordant_mates",
    "polya_clips",
    "polyt_clips",
    "both_sides",
    "distinct_starts",
    "depth",
    "signal_ratio",
    "start",
    "end",
    "family",
    "subfamily",
    "element_strand",
    "element_start",
    "element_end",
    "identity",
    "typed_clips",
    "short_typed_clips",
    "checked_mates",
    "mei_mates",
    "assembled",
    "tsd",
    "polya_length",
    "qual",
    "genotype",
    "gq",
    "alt_reads",
    "ref_reads",
    "typed_by",
    "cipos",
    "breakpoint_clips",
];

fn write_tsv(report: &Report<'_>, out: &mut dyn Write) -> io::Result<()> {
    writeln!(out, "##mei-rs {}", env!("CARGO_PKG_VERSION"))?;
    writeln!(out, "##command: {}", report.command)?;
    writeln!(out, "##sample: {}", report.sample)?;
    writeln!(
        out,
        "##positions: 1-based; a position p means the insertion is right after base p"
    )?;
    writeln!(out, "#{}", TSV_COLUMNS.join("\t"))?;
    let opt = |v: Option<u64>| v.map_or_else(|| ".".to_owned(), |v| v.to_string());
    for (i, c) in report.candidates.iter().enumerate() {
        let fields = [
            candidate_id(i),
            report.contigs[c.contig].clone(),
            c.breakpoint.to_string(),
            opt(c.left_junction),
            opt(c.right_junction),
            opt(c.tsd_length()),
            strand(c).to_owned(),
            filter(c),
            confidence(c).to_owned(),
            c.support.to_string(),
            c.left_clips.to_string(),
            c.right_clips.to_string(),
            c.left_junction_clips.to_string(),
            c.right_junction_clips.to_string(),
            c.short_junction_clips.to_string(),
            c.unmapped_mates.to_string(),
            c.discordant_mates.to_string(),
            c.polya_clips.to_string(),
            c.polyt_clips.to_string(),
            if c.both_sides() { "yes" } else { "no" }.to_owned(),
            c.distinct_starts.to_string(),
            opt(c.depth),
            c.signal_ratio()
                .map_or_else(|| ".".to_owned(), |r| format!("{:.3}", ratio(r))),
            c.start.to_string(),
            c.end.to_string(),
        ];
        let typing = c.typing.as_ref();
        let assembly = c.assembly.as_ref();
        let genotyped =
            |f: &dyn Fn(&Genotype) -> String| c.genotype.as_ref().map_or_else(|| ".".to_owned(), f);
        let typed = |f: &dyn Fn(&Typing) -> String| typing.map_or_else(|| ".".to_owned(), f);
        let typing_fields = [
            typed(&|t| t.family.name().to_owned()),
            typed(&|t| t.subfamily.clone()),
            typed(&|t| element_strand(t).to_owned()),
            typed(&|t| (t.element_start + 1).to_string()),
            typed(&|t| t.element_end.to_string()),
            typed(&|t| format!("{:.3}", ratio(t.identity))),
            typing.map_or(0, |t| t.typed_clips).to_string(),
            typing.map_or(0, |t| t.short_typed_clips).to_string(),
            c.checked_mates.to_string(),
            c.mei_mates.to_string(),
            assembly.map_or(".", Assembly::assembled).to_owned(),
            assembly.and_then(|a| a.tsd.as_deref()).map_or_else(
                || ".".to_owned(),
                |t| String::from_utf8_lossy(t).into_owned(),
            ),
            assembly
                .and_then(|a| a.tail)
                .map_or_else(|| ".".to_owned(), |t| t.length.to_string()),
            c.qual.map_or_else(|| ".".to_owned(), |q| q.to_string()),
            genotyped(&|g| g.call.name().to_owned()),
            genotyped(&|g| g.quality.to_string()),
            genotyped(&|g| g.alternative.to_string()),
            genotyped(&|g| g.reference.to_string()),
            typed(&|t| typed_by(t).to_owned()),
            cipos(c).map_or_else(|| ".".to_owned(), |(a, b)| format!("{a},{b}")),
            c.breakpoint_clips.to_string(),
        ];
        writeln!(out, "{}\t{}", fields.join("\t"), typing_fields.join("\t"))?;
    }
    Ok(())
}

fn write_bed(report: &Report<'_>, out: &mut dyn Write) -> io::Result<()> {
    for (i, c) in report.candidates.iter().enumerate() {
        // The TSD when the junctions show one, else the base before the
        // insertion.
        let (start, end) = match (c.tsd_length(), c.right_junction) {
            (Some(length), Some(right)) => (right, right + length),
            _ => (c.breakpoint.saturating_sub(1), c.breakpoint.max(1)),
        };
        writeln!(
            out,
            "{}\t{start}\t{end}\t{}\t{}\t{}\t{}",
            report.contigs[c.contig],
            candidate_id(i),
            c.support.min(1000),
            strand(c),
            if c.is_pass() {
                confidence(c).to_owned()
            } else {
                filter(c)
            }
        )?;
    }
    Ok(())
}

fn write_json(report: &Report<'_>, out: &mut dyn Write) -> io::Result<()> {
    let json = JsonReport {
        mei_rs: env!("CARGO_PKG_VERSION"),
        command: report.command,
        sample: report.sample,
        positions: "1-based; a position p means the insertion is right after base p",
        candidates: report
            .candidates
            .iter()
            .enumerate()
            .map(|(i, c)| JsonCandidate::new(i, c, report.contigs))
            .collect(),
    };
    serde_json::to_writer_pretty(&mut *out, &json)?;
    writeln!(out)
}

/// Identifier of the `i`-th candidate of the output.
#[must_use]
pub fn candidate_id(i: usize) -> String {
    format!("MEI_{}", i + 1)
}

/// `PASS`, or the failed filters separated by `;`.
fn filter(candidate: &Candidate) -> String {
    if candidate.is_pass() {
        "PASS".to_owned()
    } else {
        candidate
            .filters
            .iter()
            .map(|f| f.name())
            .collect::<Vec<_>>()
            .join(";")
    }
}

/// Confidence level, `.` for a filtered candidate.
fn confidence(candidate: &Candidate) -> &'static str {
    candidate.confidence().map_or(".", Confidence::name)
}

/// Insertion strand from the tails: poly(A) on the reference strand for a
/// plus-strand element, poly(T) for a minus-strand one.
fn strand(candidate: &Candidate) -> &'static str {
    match candidate.polya_clips.cmp(&candidate.polyt_clips) {
        std::cmp::Ordering::Greater => "+",
        std::cmp::Ordering::Less => "-",
        std::cmp::Ordering::Equal => ".",
    }
}

/// Strand of the typed element: `-` when the clips match the reverse
/// complement of the consensus.
fn element_strand(typing: &Typing) -> &'static str {
    if typing.reverse { "-" } else { "+" }
}

/// Interval of the insertion placed by its mates, relative to its position
/// (as the VCF `CIPOS`).
fn cipos(candidate: &Candidate) -> Option<(i64, i64)> {
    #[allow(clippy::cast_possible_wrap)]
    let offset = |b: u64| b as i64 - candidate.breakpoint as i64;
    candidate
        .mate_interval
        .map(|(first, last)| (offset(first), offset(last)))
}

/// What typed the element: `clips` (junction soft-clips) or `mates` (no clip
/// typed).
fn typed_by(typing: &Typing) -> &'static str {
    if typing.by_mates { "mates" } else { "clips" }
}

fn ratio(value: f64) -> f64 {
    (value * 1000.0).round() / 1000.0
}

#[derive(Serialize)]
struct JsonReport<'a> {
    mei_rs: &'a str,
    command: &'a str,
    sample: &'a str,
    positions: &'a str,
    candidates: Vec<JsonCandidate>,
}

/// Assembly of a candidate, in JSON.
#[derive(Serialize)]
struct JsonAssembly {
    assembled: &'static str,
    #[serde(skip_serializing_if = "Option::is_none")]
    left: Option<JsonContig>,
    #[serde(skip_serializing_if = "Option::is_none")]
    right: Option<JsonContig>,
    tsd: Option<String>,
    polya_length: Option<u32>,
    tail_base: Option<String>,
}

/// An assembled junction, in JSON (junction 1-based: the last base of the
/// upstream flank, or the base before the downstream flank).
#[derive(Serialize)]
struct JsonContig {
    junction: u64,
    contig: String,
    inserted: String,
}

impl JsonAssembly {
    fn new(a: &Assembly) -> Self {
        let text = |s: &[u8]| String::from_utf8_lossy(s).into_owned();
        Self {
            assembled: a.assembled(),
            left: a.left.as_ref().map(|l| JsonContig {
                junction: l.junction,
                contig: text(&l.contig),
                inserted: text(l.left_insert()),
            }),
            right: a.right.as_ref().map(|r| JsonContig {
                junction: r.junction,
                contig: text(&r.contig),
                inserted: text(r.right_insert()),
            }),
            tsd: a.tsd.as_deref().map(text),
            polya_length: a.tail.map(|t| t.length),
            tail_base: a.tail.map(|t| char::from(t.base).to_string()),
        }
    }
}

#[derive(Serialize)]
struct JsonCandidate {
    id: String,
    contig: String,
    position: u64,
    left_junction: Option<u64>,
    right_junction: Option<u64>,
    tsd_length: Option<u64>,
    strand: &'static str,
    filter: Vec<&'static str>,
    confidence: Option<&'static str>,
    support: usize,
    left_clips: usize,
    right_clips: usize,
    left_junction_clips: usize,
    right_junction_clips: usize,
    short_junction_clips: usize,
    unmapped_mates: usize,
    discordant_mates: usize,
    polya_clips: usize,
    polyt_clips: usize,
    both_sides: bool,
    distinct_starts: usize,
    depth: Option<u64>,
    signal_ratio: Option<f64>,
    start: u64,
    end: u64,
    family: Option<&'static str>,
    subfamily: Option<String>,
    element_strand: Option<&'static str>,
    element_start: Option<usize>,
    element_end: Option<usize>,
    identity: Option<f64>,
    typed_by: Option<&'static str>,
    typed_clips: usize,
    cipos: Option<(i64, i64)>,
    breakpoint_clips: usize,
    short_typed_clips: usize,
    checked_mates: usize,
    mei_mates: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    assembly: Option<JsonAssembly>,
    qual: Option<u32>,
    genotype: Option<&'static str>,
    gq: Option<u32>,
    alt_reads: Option<u32>,
    ref_reads: Option<u32>,
    signals: Vec<JsonSignal>,
}

impl JsonCandidate {
    fn new(i: usize, c: &Candidate, contigs: &[String]) -> Self {
        Self {
            id: candidate_id(i),
            contig: contigs[c.contig].clone(),
            position: c.breakpoint,
            left_junction: c.left_junction,
            right_junction: c.right_junction,
            tsd_length: c.tsd_length(),
            strand: strand(c),
            filter: if c.is_pass() {
                vec!["PASS"]
            } else {
                c.filters.iter().map(|f| f.name()).collect()
            },
            confidence: c.confidence().map(Confidence::name),
            support: c.support,
            left_clips: c.left_clips,
            right_clips: c.right_clips,
            left_junction_clips: c.left_junction_clips,
            right_junction_clips: c.right_junction_clips,
            short_junction_clips: c.short_junction_clips,
            unmapped_mates: c.unmapped_mates,
            discordant_mates: c.discordant_mates,
            polya_clips: c.polya_clips,
            polyt_clips: c.polyt_clips,
            both_sides: c.both_sides(),
            distinct_starts: c.distinct_starts,
            depth: c.depth,
            signal_ratio: c.signal_ratio().map(ratio),
            start: c.start,
            end: c.end,
            family: c.typing.as_ref().map(|t| t.family.name()),
            subfamily: c.typing.as_ref().map(|t| t.subfamily.clone()),
            element_strand: c.typing.as_ref().map(element_strand),
            element_start: c.typing.as_ref().map(|t| t.element_start + 1),
            element_end: c.typing.as_ref().map(|t| t.element_end),
            identity: c.typing.as_ref().map(|t| ratio(t.identity)),
            typed_by: c.typing.as_ref().map(typed_by),
            cipos: cipos(c),
            breakpoint_clips: c.breakpoint_clips,
            typed_clips: c.typing.as_ref().map_or(0, |t| t.typed_clips),
            short_typed_clips: c.typing.as_ref().map_or(0, |t| t.short_typed_clips),
            checked_mates: c.checked_mates,
            mei_mates: c.mei_mates,
            assembly: c.assembly.as_ref().map(JsonAssembly::new),
            qual: c.qual,
            genotype: c.genotype.map(|g| g.call.name()),
            gq: c.genotype.map(|g| g.quality),
            alt_reads: c.genotype.map(|g| g.alternative),
            ref_reads: c.genotype.map(|g| g.reference),
            signals: c
                .signals
                .iter()
                .enumerate()
                .map(|(i, s)| {
                    let mei = c
                        .mate_checks
                        .iter()
                        .find(|m| m.signal == i)
                        .and_then(|m| m.mei);
                    JsonSignal::new(s, contigs, mei)
                })
                .collect(),
        }
    }
}

#[derive(Serialize)]
struct JsonSignal {
    kind: &'static str,
    position: u64,
    side: &'static str,
    read: String,
    segment: u8,
    strand: &'static str,
    /// Aligned bases of the read, 1-based inclusive.
    read_start: u64,
    read_end: u64,
    mapq: Option<u8>,
    #[serde(skip_serializing_if = "Option::is_none")]
    clipped_sequence: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    median_quality: Option<u8>,
    #[serde(skip_serializing_if = "Option::is_none")]
    tail: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    supplementary: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    short: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    mate_mapq: Option<u8>,
    #[serde(skip_serializing_if = "Option::is_none")]
    mate_contig: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    insert_length: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    mei_mate: Option<bool>,
}

impl JsonSignal {
    fn new(s: &Signal, contigs: &[String], mei_mate: Option<bool>) -> Self {
        let mut json = Self {
            kind: s.kind.name(),
            position: s.position,
            side: match s.side {
                Side::Left => "left",
                Side::Right => "right",
            },
            read: s.read.clone(),
            segment: if s.first_segment { 1 } else { 2 },
            strand: if s.reverse { "-" } else { "+" },
            read_start: s.start + 1,
            read_end: s.end,
            mapq: s.mapq,
            clipped_sequence: None,
            median_quality: None,
            tail: None,
            supplementary: None,
            short: None,
            mate_mapq: None,
            mate_contig: None,
            insert_length: None,
            mei_mate,
        };
        match &s.kind {
            SignalKind::SoftClip {
                sequence,
                median_quality,
                tail,
                supplementary,
                short,
            } => {
                json.clipped_sequence = Some(String::from_utf8_lossy(sequence).into_owned());
                json.median_quality = Some(*median_quality);
                json.tail = tail.map(|t| format!("{}{}", char::from(t.base), t.length));
                json.supplementary = Some(*supplementary);
                json.short = Some(*short);
            }
            SignalKind::MateUnmapped => {}
            SignalKind::MateLowMapq { mapq, .. } => json.mate_mapq = Some(*mapq),
            SignalKind::MateElsewhere { mate } => {
                json.mate_contig = Some(contigs[mate.contig()].clone());
            }
            SignalKind::LargeInsert { length, .. } => json.insert_length = Some(*length),
        }
        json
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cluster::tests::{clip, mate};
    use crate::cluster::{ClusterOptions, Filter, cluster};

    fn candidates() -> Vec<Candidate> {
        let mut signals = Vec::new();
        for i in 0..3 {
            signals.push(clip(1_012, Side::Right, &format!("r{i}"), None));
            signals.push(clip(1_000, Side::Left, &format!("l{i}"), Some(b'A')));
        }
        signals.push(mate(
            1_050,
            Side::Left,
            "m",
            SignalKind::MateElsewhere {
                mate: crate::reader::signals::Locus {
                    contig: 1,
                    span: None,
                    start: 99,
                },
            },
        ));
        signals.push(clip(5_000, Side::Left, "e1", None));
        signals.push(clip(5_000, Side::Left, "e2", None));
        signals.push(clip(5_000, Side::Left, "e3", None));
        signals.sort();
        let mut candidates = cluster(signals, &ClusterOptions::default());
        candidates[0].depth = Some(40);
        candidates[0].typing = Some(Typing {
            family: crate::typer::consensus::Family::Alu,
            subfamily: "AluYa5".to_owned(),
            reverse: false,
            element_start: 0,
            element_end: 281,
            identity: 0.98765,
            typed_clips: 6,
            short_typed_clips: 0,
            by_mates: false,
        });
        candidates[1].depth = Some(100);
        candidates[1].filters = vec![Filter::ProbeEdge, Filter::LowRatio];
        candidates
    }

    fn render(format: OutputFormat) -> String {
        let candidates = candidates();
        let contigs = ["chr1".to_owned(), "chr2".to_owned()];
        let report = Report {
            sample: "s1",
            command: "mei-rs scan",
            contigs: &contigs,
            contig_lengths: &[100_000],
            reference: None,
            reference_bases: None,
            candidates: &candidates,
        };
        let mut out = Vec::new();
        write(&report, format, &mut out).unwrap();
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn tsv() {
        let tsv = render(OutputFormat::Tsv);
        let lines: Vec<&str> = tsv.lines().collect();
        assert_eq!(lines[1], "##command: mei-rs scan");
        assert_eq!(lines[2], "##sample: s1");
        assert_eq!(lines[4], format!("#{}", TSV_COLUMNS.join("\t")));
        assert_eq!(
            lines[5],
            "MEI_1\tchr1\t1000\t1012\t1000\t12\t+\tPASS\thigh\t7\t3\t3\t3\t3\t0\t0\t1\t3\t0\tyes\t3\t40\t0.175\t1000\t1050\tALU\tAluYa5\t+\t1\t281\t0.988\t6\t0\t0\t0\t.\t.\t.\t.\t.\t.\t.\t.\tclips\t.\t0"
        );
        assert_eq!(
            lines[6],
            "MEI_2\tchr1\t5000\t.\t5000\t.\t.\tprobe_edge;low_ratio\t.\t3\t3\t0\t0\t3\t0\t0\t0\t0\t0\tno\t1\t100\t0.030\t5000\t5000\t.\t.\t.\t.\t.\t.\t0\t0\t0\t0\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t0"
        );
        assert_eq!(lines.len(), 7);
        for line in &lines[4..] {
            assert_eq!(line.split('\t').count(), TSV_COLUMNS.len());
        }
    }

    #[test]
    fn bed() {
        assert_eq!(
            render(OutputFormat::Bed),
            "chr1\t1000\t1012\tMEI_1\t7\t+\thigh\nchr1\t4999\t5000\tMEI_2\t3\t.\tprobe_edge;low_ratio\n"
        );
    }

    #[test]
    fn json() {
        let json: serde_json::Value = serde_json::from_str(&render(OutputFormat::Json)).unwrap();
        assert_eq!(json["sample"], "s1");
        let first = &json["candidates"][0];
        assert_eq!(first["id"], "MEI_1");
        assert_eq!(first["filter"], serde_json::json!(["PASS"]));
        assert_eq!(first["tsd_length"], 12);
        assert_eq!(first["confidence"], "high");
        assert_eq!(first["left_junction_clips"], 3);
        assert!(json["candidates"][1]["confidence"].is_null());
        assert_eq!(first["signal_ratio"], 0.175);
        assert_eq!(first["family"], "ALU");
        assert_eq!(first["subfamily"], "AluYa5");
        assert_eq!(first["element_strand"], "+");
        assert_eq!(
            (first["element_start"].clone(), first["element_end"].clone()),
            (1.into(), 281.into())
        );
        assert_eq!(first["identity"], 0.988);
        assert_eq!(first["typed_clips"], 6);
        assert!(json["candidates"][1]["family"].is_null());
        assert_eq!(json["candidates"][1]["typed_clips"], 0);
        let signals = first["signals"].as_array().unwrap();
        assert_eq!(signals.len(), 7);
        let clip = &signals[0];
        assert_eq!(clip["kind"], "soft_clip");
        assert_eq!(clip["tail"], "A20");
        assert_eq!(clip["read_start"], 1_001);
        assert!(clip.get("mate_contig").is_none());
        let mate = signals
            .iter()
            .find(|s| s["kind"] == "mate_elsewhere")
            .unwrap();
        assert_eq!(mate["mate_contig"], "chr2");
        assert!(mate.get("clipped_sequence").is_none());
        assert_eq!(
            json["candidates"][1]["filter"],
            serde_json::json!(["probe_edge", "low_ratio"])
        );
    }
}
