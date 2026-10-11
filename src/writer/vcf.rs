//! VCF 4.2 output of the candidates.
//!
//! One record per candidate, at the base right before the insertion (VCF
//! `POS`, 1-based: the `position` of the TSV), with a symbolic `ALT`
//! (`<INS:ME:ALU>`, `<INS:ME:LINE1>`, `<INS:ME:SVA>`, or `<INS:ME>` for an
//! untyped candidate) and `REF` the reference base there (`N` without
//! `--reference`). `QUAL` is the quality score, `FILTER` the filters of the
//! candidate, and the sample column its genotype. Every `INFO`, `FORMAT`,
//! `FILTER` and `ALT` field is declared in the header.

use std::io::{self, Write};

use super::{Report, candidate_id, confidence, filter};
use crate::cluster::{Candidate, Filter};
use crate::typer::consensus::Family;

/// `FILTER` lines of the header: every filter that can be written.
const FILTERS: [(Filter, &str); 8] = [
    (
        Filter::NoJunction,
        "Fewer soft-clips at the junctions than --min-clips (short ones count only with an MEI mate)",
    ),
    (
        Filter::NoTail,
        "No soft-clip with a poly(A/T) tail (with --require-tail)",
    ),
    (
        Filter::Untyped,
        "Fewer junction soft-clips typed to an MEI consensus than --min-typed-clips",
    ),
    (
        Filter::NonMeiMates,
        "Fewer MEI mates than --min-mei-mate-fraction of the checked mates",
    ),
    (
        Filter::TandemRepeat,
        "Soft-clips explained by a longer reference homopolymer or tandem repeat",
    ),
    (
        Filter::ProbeEdge,
        "Soft-clips at a probe edge without a tail, typed clip or TSD (capture artefact)",
    ),
    (
        Filter::LowRatio,
        "Supporting reads per read of local depth below --min-signal-ratio",
    ),
    (Filter::LowQual, "Quality score below --min-qual"),
];

/// `INFO` lines of the header: ID, number, type, description.
const INFO: [(&str, &str, &str, &str); 19] = [
    ("SVTYPE", "1", "String", "Type of structural variant"),
    (
        "SVLEN",
        "1",
        "Integer",
        "Length of the element consensus covered by the reads (MEINFO end - start + 1): a lower bound of the inserted length",
    ),
    (
        "MEINFO",
        "4",
        "String",
        "Mobile element info of the form NAME,START,END,POLARITY: subfamily, part of its consensus covered by the reads (1-based), strand",
    ),
    (
        "IMPRECISE",
        "0",
        "Flag",
        "No junction assembled: junctions estimated from the soft-clips",
    ),
    (
        "LEFT_JUNCTION",
        "1",
        "Integer",
        "End of the upstream flank (1-based: the insertion is right after this base)",
    ),
    (
        "RIGHT_JUNCTION",
        "1",
        "Integer",
        "Start of the downstream flank, minus one (1-based: the downstream flank starts after this base)",
    ),
    (
        "TSD",
        "1",
        "String",
        "Target site duplication, from the assembled junctions",
    ),
    (
        "TSD_LEN",
        "1",
        "Integer",
        "Length of the target site duplication (assembled, or estimated from the soft-clips)",
    ),
    (
        "POLYA_LEN",
        "1",
        "Integer",
        "Length of the poly(A) (or poly(T)) tail, from an assembled junction",
    ),
    (
        "ASSEMBLED",
        "1",
        "String",
        "Junctions assembled: both, left, right or none",
    ),
    (
        "CONFIDENCE",
        "1",
        "String",
        "Confidence of a PASS candidate: high (tail and TSD), medium (tail), low (no tail)",
    ),
    ("SUPPORT", "1", "Integer", "Distinct supporting reads"),
    (
        "JUNCTION_CLIPS",
        "2",
        "Integer",
        "Soft-clips at the left and right junctions",
    ),
    (
        "TAIL_CLIPS",
        "2",
        "Integer",
        "Soft-clips with a poly(A) and with a poly(T) tail",
    ),
    (
        "TYPED_CLIPS",
        "1",
        "Integer",
        "Junction soft-clips typed to the element family",
    ),
    (
        "MATES",
        "2",
        "Integer",
        "Reads with an unmapped mate, and with a discordant mate (low MAPQ, another contig, large insert)",
    ),
    (
        "MEI_MATES",
        "2",
        "Integer",
        "Mates checked against the consensus sequences, and of these mates of the element family",
    ),
    (
        "DISTINCT_STARTS",
        "1",
        "Integer",
        "Distinct fragment start positions of the supporting reads",
    ),
    (
        "SIGNAL_RATIO",
        "1",
        "Float",
        "Supporting reads per read of local depth",
    ),
];

/// `FORMAT` lines of the header.
const FORMAT: [(&str, &str, &str, &str); 4] = [
    ("GT", "1", "String", "Genotype"),
    (
        "GQ",
        "1",
        "Integer",
        "Genotype quality: phred-scaled probability that the genotype is wrong",
    ),
    (
        "AD",
        "R",
        "Integer",
        "Reads of the reference allele (aligned over the whole site) and of the insertion (clipped at a junction)",
    ),
    ("DP", "1", "Integer", "Reads counted for the genotype"),
];

/// Symbolic `ALT` of a candidate.
fn alt(candidate: &Candidate) -> &'static str {
    match candidate.typing.as_ref().map(|t| t.family) {
        Some(Family::Alu) => "<INS:ME:ALU>",
        Some(Family::Line1) => "<INS:ME:LINE1>",
        Some(Family::Sva) => "<INS:ME:SVA>",
        None => "<INS:ME>",
    }
}

/// Write the report as VCF 4.2.
///
/// # Errors
///
/// Fails if writing fails.
pub fn write_vcf(report: &Report<'_>, out: &mut dyn Write) -> io::Result<()> {
    writeln!(out, "##fileformat=VCFv4.2")?;
    writeln!(out, "##source=mei-rs {}", env!("CARGO_PKG_VERSION"))?;
    writeln!(out, "##mei-rsCommand={}", report.command)?;
    if let Some(reference) = report.reference {
        writeln!(out, "##reference={reference}")?;
    }
    for (name, length) in report.contigs.iter().zip(report.contig_lengths) {
        writeln!(out, "##contig=<ID={name},length={length}>")?;
    }
    for (id, element) in [
        ("INS:ME:ALU", "an Alu element"),
        ("INS:ME:LINE1", "a LINE-1 element"),
        ("INS:ME:SVA", "an SVA element"),
        ("INS:ME", "a mobile element of unknown family"),
    ] {
        writeln!(
            out,
            "##ALT=<ID={id},Description=\"Insertion of {element}\">"
        )?;
    }
    writeln!(out, "##FILTER=<ID=PASS,Description=\"All filters passed\">")?;
    for (filter, description) in FILTERS {
        writeln!(
            out,
            "##FILTER=<ID={},Description=\"{description}\">",
            filter.name()
        )?;
    }
    for (id, number, kind, description) in INFO {
        writeln!(
            out,
            "##INFO=<ID={id},Number={number},Type={kind},Description=\"{description}\">"
        )?;
    }
    for (id, number, kind, description) in FORMAT {
        writeln!(
            out,
            "##FORMAT=<ID={id},Number={number},Type={kind},Description=\"{description}\">"
        )?;
    }
    writeln!(
        out,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
        report.sample
    )?;
    for (i, c) in report.candidates.iter().enumerate() {
        let reference = report
            .reference_bases
            .and_then(|bases| bases.get(i).copied())
            .filter(u8::is_ascii_alphabetic)
            .map_or('N', |b| char::from(b.to_ascii_uppercase()));
        writeln!(
            out,
            "{}\t{}\t{}\t{reference}\t{}\t{}\t{}\t{}\tGT:GQ:AD:DP\t{}",
            report.contigs[c.contig],
            c.breakpoint.max(1),
            candidate_id(i),
            alt(c),
            c.qual.map_or_else(|| ".".to_owned(), |q| q.to_string()),
            filter(c),
            info(c).join(";"),
            sample(c)
        )?;
    }
    Ok(())
}

/// `INFO` fields of a candidate.
fn info(c: &Candidate) -> Vec<String> {
    let mut fields = vec!["SVTYPE=INS".to_owned()];
    if let Some(t) = &c.typing {
        fields.push(format!("SVLEN={}", t.element_end - t.element_start));
        fields.push(format!(
            "MEINFO={},{},{},{}",
            t.subfamily,
            t.element_start + 1,
            t.element_end,
            if t.reverse { '-' } else { '+' }
        ));
    }
    let assembled = c.assembly.as_ref();
    if assembled.is_none_or(|a| a.left.is_none() && a.right.is_none()) {
        fields.push("IMPRECISE".to_owned());
    }
    if let Some(left) = c.left_junction {
        fields.push(format!("LEFT_JUNCTION={left}"));
    }
    if let Some(right) = c.right_junction {
        fields.push(format!("RIGHT_JUNCTION={right}"));
    }
    if let Some(tsd) = assembled.and_then(|a| a.tsd.as_deref()) {
        fields.push(format!("TSD={}", String::from_utf8_lossy(tsd)));
    }
    if let Some(length) = c.tsd_length() {
        fields.push(format!("TSD_LEN={length}"));
    }
    if let Some(tail) = assembled.and_then(|a| a.tail) {
        fields.push(format!("POLYA_LEN={}", tail.length));
    }
    if let Some(a) = assembled {
        fields.push(format!("ASSEMBLED={}", a.assembled()));
    }
    if c.is_pass() {
        fields.push(format!("CONFIDENCE={}", confidence(c)));
    }
    fields.extend([
        format!("SUPPORT={}", c.support),
        format!(
            "JUNCTION_CLIPS={},{}",
            c.left_junction_clips, c.right_junction_clips
        ),
        format!("TAIL_CLIPS={},{}", c.polya_clips, c.polyt_clips),
        format!(
            "TYPED_CLIPS={}",
            c.typing.as_ref().map_or(0, |t| t.typed_clips)
        ),
        format!("MATES={},{}", c.unmapped_mates, c.discordant_mates),
        format!("MEI_MATES={},{}", c.checked_mates, c.mei_mates),
        format!("DISTINCT_STARTS={}", c.distinct_starts),
    ]);
    if let Some(ratio) = c.signal_ratio() {
        fields.push(format!("SIGNAL_RATIO={:.3}", super::ratio(ratio)));
    }
    fields
}

/// Sample column (`GT:GQ:AD:DP`) of a candidate.
fn sample(c: &Candidate) -> String {
    let Some(g) = c.genotype else {
        return "./.:.:.:.".to_owned();
    };
    // Reference reads are counted at each junction: once per site here.
    let junctions = u32::from(c.left_junction.is_some()) + u32::from(c.right_junction.is_some());
    let reference = g.reference / junctions.max(1);
    format!(
        "{}:{}:{reference},{}:{}",
        g.call.name(),
        g.quality,
        g.alternative,
        reference + g.alternative
    )
}
