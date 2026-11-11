#![no_main]
//! The four output formats of candidates built from the input (see
//! `common`): written without error, one TSV row and one VCF record per
//! candidate, and the VCF read back by noodles.

mod common;

use std::io::BufReader;

use libfuzzer_sys::fuzz_target;
use mei_rs::cli::OutputFormat;
use mei_rs::writer::{self, Report};
use noodles::vcf;

fuzz_target!(|data: &[u8]| {
    let set = common::consensus();
    let mut input = common::Input(data);
    let anonymize = input.u8() & 1 != 0;
    let candidates = common::candidates(&mut input, set);
    let contigs = ["chr1".to_owned()];
    let reference_bases: Vec<u8> = candidates.iter().map(|_| b'A').collect();
    let report = Report {
        sample: "fuzz",
        command: "mei-rs scan",
        contigs: &contigs,
        contig_lengths: &[common::CONTIG_LENGTH],
        reference: Some("reference.fa"),
        reference_bases: Some(&reference_bases),
        candidates: &candidates,
        anonymize,
    };
    let render = |format| {
        let mut out = Vec::new();
        writer::write(&report, format, &mut out).expect("writing to memory");
        String::from_utf8(out).expect("UTF-8 output")
    };
    let tsv = render(OutputFormat::Tsv);
    assert_eq!(
        tsv.lines().filter(|l| !l.starts_with('#')).count(),
        candidates.len()
    );
    render(OutputFormat::Bed);
    render(OutputFormat::Json);

    let text = render(OutputFormat::Vcf);
    let mut reader = vcf::io::Reader::new(BufReader::new(text.as_bytes()));
    let header = reader.read_header().expect("a valid VCF header");
    let records: Vec<_> = reader
        .record_bufs(&header)
        .collect::<Result<_, _>>()
        .expect("valid VCF records");
    assert_eq!(records.len(), candidates.len());
});
