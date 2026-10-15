#![allow(missing_docs)] // `criterion_group!` generates an undocumented function.
//! Benchmarks of the stages of `scan` on the synthetic fixtures (extraction,
//! clustering, typing) and of the routines they spend their time in. Run
//! with `cargo bench`; put the results of the affected benchmarks in the
//! pull request of every performance change.

use std::fmt::Write;
use std::hint::black_box;
use std::io::Cursor;
use std::path::{Path, PathBuf};

use criterion::{Criterion, criterion_group, criterion_main};

use mei_rs::cluster::{ClusterOptions, cluster};
use mei_rs::reader::alignment::AlignmentFile;
use mei_rs::reader::extract::{ExtractOptions, extract};
use mei_rs::reader::signals::{ClipOptions, poly_tail};
use mei_rs::regions::bed;
use mei_rs::regions::targets::TargetSet;
use mei_rs::typer::align;
use mei_rs::typer::consensus::ConsensusSet;
use mei_rs::typer::{type_candidate, type_sequence};

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn stages(c: &mut Criterion) {
    let reference = data("reference.fa");
    let file = AlignmentFile::open(&data("positive.bam"), Some(&reference)).unwrap();
    let probes = bed::read(&data("targets.bed")).unwrap();
    let targets = TargetSet::new(probes, 300, &file.contigs()).unwrap();
    let options = ExtractOptions::default();
    let set = ConsensusSet::load(None).unwrap();

    c.bench_function("extract/positive.bam", |b| {
        b.iter(|| extract(black_box(&file), black_box(&targets), &options).unwrap());
    });

    let extraction = extract(&file, &targets, &options).unwrap();
    let cluster_options = ClusterOptions::default();
    c.bench_function("cluster/positive.bam", |b| {
        b.iter(|| cluster(black_box(extraction.signals.clone()), &cluster_options));
    });

    let candidates = cluster(extraction.signals.clone(), &cluster_options);
    c.bench_function("type_candidates/positive.bam", |b| {
        b.iter(|| {
            for candidate in &candidates {
                black_box(type_candidate(black_box(&set), candidate));
            }
        });
    });
}

fn routines(c: &mut Criterion) {
    let set = ConsensusSet::load(None).unwrap();
    let alu = set
        .sequences()
        .iter()
        .find(|s| s.name.starts_with("AluY"))
        .expect("an AluY consensus");
    let clip = &alu.sequence[40..190];

    c.bench_function("type_sequence/150 bp Alu clip", |b| {
        b.iter(|| type_sequence(black_box(&set), black_box(clip)));
    });
    c.bench_function("align_local/150 x 300 bp", |b| {
        b.iter(|| align::local(black_box(clip), black_box(&alu.sequence[..300])));
    });
    let tailed: Vec<u8> = [clip, &[b'A'; 30][..]].concat();
    c.bench_function("poly_tail/180 bp clip", |b| {
        b.iter(|| poly_tail(black_box(&tailed), false, &ClipOptions::default()));
    });
    let mut bed_text = String::new();
    for i in 0..50_000 {
        writeln!(
            bed_text,
            "chr{}\t{}\t{}",
            i % 22 + 1,
            i * 500,
            i * 500 + 150
        )
        .unwrap();
    }
    c.bench_function("bed_parse/50,000 probes", |b| {
        b.iter(|| {
            bed::parse(
                Cursor::new(black_box(bed_text.as_bytes())),
                Path::new("a.bed"),
            )
            .unwrap()
        });
    });
}

criterion_group!(benches, stages, routines);
criterion_main!(benches);
