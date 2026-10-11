//! End-to-end tests of `mei-rs scan` on the synthetic fixtures
//! (see tests/data/README.md).
//!
//! The expected outputs are in tests/data/expected; regenerate them with
//! `UPDATE_EXPECTED=1 cargo test --test scan` after an intended change.

use std::fs;
use std::path::Path;

use assert_cmd::Command;
use predicates::prelude::*;

const ROOT: &str = env!("CARGO_MANIFEST_DIR");

/// Run `mei-rs scan` from the repository root, quietly, and return stdout.
fn scan(args: &[&str]) -> String {
    let output = Command::cargo_bin("mei-rs")
        .unwrap()
        .current_dir(ROOT)
        .args(["-q", "scan"])
        .args(args)
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    String::from_utf8(output).unwrap()
}

/// TSV rows (without headers), as columns.
fn rows(tsv: &str) -> Vec<Vec<String>> {
    tsv.lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect()
}

/// Value of a column of a TSV row.
fn column<'a>(tsv: &str, row: &'a [String], name: &str) -> &'a str {
    let header = tsv.lines().find(|l| l.starts_with("#id")).unwrap();
    let index = header[1..].split('\t').position(|c| c == name).unwrap();
    &row[index]
}

/// Remove what depends on the version and on the command line.
fn normalize(output: &str) -> String {
    output
        .lines()
        .map(|line| {
            let trimmed = line.trim_start();
            if line.starts_with("##mei-rs ") {
                "##mei-rs {version}".to_owned()
            } else if line.starts_with("##command: ") {
                "##command: {command}".to_owned()
            } else if trimmed.starts_with("\"mei_rs\": ") {
                line.replace(env!("CARGO_PKG_VERSION"), "{version}")
            } else if trimmed.starts_with("\"command\": ") {
                format!(
                    "{}\"command\": \"{{command}}\",",
                    &line[..line.len() - trimmed.len()]
                )
            } else {
                line.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n")
}

const POSITIVE: [&str; 4] = [
    "--bam",
    "tests/data/positive.bam",
    "--bed",
    "tests/data/targets.bed",
];

#[test]
fn positive_sample_has_one_pass_candidate_per_insertion() {
    let tsv = scan(&POSITIVE);
    let rows = rows(&tsv);
    // Truth: tests/data/truth.tsv (SVA, Alu and LINE-1 insertions, each
    // between the copies of its TSD), in position order.
    let truth = fs::read_to_string(Path::new(ROOT).join("tests/data/truth.tsv")).unwrap();
    let truth: Vec<Vec<&str>> = truth
        .lines()
        .skip(1)
        .map(|l| l.split('\t').collect())
        .collect();
    assert_eq!(rows.len(), truth.len(), "{tsv}");
    for (row, truth) in rows.iter().zip(&truth) {
        let (left, right): (u64, u64) = (truth[1].parse().unwrap(), truth[2].parse().unwrap());
        let position: u64 = column(&tsv, row, "position").parse().unwrap();
        assert!(
            (right - 5..=left + 5).contains(&position),
            "position {position} not within 5 bp of {right}-{left}"
        );
        assert_eq!(column(&tsv, row, "filter"), "PASS");
        assert_eq!(column(&tsv, row, "confidence"), "high");
        assert_eq!(column(&tsv, row, "contig"), truth[0]);
        assert_eq!(column(&tsv, row, "left_junction"), truth[1]);
        assert_eq!(column(&tsv, row, "right_junction"), (right - 1).to_string());
        assert_eq!(column(&tsv, row, "tsd_length"), truth[5].len().to_string());
        assert_eq!(column(&tsv, row, "strand"), truth[10]);
        assert_eq!(column(&tsv, row, "both_sides"), "yes");
        // Typed: family, subfamily and orientation of the element.
        assert_eq!(column(&tsv, row, "family"), truth[7]);
        assert_eq!(column(&tsv, row, "subfamily"), truth[6]);
        assert_eq!(column(&tsv, row, "element_strand"), truth[10]);
        assert!(column(&tsv, row, "typed_clips").parse::<u32>().unwrap() >= 2);
    }
}

#[test]
fn assembly_resolves_the_tsd_and_the_tail() {
    let tsv = scan(&[&POSITIVE[..], &["--reference", "tests/data/reference.fa"]].concat());
    let assembled = rows(&tsv);
    let truth = fs::read_to_string(Path::new(ROOT).join("tests/data/truth.tsv")).unwrap();
    assert_eq!(assembled.len(), 3, "{tsv}");
    for (row, truth) in assembled.iter().zip(truth.lines().skip(1)) {
        let truth: Vec<&str> = truth.split('\t').collect();
        assert_eq!(column(&tsv, row, "assembled"), "both");
        let (mut tsd, mut polya) = (truth[5].to_owned(), truth[9].to_owned());
        let mut left: u64 = truth[1].parse().unwrap();
        if truth[6] == "L1HS_3end" {
            // The poly(T) of the minus-strand LINE-1 starts with `TT`, the
            // reference after its TSD too (`TTC`): the contig matches the
            // reference over them, and the same alternative haplotype reads
            // as a TSD 2 bp longer and a tail 2 bp shorter.
            tsd.push_str("TT");
            polya = "28".to_owned();
            left += 2;
        }
        assert_eq!(column(&tsv, row, "tsd"), tsd);
        assert_eq!(column(&tsv, row, "tsd_length"), tsd.len().to_string());
        assert_eq!(column(&tsv, row, "polya_length"), polya);
        assert_eq!(column(&tsv, row, "left_junction"), left.to_string());
        assert_eq!(column(&tsv, row, "confidence"), "high");
    }
    // Without the reference, nothing is assembled.
    let tsv = scan(&POSITIVE);
    for row in &rows(&tsv) {
        assert_eq!(column(&tsv, row, "assembled"), ".");
    }
}

#[test]
fn negative_sample_has_no_pass_candidate() {
    let tsv = scan(&[
        "--bam",
        "tests/data/negative.bam",
        "--bed",
        "tests/data/targets.bed",
    ]);
    assert_eq!(rows(&tsv), Vec::<Vec<String>>::new());

    // With --keep-filtered, only the probe-edge artefact shows.
    let tsv = scan(&[
        "--bam",
        "tests/data/negative.bam",
        "--bed",
        "tests/data/targets.bed",
        "--keep-filtered",
    ]);
    let rows = rows(&tsv);
    assert_eq!(rows.len(), 1);
    // No poly(A/T) tail in the artefact either.
    assert_eq!(column(&tsv, &rows[0], "filter"), "untyped;probe_edge");
    assert_eq!(column(&tsv, &rows[0], "position"), "3000");
}

#[test]
fn insertion_at_probe_edges_stays_pass() {
    // The insertion target tiled with two probes whose edges fall on the two
    // junctions: every soft-clip stops at a probe edge, but the poly(A) tail
    // and the TSD rule out a capture artefact. The artefact stays flagged.
    let tsv = scan(&[
        "--bam",
        "tests/data/positive.bam",
        "--bed",
        "tests/data/targets_tiled.bed",
        "--keep-filtered",
    ]);
    let rows = rows(&tsv);
    // The artefact, the three insertions and the chimeric junction.
    assert_eq!(rows.len(), 5, "{tsv}");
    let artefact = rows.iter().find(|r| column(&tsv, r, "position") == "3000");
    assert_eq!(
        column(&tsv, artefact.unwrap(), "filter"),
        "untyped;probe_edge"
    );
    let insertion = rows.iter().find(|r| column(&tsv, r, "position") == "10012");
    let insertion = insertion.unwrap();
    assert_eq!(column(&tsv, insertion, "filter"), "PASS");
    assert_eq!(column(&tsv, insertion, "confidence"), "high");

    // The clips are exactly at the edges: the hallmarks keep it, not the
    // tolerance.
    let exact = scan(&[
        "--bam",
        "tests/data/positive.bam",
        "--bed",
        "tests/data/targets_tiled.bed",
        "--probe-edge-tolerance",
        "0",
    ]);
    assert_eq!(
        exact.lines().filter(|l| !l.starts_with('#')).count(),
        3,
        "{exact}"
    );
}

#[test]
fn random_mates_are_not_mei_mates() {
    // The chimeric junction of `positive` (tests/data/README.md): typed Alu
    // clips, but random unmapped mates, none of them a mobile element.
    let tsv = scan(&[&POSITIVE[..], &["--keep-filtered"]].concat());
    let found = rows(&tsv);
    let chimera = found
        .iter()
        .find(|r| column(&tsv, r, "contig") == "chr2")
        .unwrap();
    assert_eq!(column(&tsv, chimera, "filter"), "non_mei_mates");
    assert_eq!(column(&tsv, chimera, "family"), "ALU");
    assert_eq!(column(&tsv, chimera, "checked_mates"), "6");
    assert_eq!(column(&tsv, chimera, "mei_mates"), "0");

    // Most mates of the Alu map to its chr2 copy with MAPQ 0: without the
    // reference only the unmapped ones among the 40 sampled are checked;
    // with it, all 40 are, and all are MEI mates.
    for (reference, checked, mei) in [(false, 4, 4), (true, 40, 40)] {
        let mut args = POSITIVE.to_vec();
        if reference {
            args.extend(["--reference", "tests/data/reference.fa"]);
        }
        let tsv = scan(&args);
        let found = rows(&tsv);
        let alu = found
            .iter()
            .find(|r| column(&tsv, r, "family") == "ALU")
            .unwrap();
        assert_eq!(column(&tsv, alu, "checked_mates"), checked.to_string());
        assert_eq!(column(&tsv, alu, "mei_mates"), mei.to_string());
    }
}

#[test]
fn reference_tract_length_variant() {
    // A heterozygous expansion of a reference poly(A) tract
    // (tests/data/README.md): without the reference, its clips look like an
    // insertion with a poly(A) tail and a 12 bp TSD.
    let tract = ["--bed", "tests/data/targets_tract.bed", "--keep-filtered"];
    let tsv = scan(&[&["--bam", "tests/data/positive.bam"][..], &tract].concat());
    let found = rows(&tsv);
    assert_eq!(found.len(), 1, "{tsv}");
    assert_eq!(column(&tsv, &found[0], "filter"), "PASS");
    assert_eq!(column(&tsv, &found[0], "confidence"), "high");
    assert_eq!(column(&tsv, &found[0], "tsd_length"), "12");

    // With the reference, the clips are the reference with a longer tract.
    for bam in ["tests/data/positive.bam", "tests/data/positive.cram"] {
        let tsv = scan(
            &[
                &["--bam", bam, "--reference", "tests/data/reference.fa"][..],
                &tract,
            ]
            .concat(),
        );
        let found = rows(&tsv);
        assert_eq!(found.len(), 1, "{tsv}");
        assert_eq!(column(&tsv, &found[0], "filter"), "tandem_repeat", "{bam}");
    }

    // The simulated insertions stay PASS.
    let tsv = scan(&[&POSITIVE[..], &["--reference", "tests/data/reference.fa"]].concat());
    let found = rows(&tsv);
    assert_eq!(found.len(), 3, "{tsv}");
    for row in &found {
        assert_eq!(column(&tsv, row, "confidence"), "high");
    }
}

#[test]
fn bam_and_cram_give_the_same_output() {
    for format in ["tsv", "bed", "json"] {
        let bam = scan(
            &[
                &POSITIVE[..],
                &["--reference", "tests/data/reference.fa", "--format", format],
            ]
            .concat(),
        );
        let cram = scan(&[
            "--bam",
            "tests/data/positive.cram",
            "--bed",
            "tests/data/targets.bed",
            "--reference",
            "tests/data/reference.fa",
            "--format",
            format,
        ]);
        assert_eq!(normalize(&bam), normalize(&cram), "{format}");
    }
}

#[test]
fn same_output_with_any_thread_count() {
    let one = scan(&[&POSITIVE[..], &["--format", "json", "--threads", "1"]].concat());
    let eight = scan(&[&POSITIVE[..], &["--format", "json", "--threads", "8"]].concat());
    assert_eq!(normalize(&one), normalize(&eight));
}

#[test]
fn expected_outputs() {
    let update = std::env::var_os("UPDATE_EXPECTED").is_some();
    for (format, extension) in [("tsv", "tsv"), ("bed", "bed"), ("json", "json")] {
        let output = normalize(&scan(
            &[&POSITIVE[..], &["--keep-filtered", "--format", format]].concat(),
        )) + "\n";
        let path = Path::new(ROOT)
            .join("tests/data/expected")
            .join(format!("positive.{extension}"));
        if update {
            fs::write(&path, &output).unwrap();
        }
        let expected = fs::read_to_string(&path).unwrap();
        assert_eq!(
            output,
            expected,
            "{} differs (UPDATE_EXPECTED=1 to update)",
            path.display()
        );
    }
}

#[test]
fn output_file() {
    let dir = std::env::temp_dir().join(format!("mei-rs-output-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let path = dir.join("candidates.bed");
    let stdout = scan(
        &[
            &POSITIVE[..],
            &["--format", "bed", "-o", path.to_str().unwrap()],
        ]
        .concat(),
    );
    assert_eq!(stdout, "");
    let bed = fs::read_to_string(&path).unwrap();
    assert_eq!(bed.lines().count(), 3);
    assert!(bed.starts_with("chr1\t6000\t6010\tMEI_1\t"), "{bed}");
    fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn input_errors() {
    let dir = std::env::temp_dir().join(format!("mei-rs-errors-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let bed = dir.join("nochr.bed");
    fs::write(&bed, "1\t9850\t10150\n").unwrap();
    for (args, message) in [
        (
            vec![
                "--bam",
                "tests/data/positive.bam",
                "--bed",
                bed.to_str().unwrap(),
            ],
            "named differently",
        ),
        (
            vec![
                "--bam",
                "tests/data/positive.cram",
                "--bed",
                "tests/data/targets.bed",
            ],
            "--reference is required",
        ),
    ] {
        Command::cargo_bin("mei-rs")
            .unwrap()
            .current_dir(ROOT)
            .arg("scan")
            .args(args)
            .assert()
            .code(1)
            .stdout("")
            .stderr(predicate::str::contains(message));
    }
    fs::remove_dir_all(&dir).unwrap();
}
