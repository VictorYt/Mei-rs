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
fn positive_sample_has_one_pass_candidate_at_the_insertion() {
    let tsv = scan(&POSITIVE);
    let rows = rows(&tsv);
    assert_eq!(rows.len(), 1, "{tsv}");
    let row = &rows[0];
    // Truth: tests/data/truth.tsv (insertion between the TSD copies of
    // chr1:10001-10012).
    let truth = fs::read_to_string(Path::new(ROOT).join("tests/data/truth.tsv")).unwrap();
    let truth: Vec<&str> = truth.lines().nth(1).unwrap().split('\t').collect();
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
    assert_eq!(column(&tsv, row, "strand"), truth[9]);
    assert_eq!(column(&tsv, row, "both_sides"), "yes");
    assert!(column(&tsv, row, "polya_clips").parse::<u32>().unwrap() > 0);
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
    assert_eq!(column(&tsv, &rows[0], "filter"), "no_tail;probe_edge");
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
    assert_eq!(rows.len(), 2, "{tsv}");
    let artefact = rows.iter().find(|r| column(&tsv, r, "position") == "3000");
    assert_eq!(
        column(&tsv, artefact.unwrap(), "filter"),
        "no_tail;probe_edge"
    );
    let insertion = rows.iter().find(|r| column(&tsv, r, "position") != "3000");
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
        1,
        "{exact}"
    );
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
    assert_eq!(bed.lines().count(), 1);
    assert!(bed.starts_with("chr1\t10000\t10012\tMEI_1\t"), "{bed}");
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
