//! End-to-end tests of the command-line interface.

use assert_cmd::Command;
use predicates::prelude::*;

fn mei_rs() -> Command {
    Command::cargo_bin("mei-rs").unwrap()
}

#[test]
fn help_lists_scan() {
    mei_rs()
        .arg("--help")
        .assert()
        .success()
        .stdout(predicate::str::contains("scan"));
}

#[test]
fn scan_help_lists_every_option() {
    let assert = mei_rs().args(["scan", "--help"]).assert().success();
    let stdout = String::from_utf8(assert.get_output().stdout.clone()).unwrap();
    for option in [
        "--bam",
        "--bed",
        "--reference",
        "--padding",
        "--region",
        "--output",
        "--format",
        "--threads",
        "--verbose",
        "--quiet",
    ] {
        assert!(stdout.contains(option), "{option} missing from scan --help");
    }
}

#[test]
fn version() {
    mei_rs()
        .arg("--version")
        .assert()
        .success()
        .stdout(predicate::str::starts_with(concat!(
            "mei-rs ",
            env!("CARGO_PKG_VERSION")
        )));
}

#[test]
fn invalid_option_fails_with_usage_error() {
    mei_rs()
        .args([
            "scan", "--bam", "s.bam", "--bed", "t.bed", "--format", "vcf",
        ])
        .assert()
        .code(2)
        .stderr(predicate::str::contains("invalid value 'vcf'"));
}

fn data(name: &str) -> String {
    format!("{}/tests/data/{name}", env!("CARGO_MANIFEST_DIR"))
}

#[test]
fn scan_checks_the_inputs_then_stops() {
    for (bam, reference) in [
        ("positive.bam", None),
        ("positive.cram", Some("reference.fa")),
    ] {
        let mut cmd = mei_rs();
        cmd.args(["scan", "--bam", &data(bam), "--bed", &data("targets.bed")]);
        if let Some(reference) = reference {
            cmd.args(["--reference", &data(reference)]);
        }
        cmd.assert()
            .code(1)
            .stderr(predicate::str::contains("sample positive"))
            .stderr(predicate::str::contains(
                "3 targets (800 bp), 3 regions to scan (2600 bp with 300 bp of padding)",
            ))
            .stderr(predicate::str::contains(
                "signal extraction is not implemented yet",
            ));
    }
}

#[test]
fn scan_input_errors() {
    let bed = data("targets.bed");
    for (args, message) in [
        (
            vec!["--bam", "missing.bam", "--bed", &bed],
            "cannot read missing.bam",
        ),
        (
            vec!["--bam", &data("positive.bam"), "--bed", "missing.bed"],
            "cannot read BED file missing.bed",
        ),
        (
            vec!["--bam", &data("positive.cram"), "--bed", &bed],
            "--reference is required",
        ),
    ] {
        mei_rs()
            .arg("scan")
            .args(args)
            .assert()
            .code(1)
            .stderr(predicate::str::contains(message));
    }
}

#[test]
fn scan_region() {
    mei_rs()
        .args([
            "scan",
            "--bam",
            &data("positive.bam"),
            "--bed",
            &data("targets.bed"),
        ])
        .args(["--region", "chr1:1-100"])
        .assert()
        .code(1)
        .stderr(predicate::str::contains(
            "--region chr1:1-100 holds no target",
        ));
}
