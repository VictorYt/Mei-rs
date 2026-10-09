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

#[test]
fn scan_is_not_implemented_yet() {
    mei_rs()
        .args(["scan", "--bam", "s.bam", "--bed", "t.bed"])
        .assert()
        .code(1)
        .stderr(predicate::str::contains("not implemented yet"));
}
