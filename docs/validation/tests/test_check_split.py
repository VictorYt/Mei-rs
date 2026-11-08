"""Tests of `check_split.py`."""

import pathlib

import check_split

SPLIT = pathlib.Path(__file__).resolve().parents[1] / "split.tsv"


def row(sample, set_, family):
    """Make one row of `split.tsv`."""
    return {"sample": sample, "set": set_, "family_id": family}


def test_published_split_has_no_leak():
    """The split of the repository passes."""
    assert check_split.split_leaks(check_split.read_split(SPLIT)) == []


def test_family_across_tuning_and_test():
    """A trio with the mother in the test set fails."""
    rows = [row("child", "tuning", "F1"), row("father", "tuning", "F1"), row("mother", "test", "F1")]
    problems = check_split.split_leaks(rows)
    assert len(problems) == 1 and "F1" in problems[0] and "mother (test)" in problems[0]


def test_related_sample_needs_kin_in_tuning():
    """A `related` sample is left out because of a relative in the tuning set."""
    assert check_split.split_leaks([row("a", "tuning", "F1"), row("b", "related", "F1")]) == []
    assert len(check_split.split_leaks([row("b", "related", "F2")])) == 1


def test_unknown_set():
    """A misspelt set is reported."""
    assert "unknown set" in check_split.split_leaks([row("a", "tunning", "F1")])[0]


def test_pedigree_mismatch():
    """A family that differs from the pedigree is reported."""
    rows = [row("a", "tuning", "F1"), row("b", "test", "b")]
    assert check_split.pedigree_mismatches(rows, {"a": "F1", "b": "F2"}) == ["b: family_id b, pedigree F2"]
