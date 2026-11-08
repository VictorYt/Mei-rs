"""Tests of the matching of `evaluate.py`."""

import evaluate


def truth(id_, pos, family="Alu", genotypes=None):
    """Make a truth MEI on chr1; sample S carries it unless `genotypes` says otherwise."""
    return {
        "id": id_,
        "contig": "1",
        "pos": pos,
        "family": family,
        "distance": 0,
        "genotypes": genotypes or {"S": "1|0", "O": "0|0"},
    }


def call(*points, family="Alu"):
    """Make a PASS call on chr1 at `points` (position, junctions)."""
    return {"contig": "1", "points": list(points), "family": family, "pass": True, "imprecise": False}


def statuses(calls, truths, tolerance=20):
    """Match the calls of sample S and return their statuses."""
    return [s for _, s, _, _ in evaluate.evaluate("S", calls, evaluate.TruthIndex(truths), tolerance)]


def test_one_call_one_truth():
    """A call within the tolerance of an MEI of the sample is true."""
    assert statuses([call(1000)], [truth("a", 1010)]) == ["tp"]
    assert statuses([call(1000)], [truth("a", 1021)]) == ["fp"]


def test_junction_within_tolerance():
    """The left or right junction may be the point that matches."""
    assert statuses([call(1000, 1000, 1015)], [truth("a", 1030)]) == ["tp"]


def test_second_call_is_duplicate():
    """Two calls of one MEI: the closest is true, the other a duplicate."""
    rows = evaluate.evaluate("S", [call(1012), call(1003)], evaluate.TruthIndex([truth("a", 1000)]), 20)
    assert [s for _, s, _, _ in rows] == ["duplicate", "tp"]
    assert evaluate.found(rows) == {"a"}
    assert evaluate.count(rows)["false"] == 1


def test_one_call_credits_one_truth():
    """A call between two MEIs credits only the closest; the second needs its own call."""
    truths = [truth("a", 1000), truth("b", 1015)]
    rows = evaluate.evaluate("S", [call(1012)], evaluate.TruthIndex(truths), 20)
    assert evaluate.found(rows) == {"b"}
    rows = evaluate.evaluate("S", [call(1012), call(995)], evaluate.TruthIndex(truths), 20)
    assert evaluate.found(rows) == {"a", "b"}


def test_closest_pairs_first():
    """Greedy matching takes the closest pair, then matches the other call to the other MEI."""
    truths = [truth("a", 1000), truth("b", 1030)]
    rows = evaluate.evaluate("S", [call(1018), call(1012)], evaluate.TruthIndex(truths), 20)
    assert [(s, t["id"]) for _, s, t, _ in rows] == [("tp", "b"), ("tp", "a")]


def test_wrong_family():
    """An Alu call on a LINE-1 is false; an untyped call or an `other` MEI is not checked."""
    assert statuses([call(1000)], [truth("a", 1000, "L1")]) == ["wrong_family"]
    assert statuses([call(1000, family=None)], [truth("a", 1000, "L1")]) == ["tp"]
    assert statuses([call(1000, family="SVA")], [truth("a", 1000, "other")]) == ["tp"]


def test_wrong_family_does_not_find():
    """A call of the wrong family does not count in the recall."""
    rows = evaluate.evaluate("S", [call(1000)], evaluate.TruthIndex([truth("a", 1000, "L1")]), 20)
    assert evaluate.found(rows) == set()


def test_other_sample_and_unknown():
    """A call near an MEI of another sample, or of a haplotype not genotyped, is neither true nor false."""
    assert statuses([call(1000)], [truth("a", 1000, genotypes={"S": "0|0", "O": "1|0"})]) == ["other_sample"]
    assert statuses([call(1000)], [truth("a", 1000, genotypes={"S": ".|0", "O": "1|0"})]) == ["unknown"]
    assert statuses([call(1000)], [truth("a", 1000, genotypes={"S": "0|0", "O": "0|0"})]) == ["fp"]


def test_tolerance():
    """A larger tolerance matches farther calls."""
    assert statuses([call(1000)], [truth("a", 1050)], tolerance=50) == ["tp"]
    assert statuses([call(1000)], [truth("a", 1050)], tolerance=20) == ["fp"]
