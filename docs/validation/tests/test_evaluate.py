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


def per_sample(samples, truths, calls_of):
    """Build the input of `Loci` and `outcomes`: PASS rows of every sample."""
    index = evaluate.TruthIndex(truths)
    out = []
    for s in samples:
        carried = [t for t in truths if evaluate.carries(t["genotypes"][s])]
        out.append((s, carried, {"PASS": evaluate.evaluate(s, calls_of[s], index, 20)}))
    return out


def test_loci_join_samples_at_one_site():
    """An MEI and the false calls of other samples nearby are one locus; far calls are another."""
    genotypes = {"A": "1|0", "B": "0|0"}
    data = per_sample(["A", "B"], [truth("a", 1000, genotypes=genotypes)], {"A": [call(1000)], "B": [call(1050)]})
    loci = evaluate.Loci(data)
    assert loci.count == 1
    cells = evaluate.outcomes(data, loci, "PASS")
    assert cells == {("A", 0): [1, 0, 0], ("B", 0): [0, 1, 0]}
    data = per_sample(["A", "B"], [truth("a", 1000, genotypes=genotypes)], {"A": [], "B": [call(5000)]})
    loci = evaluate.Loci(data)
    assert loci.count == 2
    assert evaluate.outcomes(data, loci, "PASS") == {("A", 0): [0, 0, 1], ("B", 1): [0, 1, 0]}


def test_two_way_interval_wider_for_a_recurrent_locus():
    """A false call repeated at one locus in every sample widens the interval once loci are resampled."""
    cells = []
    for i in range(10):
        s = f"S{i}"
        cells.append((s, s, 0, [0, 1, 0]))  # the same artefact in every sample
        cells += [(s, s, 1 + i * 3 + j, [1, 0, 0]) for j in range(3)]
    groups = [f"S{i}" for i in range(10)]
    samples_only = evaluate.bootstrap(cells, groups, 31, 500)
    two_way = evaluate.bootstrap(cells, groups, 31, 500, by_locus=True)
    assert samples_only["precision"] == ("0.750", "0.750")
    low, high = (float(x) for x in two_way["precision"])
    assert low < 0.75 < high


def test_families_resampled_together():
    """Two samples of one family move together: with one family, the interval has no width."""
    cells = [("F", "A", 0, [1, 0, 0]), ("F", "B", 1, [0, 1, 0])]
    assert evaluate.bootstrap(cells, ["F"], 2, 200)["precision"] == ("0.500", "0.500")


def test_clopper_pearson_lower_bound():
    """59 of 59 found is the classic 95% claim; known values of the exact bound."""
    assert evaluate.clopper_pearson_lower(59, 59) > 0.95
    assert evaluate.clopper_pearson_lower(58, 58) < 0.95
    assert abs(evaluate.clopper_pearson_lower(5, 10) - 0.2224) < 1e-3
    assert evaluate.clopper_pearson_lower(0, 10) == 0.0


def test_bins():
    """Half-open bins named by their bounds."""
    assert [evaluate.bins(v, [100, 150]) for v in (99, 100, 149, 150)] == ["< 100", "100-149", "100-149", ">= 150"]


def test_haploid_calls_agree_with_any_carrier():
    """A haploid `1` on the X of a male agrees with the truth whatever its zygosity."""
    assert evaluate.concordant("0/1", "0/1")
    assert not evaluate.concordant("1/1", "0/1")
    assert evaluate.concordant("0/1", "1") and evaluate.concordant("1/1", "1")
    assert not evaluate.concordant("0/1", "0")
