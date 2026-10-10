#!/usr/bin/env python3
"""Precision and recall of `mei-rs scan` candidates against the HGSVC3 MEI calls.

Truth set: the HGSVC3 mobile element insertions called against GRCh38 from
long-read assemblies (`MEI_Callset_GRCh38.ALL.20241211.csv.gz`, see README.md).
A truth MEI belongs to a sample when one of its haplotypes carries it.

Coordinates:
  * truth `POS` is VCF 1-based: the insertion is right after base `POS`;
  * `mei-rs scan` TSV positions are 1-based with the same meaning;
  * the BED targets are 0-based half-open: a truth MEI is kept when its
    insertion point (boundary `POS`, between bases `POS` and `POS + 1`) lies in
    a target extended by `--padding` bases on each side.

A candidate matches a truth MEI of the same sample when `position`,
`left_junction` or `right_junction` lies within `--tolerance` bases of `POS`
(the junctions of an insertion with a TSD are a TSD length apart).

Precision counts a candidate as a false positive only when no truth MEI of the
sample is near it and the sample is genotyped (no missing haplotype) at every
truth MEI within the tolerance; candidates near an MEI carried by another
HGSVC sample are reported apart (`other_sample`), since they point to genotype
errors in either call set or to real population variants.

Only the Python standard library is used.
"""

import argparse
import bisect
import collections
import csv
import gzip
import sys

FAMILIES = {
    "SINE/Alu": "Alu",
    "LINE/L1": "L1",
    "Retroposon/SVA": "SVA",
}


def open_text(path):
    """Open a plain or gzipped text file."""
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_targets(path, padding):
    """Merged padded targets per contig, as sorted 0-based half-open intervals."""
    raw = collections.defaultdict(list)
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            contig, start, end = line.split("\t")[:3]
            raw[contig].append((max(0, int(start) - padding), int(end) + padding))
    merged = {}
    for contig, intervals in raw.items():
        intervals.sort()
        out = [list(intervals[0])]
        for start, end in intervals[1:]:
            if start <= out[-1][1]:
                out[-1][1] = max(out[-1][1], end)
            else:
                out.append([start, end])
        merged[contig] = ([s for s, _ in out], [e for _, e in out])
    return merged


def in_targets(targets, contig, boundary):
    """Whether the boundary `boundary` (between bases) lies in a padded target."""
    if contig not in targets:
        return False
    starts, ends = targets[contig]
    i = bisect.bisect_right(starts, boundary) - 1
    return i >= 0 and boundary <= ends[i]


def target_distance(targets, contig, boundary):
    """Distance from the boundary `boundary` to the nearest unpadded target (0 inside)."""
    starts, ends = targets.get(contig, ([], []))
    i = bisect.bisect_right(starts, boundary) - 1
    best = None
    for j in (i, i + 1):
        if 0 <= j < len(starts):
            d = max(0, starts[j] - boundary, boundary - ends[j])
            best = d if best is None else min(best, d)
    return best


def load_truth(path, samples, targets, raw_targets):
    """Truth MEIs in the targets: list of dicts with the genotype of every sample."""
    truth = []
    with open_text(path) as handle:
        reader = csv.reader(handle)
        header = next(reader)
        col = {name: i for i, name in enumerate(header)}
        missing = [s for s in samples if s not in col]
        if missing:
            sys.exit(f"samples not in the truth set: {', '.join(missing)}")
        others = [h for h in header[col["FORMAT"] + 1 : col["Caller_Count"]]]
        for row in reader:
            contig, pos = row[col["CHROM"]], int(row[col["POS"]])
            if not in_targets(targets, contig, pos):
                continue
            genotypes = {s: row[col[s]] for s in others}
            truth.append(
                {
                    "id": row[col["ID"]],
                    "contig": contig,
                    "pos": pos,
                    "family": FAMILIES.get(row[col["TE_Designation"]], "other"),
                    "distance": target_distance(raw_targets, contig, pos),
                    "genotypes": genotypes,
                }
            )
    return truth


def carries(genotype):
    """A haplotype carries the insertion."""
    return "1" in genotype.split("|")


def missing(genotype):
    """At least one haplotype is not genotyped."""
    return "." in genotype.split("|")[:2]


def load_candidates(path):
    """Candidates of a `mei-rs scan` TSV (written with `--keep-filtered`)."""
    candidates = []
    with open(path) as handle:
        columns = None
        for line in handle:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                columns = [c.lstrip("#") for c in fields]
                continue
            row = dict(zip(columns, fields))
            row["points"] = [
                int(row[k])
                for k in ("position", "left_junction", "right_junction")
                if row[k] != "."
            ]
            candidates.append(row)
    return candidates


def tier(candidate):
    """`high`, `medium`, `low` for PASS candidates, else the first filter."""
    if candidate["filter"] == "PASS":
        return candidate["confidence"]
    return candidate["filter"].split(";")[0]


def evaluate(sample, candidates, truth, tolerance):
    """Match the candidates of one sample against the truth MEIs."""
    by_contig = collections.defaultdict(list)
    for t in truth:
        by_contig[t["contig"]].append(t)
    for items in by_contig.values():
        items.sort(key=lambda t: t["pos"])
    positions = {c: [t["pos"] for t in items] for c, items in by_contig.items()}

    def near(candidate):
        items = by_contig.get(candidate["contig"], [])
        pos = positions.get(candidate["contig"], [])
        lo = bisect.bisect_left(pos, min(candidate["points"]) - tolerance)
        hi = bisect.bisect_right(pos, max(candidate["points"]) + tolerance)
        return [
            t
            for t in items[lo:hi]
            if any(abs(p - t["pos"]) <= tolerance for p in candidate["points"])
        ]

    found = collections.defaultdict(set)  # truth id -> tiers of matching candidates
    rows = []
    for c in candidates:
        hits = near(c)
        own = [t for t in hits if carries(t["genotypes"][sample])]
        if own:
            status = "tp"
            for t in own:
                found[t["id"]].add(tier(c))
        elif any(missing(t["genotypes"][sample]) for t in hits):
            status = "unknown"
        elif any(any(carries(g) for g in t["genotypes"].values()) for t in hits):
            status = "other_sample"
        else:
            status = "fp"
        rows.append((c, status))
    carried = [t for t in truth if carries(t["genotypes"][sample])]
    return rows, carried, found


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--truth", required=True, help="MEI_Callset_GRCh38 CSV(.gz)")
    parser.add_argument("--bed", required=True, help="capture targets (BED)")
    parser.add_argument("--padding", type=int, default=300, help="as in `scan` [300]")
    parser.add_argument("--tolerance", type=int, default=20, help="breakpoint distance [20]")
    parser.add_argument(
        "--fn-out", help="write the missed truth MEIs (PASS level) to this TSV"
    )
    parser.add_argument(
        "--fp-out", help="write the PASS false positives to this TSV"
    )
    parser.add_argument(
        "runs", nargs="+", metavar="SAMPLE=TSV", help="scan output of each sample"
    )
    args = parser.parse_args()

    runs = [r.split("=", 1) for r in args.runs]
    targets = load_targets(args.bed, args.padding)
    raw_targets = load_targets(args.bed, 0)
    truth = load_truth(args.truth, [s for s, _ in runs], targets, raw_targets)
    strata = [(0, 0, "in target"), (1, 100, "1-100 bp"), (101, args.padding, "padding >100 bp")]
    strata_total = collections.Counter()
    strata_any = collections.Counter()
    strata_pass = collections.Counter()

    tiers = ["high", "medium", "low"]
    total = collections.Counter()
    family_total = collections.Counter()
    family_found = collections.Counter()
    fn_rows, fp_rows = [], []
    print(
        f"# truth MEIs in the padded targets (any HGSVC3 sample): {len(truth)}; "
        f"tolerance {args.tolerance} bp"
    )
    print(
        "sample\ttruth\tfound_any\tfound_pass\tpass\ttp\tfp\tother_sample\tunknown"
        "\tprecision\trecall_pass\trecall_any\t"
        + "\t".join(f"tp_{t}/{t}" for t in tiers)
    )
    for sample, path in runs:
        candidates = load_candidates(path)
        rows, carried, found = evaluate(sample, candidates, truth, args.tolerance)
        passing = [(c, s) for c, s in rows if c["filter"] == "PASS"]
        status = collections.Counter(s for _, s in passing)
        found_pass = {
            i for i, ts in found.items() if ts & set(tiers)
        }
        per_tier = []
        for t in tiers:
            in_tier = [s for c, s in passing if c["confidence"] == t]
            per_tier.append(f"{in_tier.count('tp')}/{len(in_tier)}")
            total[f"tp_{t}"] += in_tier.count("tp")
            total[f"n_{t}"] += len(in_tier)
        for t in carried:
            for lo, hi, name in strata:
                if lo <= t["distance"] <= hi:
                    strata_total[name] += 1
                    strata_any[name] += t["id"] in found
                    strata_pass[name] += t["id"] in found_pass
            family_total[t["family"]] += 1
            family_found[t["family"]] += t["id"] in found_pass
            if t["id"] not in found_pass and args.fn_out:
                fn_rows.append(
                    (sample, t["id"], t["family"], ",".join(sorted(found.get(t["id"], []))) or ".")
                )
        for c, s in passing:
            if s == "fp" and args.fp_out:
                fp_rows.append((sample, c))
        values = {
            "truth": len(carried),
            "found_any": len(found),
            "found_pass": len(found_pass),
            "pass": len(passing),
            "tp": status["tp"],
            "fp": status["fp"],
            "other_sample": status["other_sample"],
            "unknown": status["unknown"],
        }
        total.update(values)
        print(format_row(sample, values, per_tier))
    print(
        format_row(
            "all",
            total,
            [f"{total[f'tp_{t}']}/{total[f'n_{t}']}" for t in tiers],
        )
    )
    print("\n# recall by distance to the target\tfound_any\tfound_pass")
    for _, _, name in strata:
        if strata_total[name]:
            print(
                f"{name}\t{strata_any[name]}/{strata_total[name]}"
                f"\t{strata_pass[name]}/{strata_total[name]}"
            )
    print("\n# recall (PASS) by family")
    for family in ("Alu", "L1", "SVA", "other"):
        if family_total[family]:
            print(
                f"{family}\t{family_found[family]}/{family_total[family]}"
                f"\t{family_found[family] / family_total[family]:.3f}"
            )
    if args.fn_out:
        with open(args.fn_out, "w") as out:
            out.write("sample\ttruth_id\tfamily\tcandidate_tiers\n")
            for row in fn_rows:
                out.write("\t".join(row) + "\n")
    if args.fp_out:
        with open(args.fp_out, "w") as out:
            keys = ["contig", "position", "confidence", "support", "tsd_length", "strand",
                    "left_junction_clips", "right_junction_clips", "polya_clips",
                    "polyt_clips", "unmapped_mates", "discordant_mates", "depth"]
            out.write("sample\t" + "\t".join(keys) + "\n")
            for sample, c in fp_rows:
                out.write(sample + "\t" + "\t".join(c[k] for k in keys) + "\n")


def format_row(name, v, per_tier):
    """One line of the results table."""
    decided = v["tp"] + v["fp"]
    precision = f"{v['tp'] / decided:.3f}" if decided else "."
    recall_pass = f"{v['found_pass'] / v['truth']:.3f}" if v["truth"] else "."
    recall_any = f"{v['found_any'] / v['truth']:.3f}" if v["truth"] else "."
    return "\t".join(
        [name]
        + [str(v[k]) for k in ("truth", "found_any", "found_pass", "pass", "tp", "fp",
                               "other_sample", "unknown")]
        + [precision, recall_pass, recall_any]
        + per_tier
    )


if __name__ == "__main__":
    main()
