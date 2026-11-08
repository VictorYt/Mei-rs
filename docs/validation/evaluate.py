#!/usr/bin/env python3
"""Precision and recall of MEI calls against the HGSVC3 MEI calls.

Truth set: the HGSVC3 mobile element insertions called against GRCh38 from
long-read assemblies (`MEI_Callset_GRCh38.ALL.20241211.csv.gz`, see README.md).
A truth MEI belongs to a sample when one of its haplotypes carries it.

Coordinates:
  * truth `POS` is VCF 1-based: the insertion is right after base `POS`;
  * `mei-rs scan` positions (TSV and VCF) are 1-based with the same meaning;
  * the BED targets are 0-based half-open: a truth MEI is kept when its
    insertion point (boundary `POS`, between bases `POS` and `POS + 1`) lies in
    a target extended by `--padding` bases on each side.

Truth positions can be lifted to another assembly with `--chain` (a UCSC
chain file, e.g. `hg38ToHg19.over.chain.gz` for GRCh37 exomes). Contig names
are compared without a leading `chr` (`chr1` and `1` are the same contig).

Calls are read from a `mei-rs scan` TSV (written with `--keep-filtered`) or
from a VCF (`.vcf` or `.vcf.gz`) of `mei-rs` or of another caller: `POS`, and
the `LEFT_JUNCTION` and `RIGHT_JUNCTION` INFO fields when present; a record
passes when its `FILTER` is `PASS` or `.`.

Calls outside the padded targets are left out (callers that scan the whole
alignment file call insertions off target too).

A call matches a truth MEI of the same sample when its position, left junction
or right junction lies within `--tolerance` bases of `POS` (the junctions of
an insertion with a TSD are a TSD length apart).

Precision counts a call as a false positive only when no truth MEI of the
sample is near it and the sample is genotyped (no missing haplotype) at every
truth MEI within the tolerance; calls near an MEI carried by another HGSVC
sample are reported apart (`other_sample`), since they point to genotype
errors in either call set or to real population variants.

Only the Python standard library is used.
"""

import argparse
import bisect
import collections
import csv
import gzip
import random
import sys

FAMILIES = {
    "SINE/Alu": "Alu",
    "LINE/L1": "L1",
    "Retroposon/SVA": "SVA",
}
# Family of a call: `family` column of the TSV, last field of a symbolic ALT (`<INS:ME:ALU>`).
CALL_FAMILIES = {"ALU": "Alu", "LINE1": "L1", "SVA": "SVA"}
TIERS = ["high", "medium", "low"]


def open_text(path):
    """Open a plain or gzipped text file."""
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def contig_name(name):
    """Contig name without a leading `chr` (GRCh38 `chr1`, GRCh37 `1`)."""
    return name[3:] if name.startswith("chr") else name


class Chain:
    """Liftover of single bases with a UCSC chain file.

    Positions are 0-based on both assemblies; a base covered by several chains
    takes the one with the highest score.
    """

    def __init__(self, path):
        """Read the chain file `path`."""
        # Per source contig: (start, end, score, target contig, target size,
        # strand, block starts, [(block end, target start)]).
        self.chains = collections.defaultdict(list)
        with open_text(path) as handle:
            chain = None
            for line in handle:
                fields = line.split()
                if not fields:
                    continue
                if fields[0] == "chain":
                    t_pos, q_pos = int(fields[5]), int(fields[10])
                    chain = (
                        int(fields[5]),
                        int(fields[6]),
                        int(fields[1]),
                        contig_name(fields[7]),
                        int(fields[8]),
                        fields[9],
                        [],
                        [],
                    )
                    self.chains[contig_name(fields[2])].append(chain)
                    continue
                size = int(fields[0])
                chain[6].append(t_pos)
                chain[7].append((t_pos + size, q_pos))
                if len(fields) == 3:
                    t_pos += size + int(fields[1])
                    q_pos += size + int(fields[2])

    def lift(self, contig, base):
        """Contig, 0-based position and strand of the 0-based `base`, or None."""
        best = None
        for start, end, score, q_contig, q_size, strand, starts, blocks in self.chains.get(
            contig, []
        ):
            if not start <= base < end or (best is not None and score <= best[0]):
                continue
            i = bisect.bisect_right(starts, base) - 1
            if i < 0 or base >= blocks[i][0]:
                continue
            q = blocks[i][1] + base - starts[i]
            best = (score, q_contig, q if strand == "+" else q_size - 1 - q, strand)
        return None if best is None else best[1:]


def lift_insertion(chain, contig, pos):
    """Contig and 1-based `POS` of an insertion after base `pos` on the other assembly."""
    lifted = chain.lift(contig, pos - 1)
    if lifted is None:
        return None
    q_contig, q, strand = lifted
    # On the reverse strand, the insertion lies before the lifted base.
    return q_contig, q + 1 if strand == "+" else q


def load_targets(path, padding):
    """Read the targets: merged padded intervals per contig, sorted, 0-based half-open."""
    raw = collections.defaultdict(list)
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            contig, start, end = line.split("\t")[:3]
            raw[contig_name(contig)].append(
                (max(0, int(start) - padding), int(end) + padding)
            )
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


def load_truth(path, samples, targets, raw_targets, chain):
    """Truth MEIs in the targets: list of dicts with the genotype of every sample."""
    truth = []
    unlifted = 0
    with open_text(path) as handle:
        reader = csv.reader(handle)
        header = next(reader)
        col = {name: i for i, name in enumerate(header)}
        absent = [s for s in samples if s not in col]
        if absent:
            sys.exit(f"samples not in the truth set: {', '.join(absent)}")
        others = header[col["FORMAT"] + 1 : col["Caller_Count"]]
        for row in reader:
            contig, pos = contig_name(row[col["CHROM"]]), int(row[col["POS"]])
            if chain is not None:
                lifted = lift_insertion(chain, contig, pos)
                if lifted is None:
                    unlifted += 1
                    continue
                contig, pos = lifted
            if not in_targets(targets, contig, pos):
                continue
            truth.append(
                {
                    "id": row[col["ID"]],
                    "contig": contig,
                    "pos": pos,
                    "family": FAMILIES.get(row[col["TE_Designation"]], "other"),
                    "distance": target_distance(raw_targets, contig, pos),
                    "genotypes": {s: row[col[s]] for s in others},
                }
            )
    return truth, unlifted


def carries(genotype):
    """Tell whether a haplotype carries the insertion."""
    return "1" in genotype.split("|")


def missing(genotype):
    """At least one haplotype is not genotyped."""
    return "." in genotype.split("|")[:2]


def unphased(genotype):
    """Truth genotype as `0/1` or `1/1`."""
    return "1/1" if genotype.split("|")[:2] == ["1", "1"] else "0/1"


def load_tsv(path):
    """Read the calls of a `mei-rs scan` TSV (written with `--keep-filtered`)."""
    calls = []
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
            calls.append(
                {
                    "contig": contig_name(row["contig"]),
                    "points": [
                        int(row[k])
                        for k in ("position", "left_junction", "right_junction")
                        if row[k] != "."
                    ],
                    "pass": row["filter"] == "PASS",
                    "confidence": row["confidence"],
                    "qual": int(row["qual"]) if row.get("qual", ".") != "." else None,
                    "genotype": row.get("genotype", "."),
                    "family": CALL_FAMILIES.get(row.get("family", ".")),
                    "imprecise": row.get("cipos", ".") != ".",
                    "row": row,
                }
            )
    return calls


def load_vcf(path):
    """Read the calls of a VCF (first sample column for the genotype)."""
    calls = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info = dict(
                f.split("=", 1) if "=" in f else (f, "") for f in fields[7].split(";")
            )
            points = [int(fields[1])] + [
                int(info[k]) for k in ("LEFT_JUNCTION", "RIGHT_JUNCTION") if k in info
            ]
            genotype = "."
            if len(fields) > 9:
                keys, values = fields[8].split(":"), fields[9].split(":")
                genotype = dict(zip(keys, values)).get("GT", ".").replace("|", "/")
            calls.append(
                {
                    "contig": contig_name(fields[0]),
                    "points": points,
                    "pass": fields[6] in ("PASS", "."),
                    "confidence": info.get("CONFIDENCE", "."),
                    "qual": None if fields[5] == "." else float(fields[5]),
                    "genotype": genotype,
                    "family": CALL_FAMILIES.get(fields[4].strip("<>").split(":")[-1]),
                    "imprecise": "CIPOS" in info or "IMPRECISE" in info,
                    "row": {"contig": fields[0], "position": fields[1], "info": fields[7]},
                }
            )
    return calls


def load_calls(path):
    """Read the calls of a TSV or VCF file."""
    return load_vcf(path) if path.endswith((".vcf", ".vcf.gz")) else load_tsv(path)


class TruthIndex:
    """Truth MEIs sorted by position on each contig, for the searches within a tolerance."""

    def __init__(self, truth):
        """Index the truth MEIs `truth`."""
        self.by_contig = collections.defaultdict(list)
        for t in truth:
            self.by_contig[t["contig"]].append(t)
        for items in self.by_contig.values():
            items.sort(key=lambda t: t["pos"])
        self.positions = {c: [t["pos"] for t in items] for c, items in self.by_contig.items()}

    def near(self, call, tolerance):
        """List the truth MEIs within `tolerance` of the position or a junction of `call`."""
        items = self.by_contig.get(call["contig"], [])
        pos = self.positions.get(call["contig"], [])
        lo = bisect.bisect_left(pos, min(call["points"]) - tolerance)
        hi = bisect.bisect_right(pos, max(call["points"]) + tolerance)
        return [t for t in items[lo:hi] if distance(call, t) <= tolerance]


def distance(call, t):
    """Distance from the nearest point (position, junctions) of a call to a truth MEI."""
    return min(abs(p - t["pos"]) for p in call["points"])


def same_family(call, t):
    """Tell whether the family of a call agrees with a truth MEI (an untyped call or an `other` MEI always does)."""
    return call["family"] is None or t["family"] == "other" or call["family"] == t["family"]


def evaluate(sample, calls, index, tolerance):
    """Match the calls of one sample one to one with the truth MEIs it carries.

    Pairs of a call and a truth MEI of the sample within `tolerance` are taken
    closest first, each call and each MEI at most once. A matched pair is `tp`
    when the families agree, `wrong_family` otherwise; a call left with only
    matched MEIs of the sample nearby is a `duplicate`. Both are false calls.
    A call near no MEI of the sample is `unknown` when the sample is not
    genotyped at an MEI nearby, `other_sample` when another sample carries one,
    `fp` otherwise. Returns `(call, status, matched MEI or None, MEIs of the
    sample nearby)` for each call.
    """
    nearby = []
    pairs = []
    for i, c in enumerate(calls):
        hits = index.near(c, tolerance)
        own = [t for t in hits if carries(t["genotypes"][sample])]
        nearby.append((hits, own))
        pairs += [(distance(c, t), i, t["id"], t) for t in own]
    matched = {}
    used = set()
    for _, i, truth_id, t in sorted(pairs, key=lambda p: p[:3]):
        if i not in matched and truth_id not in used:
            matched[i] = t
            used.add(truth_id)
    rows = []
    for i, c in enumerate(calls):
        hits, own = nearby[i]
        t = matched.get(i)
        if t is not None:
            status = "tp" if same_family(c, t) else "wrong_family"
        elif own:
            status = "duplicate"
        elif any(missing(h["genotypes"][sample]) for h in hits):
            status = "unknown"
        elif any(any(carries(g) for g in h["genotypes"].values()) for h in hits):
            status = "other_sample"
        else:
            status = "fp"
        rows.append((c, status, t, own))
    return rows


FALSE = ("fp", "wrong_family", "duplicate")


def found(rows):
    """Ids of the truth MEIs matched by a call of the right family."""
    return {t["id"] for _, s, t, _ in rows if s == "tp"}


def count(rows):
    """Count the calls of each status, and in `false` the false calls of every kind."""
    status = collections.Counter(s for _, s, _, _ in rows)
    status["false"] = sum(status[s] for s in FALSE)
    return status


def main():
    """Evaluate the calls given on the command line."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--truth", required=True, help="MEI_Callset_GRCh38 CSV(.gz)")
    parser.add_argument("--bed", required=True, help="capture targets (BED)")
    parser.add_argument(
        "--chain", help="UCSC chain file lifting the truth to the assembly of the calls"
    )
    parser.add_argument("--padding", type=int, default=300, help="as in `scan` [300]")
    parser.add_argument(
        "--tolerance", type=int, default=20, help="breakpoint distance [20]"
    )
    parser.add_argument(
        "--min-qual",
        default="",
        help="comma-separated QUAL thresholds: precision and recall of the PASS calls "
        "with at least this QUAL",
    )
    parser.add_argument(
        "--fn-out", help="write the missed truth MEIs (PASS level) to this TSV"
    )
    parser.add_argument(
        "--fp-out",
        help="write the false PASS calls to this TSV, with a `status` column: fp, "
        "wrong_family (matched to an MEI of another family) or duplicate (a second "
        "call of an MEI already matched)",
    )
    parser.add_argument(
        "--tp-out",
        help="write the PASS calls that are not false to this TSV, with a "
        "`status` column: tp, other_sample (near an MEI of another sample) or "
        "unknown (a truth MEI nearby is not genotyped in the sample), and the truth "
        "genotype (0/1 or 1/1) of the matched MEI",
    )
    parser.add_argument(
        "--tolerances",
        default="",
        help="comma-separated breakpoint distances: precision and recall of the PASS "
        "calls at each, precise and imprecise (mate-only) calls apart",
    )
    parser.add_argument(
        "--bootstrap",
        type=int,
        default=2000,
        help="resamples of the samples for the confidence intervals; 0 to skip [2000]",
    )
    parser.add_argument(
        "runs",
        nargs="+",
        metavar="SAMPLE=FILE",
        help="calls of each sample (TSV or VCF); SAMPLE is the truth set column "
        "(NA24385 for HG002)",
    )
    args = parser.parse_args()

    runs = [r.split("=", 1) for r in args.runs]
    targets = load_targets(args.bed, args.padding)
    raw_targets = load_targets(args.bed, 0)
    chain = Chain(args.chain) if args.chain else None
    truth, unlifted = load_truth(
        args.truth, [s for s, _ in runs], targets, raw_targets, chain
    )
    thresholds = [float(q) for q in args.min_qual.split(",") if q]
    strata = [
        (0, 0, "in target"),
        (1, 100, "1-100 bp"),
        (101, args.padding, "padding >100 bp"),
    ]
    strata_total = collections.Counter()
    strata_any = collections.Counter()
    strata_pass = collections.Counter()
    total = collections.Counter()
    family_total = collections.Counter()
    family_found = collections.Counter()
    genotypes = collections.Counter()
    by_qual = collections.Counter()
    fn_rows, fp_rows, tp_rows = [], [], []
    per_sample = []
    lifted = f"; {unlifted} not lifted" if chain is not None else ""
    print(
        f"# truth MEIs in the padded targets (any HGSVC3 sample): {len(truth)}; "
        f"tolerance {args.tolerance} bp{lifted}"
    )
    print(
        "sample\ttruth\tfound_any\tfound_pass\tpass\ttp\tfp\twrong_family\tduplicate\tother_sample\tunknown"
        "\tprecision\trecall_pass\trecall_any\t"
        + "\t".join(f"tp_{t}/{t}" for t in TIERS)
    )
    index = TruthIndex(truth)
    sample_calls = []
    for sample, path in runs:
        calls = [
            c
            for c in load_calls(path)
            if any(in_targets(targets, c["contig"], p) for p in c["points"])
        ]
        sample_calls.append((sample, calls))
        carried = [t for t in truth if carries(t["genotypes"][sample])]
        found_any = found(evaluate(sample, calls, index, args.tolerance))
        passing = evaluate(sample, [c for c in calls if c["pass"]], index, args.tolerance)
        found_pass = found(passing)
        status = count(passing)
        per_tier = []
        for t in TIERS:
            in_tier = [s for c, s, _, _ in passing if c["confidence"] == t]
            per_tier.append(f"{in_tier.count('tp')}/{len(in_tier)}")
            total[f"tp_{t}"] += in_tier.count("tp")
            total[f"n_{t}"] += len(in_tier)
        for c, s, t, _ in passing:
            if s == "tp" and c["genotype"] != ".":
                genotypes[(unphased(t["genotypes"][sample]), c["genotype"])] += 1
        quals = {}
        for q in thresholds:
            kept = evaluate(sample, [c for c in calls if passes(c, q)], index, args.tolerance)
            kept_status = count(kept)
            quals[q] = (kept_status["tp"], kept_status["false"], len(found(kept)))
            by_qual[(q, "tp")] += quals[q][0]
            by_qual[(q, "fp")] += quals[q][1]
            by_qual[(q, "found")] += quals[q][2]
        if args.fn_out:
            near_truth = collections.defaultdict(set)
            for c, _, _, own in evaluate(sample, calls, index, args.tolerance):
                for o in own:
                    near_truth[o["id"]].add(c["confidence"] if c["pass"] else c["row"].get("filter", "filtered"))
        for t in carried:
            for lo, hi, name in strata:
                if lo <= t["distance"] <= hi:
                    strata_total[name] += 1
                    strata_any[name] += t["id"] in found_any
                    strata_pass[name] += t["id"] in found_pass
            family_total[t["family"]] += 1
            family_found[t["family"]] += t["id"] in found_pass
            if t["id"] not in found_pass and args.fn_out:
                fn_rows.append((sample, t["id"], t["family"], ",".join(sorted(near_truth[t["id"]])) or "."))
        for c, s, t, _ in passing:
            if s in FALSE:
                fp_rows.append((sample, s, c))
            else:
                tp_rows.append((sample, s, c, unphased(t["genotypes"][sample]) if t else "."))
        values = {
            "truth": len(carried),
            "found_any": len(found_any),
            "found_pass": len(found_pass),
            "pass": len(passing),
            "tp": status["tp"],
            "fp": status["fp"],
            "wrong_family": status["wrong_family"],
            "duplicate": status["duplicate"],
            "other_sample": status["other_sample"],
            "unknown": status["unknown"],
            "false": status["false"],
        }
        total.update(values)
        per_sample.append((values, quals))
        print(format_row(sample, values, per_tier))
    print(
        format_row(
            "all", total, [f"{total[f'tp_{t}']}/{total[f'n_{t}']}" for t in TIERS]
        )
    )
    if args.bootstrap and len(per_sample) > 1:
        print_intervals(per_sample, thresholds, args.bootstrap)
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
    if genotypes:
        right = sum(n for (t, c), n in genotypes.items() if t == c)
        print(
            f"\n# genotypes of the true PASS calls (truth, called): "
            f"{right}/{sum(genotypes.values())} right"
        )
        for (t, c), n in sorted(genotypes.items()):
            print(f"{t}\t{c}\t{n}")
    if thresholds:
        print("\n# PASS calls with QUAL >= threshold\ttp\tfalse\tprecision\trecall_pass")
        for q in thresholds:
            tp, fp = by_qual[(q, "tp")], by_qual[(q, "fp")]
            precision = f"{tp / (tp + fp):.3f}" if tp + fp else "."
            recall = by_qual[(q, "found")] / total["truth"] if total["truth"] else 0
            print(f"{q:g}\t{tp}\t{fp}\t{precision}\t{recall:.3f}")
    tolerances = [int(t) for t in args.tolerances.split(",") if t]
    if tolerances:
        print_tolerances(sample_calls, truth, index, tolerances)
    if args.fn_out:
        with open(args.fn_out, "w") as out:
            out.write("sample\ttruth_id\tfamily\tcall_tiers\n")
            for row in fn_rows:
                out.write("\t".join(row) + "\n")
    if args.tp_out:
        with open(args.tp_out, "w") as out:
            keys = list(tp_rows[0][2]["row"]) if tp_rows else []
            out.write("sample\tstatus\ttruth_genotype\t" + "\t".join(keys) + "\n")
            for sample, status, c, truth_genotype in tp_rows:
                out.write(
                    f"{sample}\t{status}\t{truth_genotype}\t"
                    + "\t".join(c["row"].get(k, ".") for k in keys)
                    + "\n"
                )
    if args.fp_out:
        with open(args.fp_out, "w") as out:
            keys = list(fp_rows[0][2]["row"]) if fp_rows else []
            out.write("sample\tstatus\t" + "\t".join(keys) + "\n")
            for sample, status, c in fp_rows:
                out.write(f"{sample}\t{status}\t" + "\t".join(c["row"].get(k, ".") for k in keys) + "\n")


def passes(call, q):
    """Tell whether a call is PASS with a QUAL of at least `q`."""
    return call["pass"] and call["qual"] is not None and call["qual"] >= q


def print_tolerances(sample_calls, truth, index, tolerances):
    """Print precision and recall of the PASS calls at each tolerance, precise and imprecise calls apart.

    Imprecise calls are the mate-only insertions (a `cipos` interval in the
    TSV, `CIPOS` in a VCF); recall counts the MEIs matched by either kind.
    """
    print("\n# PASS calls by tolerance (bp)\tcalls\ttp\tfalse\twrong_family\tduplicate\tprecision\trecall_pass")
    carried = sum(carries(t["genotypes"][s]) for s, _ in sample_calls for t in truth)
    for tolerance in tolerances:
        status = {kind: collections.Counter() for kind in ("all", "precise", "imprecise")}
        n_found = 0
        for sample, calls in sample_calls:
            rows = evaluate(sample, [c for c in calls if c["pass"]], index, tolerance)
            n_found += len(found(rows))
            for c, s, _, _ in rows:
                for kind in ("all", "imprecise" if c["imprecise"] else "precise"):
                    status[kind][s] += 1
                    status[kind]["calls"] += 1
        for kind, v in status.items():
            false = sum(v[s] for s in FALSE)
            precision = f"{v['tp'] / (v['tp'] + false):.3f}" if v["tp"] + false else "."
            recall = f"{n_found / carried:.3f}" if kind == "all" and carried else "."
            print(
                f"{tolerance} {kind}\t{v['calls']}\t{v['tp']}\t{false}\t{v['wrong_family']}"
                f"\t{v['duplicate']}\t{precision}\t{recall}"
            )


def print_intervals(per_sample, thresholds, resamples):
    """95% bootstrap intervals of precision and recall, resampling the samples.

    The samples are the unit (calls of one sample are not independent: depth,
    library and ancestry are shared), so the intervals are wide when few
    samples are evaluated. The seed is fixed: the same input gives the same
    intervals.
    """
    rng = random.Random(1)
    n = len(per_sample)

    def interval(stat):
        draws = []
        for _ in range(resamples):
            pick = [per_sample[rng.randrange(n)] for _ in range(n)]
            value = stat(pick)
            if value is not None:
                draws.append(value)
        draws.sort()
        if not draws:
            return ".", ".", "."
        return (
            f"{draws[int(0.025 * len(draws))]:.3f}",
            f"{draws[min(int(0.975 * len(draws)), len(draws) - 1)]:.3f}",
            f"{len(draws)}",
        )

    def ratio(num, den):
        return None if den == 0 else num / den

    print(f"\n# 95% bootstrap intervals over {n} samples ({resamples} resamples)")
    print("# level\tmetric\tlow\thigh")
    stats = [
        ("PASS", "precision", lambda p: ratio(sum(v["tp"] for v, _ in p), sum(v["tp"] + v["false"] for v, _ in p))),
        ("PASS", "recall", lambda p: ratio(sum(v["found_pass"] for v, _ in p), sum(v["truth"] for v, _ in p))),
    ]
    for q in thresholds:
        stats.append(
            (
                f"QUAL>={q:g}",
                "precision",
                lambda p, q=q: ratio(sum(x[q][0] for _, x in p), sum(x[q][0] + x[q][1] for _, x in p)),
            )
        )
        stats.append(
            (
                f"QUAL>={q:g}",
                "recall",
                lambda p, q=q: ratio(sum(x[q][2] for _, x in p), sum(v["truth"] for v, _ in p)),
            )
        )
    for level, metric, stat in stats:
        low, high, _ = interval(stat)
        print(f"{level}\t{metric}\t{low}\t{high}")


def format_row(name, v, per_tier):
    """One line of the results table."""
    decided = v["tp"] + v["false"]
    precision = f"{v['tp'] / decided:.3f}" if decided else "."
    recall_pass = f"{v['found_pass'] / v['truth']:.3f}" if v["truth"] else "."
    recall_any = f"{v['found_any'] / v['truth']:.3f}" if v["truth"] else "."
    return "\t".join(
        [name]
        + [
            str(v[k])
            for k in (
                "truth",
                "found_any",
                "found_pass",
                "pass",
                "tp",
                "fp",
                "wrong_family",
                "duplicate",
                "other_sample",
                "unknown",
            )
        ]
        + [precision, recall_pass, recall_any]
        + per_tier
    )


if __name__ == "__main__":
    main()
