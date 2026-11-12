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
import math
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
                    "carriers": sum(carries(row[col[s]]) for s in others),
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


def concordant(truth, called):
    """Whether a call agrees with the unphased truth genotype of a true insertion.

    A diploid call must equal it; a haploid `1` (the X or Y of a male, #127) agrees with any
    genotype that carries the insertion, and a haploid `0` never does.
    """
    if called in ("0", "1"):
        return called == "1"
    return truth == called


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
                    "gq": int(row["gq"]) if row.get("gq", ".") != "." else None,
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
            genotype, gq = ".", None
            if len(fields) > 9:
                sample = dict(zip(fields[8].split(":"), fields[9].split(":")))
                genotype = sample.get("GT", ".").replace("|", "/")
                gq = int(sample["GQ"]) if sample.get("GQ", ".") != "." else None
            calls.append(
                {
                    "contig": contig_name(fields[0]),
                    "points": points,
                    "pass": fields[6] in ("PASS", "."),
                    "confidence": info.get("CONFIDENCE", "."),
                    "qual": None if fields[5] == "." else float(fields[5]),
                    "genotype": genotype,
                    "gq": gq,
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
        help="resamples of the families and loci for the confidence intervals; 0 to skip [2000]",
    )
    parser.add_argument(
        "--annotation",
        help="GENCODE GTF(.gz) of the assembly of the calls: recall by distance to the "
        "coding sequences and exons (strata, with --strata-qual)",
    )
    parser.add_argument(
        "--strata-qual",
        type=float,
        default=10,
        help="QUAL threshold of the stratified recall: that of the default profile [10]",
    )
    parser.add_argument(
        "--families",
        help="split.tsv: the `family_id` of each sample, resampled together in the "
        "bootstrap (default: every sample is its own family)",
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
        levels = {"PASS": passing}
        for q in thresholds:
            kept = evaluate(sample, [c for c in calls if passes(c, q)], index, args.tolerance)
            levels[f"QUAL>={q:g}"] = kept
            kept_status = count(kept)
            by_qual[(q, "tp")] += kept_status["tp"]
            by_qual[(q, "fp")] += kept_status["false"]
            by_qual[(q, "found")] += len(found(kept))
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
        per_sample.append((sample, carried, levels))
        print(format_row(sample, values, per_tier))
    print(
        format_row(
            "all", total, [f"{total[f'tp_{t}']}/{total[f'n_{t}']}" for t in TIERS]
        )
    )
    loci = Loci(per_sample)
    print_units(per_sample, loci)
    families = read_families(args.families) if args.families else {}
    if args.bootstrap and len(per_sample) > 1:
        print_intervals(per_sample, loci, families, args.bootstrap)
    annotation = load_annotation(args.annotation) if args.annotation else None
    samples = read_samples(args.families) if args.families else {}
    print_strata(per_sample, annotation, samples, args.strata_qual, loci, families, args.bootstrap)
    print_gq(per_sample, args.strata_qual)
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
        right = sum(n for (t, c), n in genotypes.items() if concordant(t, c))
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


# Truth MEIs and false calls closer than this (bp) are one locus across samples.
LOCUS_GAP = 100


class Loci:
    """Loci of the evaluation: truth MEIs and false PASS calls of all samples, clustered by position.

    An inherited insertion or a recurrent artefact is one locus seen in many
    samples; the clusters join positions less than `LOCUS_GAP` apart on a contig.
    """

    def __init__(self, per_sample):
        """Cluster the carried MEIs and the false calls of every level of `per_sample`."""
        points = set()
        for _, carried, levels in per_sample:
            points.update((t["contig"], t["pos"]) for t in carried)
            for rows in levels.values():
                points.update((c["contig"], min(c["points"])) for c, s, _, _ in rows if s in FALSE)
        self.index = {}
        previous, locus = None, -1
        for contig, pos in sorted(points):
            if previous is None or previous[0] != contig or pos - previous[1] > LOCUS_GAP:
                locus += 1
            self.index[(contig, pos)] = locus
            previous = (contig, pos)
        self.count = locus + 1

    def of_truth(self, t):
        """Locus of a truth MEI."""
        return self.index[(t["contig"], t["pos"])]

    def of_call(self, c):
        """Locus of a false call."""
        return self.index[(c["contig"], min(c["points"]))]


def outcomes(per_sample, loci, level):
    """Count the true calls, false calls and missed MEIs per sample and locus at one level.

    Returns `{(sample, locus): [tp, false, fn]}`; a true call counts at the
    locus of its MEI, so that `tp + fn` is the number of MEIs carried there.
    """
    cells = collections.defaultdict(lambda: [0, 0, 0])
    for sample, carried, levels in per_sample:
        rows = levels[level]
        hit = found(rows)
        for t in carried:
            cells[(sample, loci.of_truth(t))][0 if t["id"] in hit else 2] += 1
        for c, s, _, _ in rows:
            if s in FALSE:
                cells[(sample, loci.of_call(c))][1] += 1
    return cells


def read_samples(path):
    """Read the rows of a `split.tsv`, by sample."""
    with open(path) as handle:
        return {r["sample"]: r for r in csv.DictReader(handle, delimiter="\t")}


def load_annotation(path):
    """Read the merged exons and coding sequences of a GTF, per contig (0-based half-open, sorted)."""
    raw = {"exon": collections.defaultdict(list), "CDS": collections.defaultdict(list)}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 5)
            if fields[2] in raw:
                raw[fields[2]][contig_name(fields[0])].append((int(fields[3]) - 1, int(fields[4])))
    merged = {}
    for feature, contigs in raw.items():
        merged[feature] = {}
        for contig, intervals in contigs.items():
            intervals.sort()
            out = [list(intervals[0])]
            for start, end in intervals[1:]:
                if start <= out[-1][1]:
                    out[-1][1] = max(out[-1][1], end)
                else:
                    out.append([start, end])
            merged[feature][contig] = ([a for a, _ in out], [b for _, b in out])
    return merged


def clopper_pearson_lower(found, total, alpha=0.05):
    """Return the one-sided exact lower bound of a proportion (Clopper-Pearson).

    It is the sensitivity that can be claimed with confidence `1 - alpha`
    from `found` detections among `total`: the `p` at which observing
    `found` or more has probability `alpha`.
    """
    if found == 0:
        return 0.0

    def log_comb(k):
        return math.lgamma(total + 1) - math.lgamma(k + 1) - math.lgamma(total - k + 1)

    def upper_tail(p):
        terms = (log_comb(k) + k * math.log(p) + (total - k) * math.log1p(-p) for k in range(found, total + 1))
        return sum(math.exp(term) for term in terms)

    low, high = 0.0, 1.0
    for _ in range(60):
        middle = (low + high) / 2
        if upper_tail(middle) < alpha:
            low = middle
        else:
            high = middle
    return low


def regularized_beta(x, a, b):
    """Regularized incomplete beta function I_x(a, b), by its continued fraction (Lentz)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    if x > (a + 1) / (a + b + 2):
        return 1.0 - regularized_beta(1.0 - x, b, a)
    front = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log1p(-x)) / a
    tiny = 1e-300
    c, d = 1.0, 1.0 - (a + b) * x / (a + 1)
    d = 1.0 / (d if abs(d) > tiny else tiny)
    f = d
    for m in range(1, 300):
        for numerator in (
            m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m)),
            -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1)),
        ):
            d = 1.0 + numerator * d
            d = 1.0 / (d if abs(d) > tiny else tiny)
            c = 1.0 + numerator / c
            c = c if abs(c) > tiny else tiny
            f *= c * d
        if abs(c * d - 1.0) < 1e-12:
            break
    return front * f


def clustered_lower(found, total, variance, sites, alpha=0.05):
    """Return a one-sided lower bound of a proportion of clustered units (Korn and Graubard 1998).

    The carrier pairs of one site, and of one family, are not independent: the
    Clopper-Pearson bound is taken on the effective number of pairs,
    `p (1 - p) / variance` with `variance` the bootstrap variance of the
    proportion over families and loci, at most `total`. When the proportion
    is 0 or 1 the variance tells nothing, and the number of distinct `sites`
    is used instead. The bound is the `alpha` quantile of a beta distribution,
    which allows a non-integer number of successes.
    """
    if total == 0 or found == 0:
        return 0.0
    p = found / total
    n_eff = min(total, p * (1 - p) / variance) if 0 < p < 1 and variance > 0 else min(total, sites)
    a, b = p * n_eff, n_eff * (1 - p) + 1
    low, high = 0.0, 1.0
    for _ in range(60):
        middle = (low + high) / 2
        if regularized_beta(middle, a, b) < alpha:
            low = middle
        else:
            high = middle
    return low


def stratum_variances(cells, groups, n_loci, resamples):
    """Bootstrap variance of the recall of each stratum, over families and loci.

    `cells` maps a stratum to its `(group, locus, found)` carrier pairs; each
    resample draws the groups and the loci with replacement (as `bootstrap`),
    the same draws for every stratum. The seed is fixed.
    """
    rng = random.Random(1)
    sums = collections.defaultdict(lambda: [0, 0.0, 0.0])
    for _ in range(resamples):
        group_weight = collections.Counter(rng.choice(groups) for _ in groups)
        locus_weight = collections.Counter(rng.randrange(n_loci) for _ in range(n_loci))
        for stratum, pairs in cells.items():
            hit = total = 0
            for group, locus, found_it in pairs:
                weight = group_weight[group] * locus_weight[locus]
                hit += weight * found_it
                total += weight
            if total:
                value = hit / total
                acc = sums[stratum]
                acc[0] += 1
                acc[1] += value
                acc[2] += value * value
    return {k: max(0.0, sq / n - (sm / n) ** 2) for k, (n, sm, sq) in sums.items() if n > 1}


def bins(value, edges, unit=""):
    """Name of the bin of `value` among the half-open bins split at `edges`."""
    for low, high in zip([None] + edges, edges + [None]):
        if high is None or value < high:
            if low is None:
                return f"< {high}{unit}"
            return f">= {low}{unit}" if high is None else f"{low}-{high - 1}{unit}"
    return "."


def print_strata(per_sample, annotation, samples, strata_qual, loci=None, families=None, resamples=0):
    """Print the recall per stratum at the default QUAL, with one-sided 95% lower bounds.

    The unit is the carrier pair (a truth MEI in one sample). Strata: distance
    to the target, to a coding sequence and to an exon (with `annotation`),
    family, number of HGSVC3 carriers, median candidate depth and read length
    of the sample (with `samples`, from `split.tsv`). The genotype errors of
    the true calls are counted by depth. Each stratum gets its number of
    distinct sites and two lower bounds: Clopper-Pearson on the pairs, which
    treats them as independent and is optimistic for strata of few sites,
    and the clustered bound (`clustered_lower`, with `resamples` two-way
    bootstrap draws over the `families` and the `loci`, #153), the
    sensitivity that can be claimed for the stratum.
    """
    level = f"QUAL>={strata_qual:g}"
    if level not in per_sample[0][2]:
        return
    strata = collections.defaultdict(lambda: [0, 0, 0])
    sites = collections.defaultdict(set)
    cells = collections.defaultdict(list)
    families = families or {}
    genotypes = collections.defaultdict(lambda: [0, 0])
    kinds = collections.Counter()
    for sample, carried, levels in per_sample:
        rows = levels[level]
        hit = found(rows)
        info = samples.get(sample, {})
        for c, s, t, _ in rows:
            if s == "tp":
                kinds["imprecise (mate-only)" if c["imprecise"] else "precise"] += 1
                if c["genotype"] != ".":
                    depth = bins(int(info["median_candidate_depth"]), [100, 150]) if info else "all"
                    genotypes[depth][0] += concordant(unphased(t["genotypes"][sample]), c["genotype"])
                    genotypes[depth][1] += 1
        for t in carried:
            if t["distance"] == 0:
                where = "in target"
            else:
                where = "1-100 bp from it" if t["distance"] <= 100 else "101-300 bp from it"
            names = [
                ("target", where),
                ("family", t["family"]),
                ("HGSVC3 carriers (of 65)", bins(t["carriers"], [4, 17, 33])),
            ]
            if annotation is not None:
                cds = target_distance(annotation["CDS"], t["contig"], t["pos"])
                exon = target_distance(annotation["exon"], t["contig"], t["pos"])
                names.append(("coding sequence", "." if cds is None else ("CDS +/- 20 bp" if cds <= 20 else "farther")))
                names.append(("exon", "." if exon is None else ("in an exon" if exon == 0 else "outside")))
            if info:
                names.append(("median candidate depth", bins(int(info["median_candidate_depth"]), [100, 150])))
                names.append(("read length", f"{info['read_length']} bp"))
            for group, name in names:
                cell = strata[(group, name)]
                cell[0] += t["id"] in hit
                cell[1] += 1
                sites[(group, name)].add(t["id"])
                if loci is not None:
                    cells[(group, name)].append((families.get(sample, sample), loci.of_truth(t), t["id"] in hit))
    variances = {}
    if loci is not None and resamples:
        groups = sorted({families.get(sample, sample) for sample, _, _ in per_sample})
        variances = stratum_variances(cells, groups, loci.count, resamples)
    print(
        f"\n# recall at {level} by stratum (carrier pairs)\tfound\ttotal\tsites\trecall\tlower_95\tlower_95_clustered"
    )
    for key, (n_found, total, _) in sorted(strata.items()):
        group, name = key
        clustered = "."
        if key in variances:
            clustered = f"{clustered_lower(n_found, total, variances[key], len(sites[key])):.3f}"
        print(
            f"{group}: {name}\t{n_found}\t{total}\t{len(sites[key])}\t{n_found / total:.3f}"
            f"\t{clopper_pearson_lower(n_found, total):.3f}\t{clustered}"
        )
    print(f"\n# true calls at {level} by kind\t" + "\t".join(f"{k} {n}" for k, n in sorted(kinds.items())))
    if genotypes:
        print(f"\n# genotypes of the true calls at {level} by median candidate depth\tright\tcalls")
        for depth, (right, n) in sorted(genotypes.items()):
            print(f"{depth}\t{right}\t{n}")


GQ_THRESHOLDS = (0, 10, 20, 30)
GQ_BINS = [10, 20, 30, 50]


def gq_table(rows):
    """Genotype quality against the truth genotypes, from `(status, truth, genotype, gq)` rows.

    Returns two lists. Kept: for each GQ threshold of `GQ_THRESHOLDS`, the
    true calls kept, the false calls kept and the genotypes right among the
    true calls kept (calls of unknown status are left out). Calibration: for each GQ bin of `GQ_BINS`, the true
    calls, the wrong genotypes, the observed error rate and the error rate
    the GQ claims (mean of 10^(-GQ/10)). Calls without genotype count as
    GQ 0 and as wrong.
    """
    true = [(t, g, q or 0) for s, t, g, q in rows if s == "tp"]
    false = [q or 0 for s, _, _, q in rows if s in FALSE]
    kept = []
    for threshold in GQ_THRESHOLDS:
        held = [(t, g) for t, g, q in true if q >= threshold]
        right = sum(g != "." and concordant(t, g) for t, g in held)
        kept.append((threshold, len(held), len(true), sum(q >= threshold for q in false), len(false), right))
    groups = collections.defaultdict(list)
    for t, g, q in true:
        groups[bins(q, GQ_BINS)].append((g != "." and concordant(t, g), q))
    calibration = []
    for name, members in sorted(groups.items(), key=lambda kv: min(q for _, q in kv[1])):
        wrong = sum(not ok for ok, _ in members)
        claimed = sum(10 ** (-q / 10) for _, q in members) / len(members)
        calibration.append((name, len(members), wrong, wrong / len(members), claimed))
    return kept, calibration


def print_gq(per_sample, strata_qual):
    """Print the GQ table (`gq_table`) of the calls at the default QUAL (#156)."""
    level = f"QUAL>={strata_qual:g}"
    if level not in per_sample[0][2]:
        return
    rows = [
        (s, unphased(t["genotypes"][sample]) if s == "tp" else ".", c["genotype"], c.get("gq"))
        for sample, _, levels in per_sample
        for c, s, t, _ in levels[level]
    ]
    if not any(q is not None for _, _, _, q in rows):
        return
    kept, calibration = gq_table(rows)
    print(f"\n# calls at {level} kept at GQ >= threshold\ttrue\tfalse\tgenotypes_right")
    for threshold, held, n_true, false, n_false, right in kept:
        concordance = f"{right / held:.3f}" if held else "."
        print(f"{threshold}\t{held}/{n_true}\t{false}/{n_false}\t{right}/{held} ({concordance})")
    print(f"\n# GQ calibration of the true calls at {level}\tcalls\twrong\twrong_rate\tclaimed_rate")
    for name, n, wrong, observed, claimed in calibration:
        print(f"{name}\t{n}\t{wrong}\t{observed:.3f}\t{claimed:.3f}")


def read_families(path):
    """Read the `family_id` of each sample of a `split.tsv`."""
    with open(path) as handle:
        return {r["sample"]: r["family_id"] for r in csv.DictReader(handle, delimiter="\t")}


def print_units(per_sample, loci):
    """Print the carrier pairs and distinct sites behind the counts, and the recall per site.

    A truth MEI is a site; each sample that carries it is a carrier pair.
    Recall per pair is the usual recall; per site, a site counts as found
    when at least one carrier has a PASS call of it.
    """
    pairs = sum(len(carried) for _, carried, _ in per_sample)
    sites = {t["id"] for _, carried, _ in per_sample for t in carried}
    found_sites = {i for _, _, levels in per_sample for i in found(levels["PASS"])}
    every = collections.Counter()
    for _, carried, levels in per_sample:
        hit = found(levels["PASS"])
        every.update(t["id"] for t in carried if t["id"] not in hit)
    false_loci = {loci.of_call(c) for _, _, levels in per_sample for c, s, _, _ in levels["PASS"] if s in FALSE}
    false_calls = sum(count(levels["PASS"])["false"] for _, _, levels in per_sample)
    print(f"\n# units: {pairs} carrier pairs at {len(sites)} truth sites; {false_calls} false PASS calls "
          f"at {len(false_loci)} loci")
    if sites:
        print(f"sites found (PASS) in at least one carrier\t{len(found_sites)}/{len(sites)}"
              f"\t{len(found_sites) / len(sites):.3f}")
        print(f"sites found (PASS) in every carrier\t{len(sites) - len(every)}/{len(sites)}"
              f"\t{(len(sites) - len(every)) / len(sites):.3f}")


def print_intervals(per_sample, loci, families, resamples):
    """Print 95% two-way bootstrap intervals of precision and recall.

    The families of the samples (`families`; an unlisted sample is its own
    family) and the loci are resampled independently, and each sample-locus
    cell is weighted by the product of the two counts (pigeonhole bootstrap,
    Owen 2007): relatives share insertions, and one locus is called or missed
    in many samples, so a bootstrap over the samples alone treats dependent
    calls as independent and gives intervals that are too narrow. That one is
    printed too, for comparison. The seed is fixed: the same input gives the
    same intervals.
    """
    family_of = {sample: families.get(sample, sample) for sample, _, _ in per_sample}
    groups = sorted(set(family_of.values()))
    samples = [sample for sample, _, _ in per_sample]
    print(
        f"\n# 95% bootstrap intervals ({resamples} resamples): families x loci ({len(groups)} families, "
        f"{loci.count} loci), and samples alone ({len(samples)})"
    )
    print("# level\tmetric\tvalue\tlow\thigh\tlow_samples\thigh_samples")
    for level in per_sample[0][2]:
        cells = [(family_of[s], s, locus, v) for (s, locus), v in outcomes(per_sample, loci, level).items()]
        two_way = bootstrap(cells, groups, loci.count, resamples, by_locus=True)
        by_sample = bootstrap([(s, s, locus, v) for _, s, locus, v in cells], samples, loci.count, resamples)
        tp = sum(v[0] for *_, v in cells)
        false = sum(v[1] for *_, v in cells)
        fn = sum(v[2] for *_, v in cells)
        for metric, value in (("precision", ratio(tp, tp + false)), ("recall", ratio(tp, tp + fn))):
            print(
                "\t".join(
                    [level, metric, "." if value is None else f"{value:.3f}"]
                    + list(two_way[metric])
                    + list(by_sample[metric])
                )
            )


def ratio(num, den):
    """Divide, or return None when `den` is 0."""
    return None if den == 0 else num / den


def bootstrap(cells, groups, n_loci, resamples, by_locus=False):
    """Draw the 2.5th and 97.5th percentiles of precision and recall.

    `cells` are `(group, sample, locus, [tp, false, fn])`; every resample draws
    the groups with replacement and, with `by_locus`, the loci too.
    """
    rng = random.Random(1)
    draws = {"precision": [], "recall": []}
    for _ in range(resamples):
        group_weight = collections.Counter(rng.choice(groups) for _ in groups)
        locus_weight = collections.Counter(rng.randrange(n_loci) for _ in range(n_loci)) if by_locus else None
        tp = false = fn = 0
        for group, _, locus, (c_tp, c_false, c_fn) in cells:
            weight = group_weight[group] * (locus_weight[locus] if by_locus else 1)
            if weight:
                tp += weight * c_tp
                false += weight * c_false
                fn += weight * c_fn
        for metric, value in (("precision", ratio(tp, tp + false)), ("recall", ratio(tp, tp + fn))):
            if value is not None:
                draws[metric].append(value)
    out = {}
    for metric, values in draws.items():
        values.sort()
        out[metric] = (
            (f"{values[int(0.025 * len(values))]:.3f}", f"{values[min(int(0.975 * len(values)), len(values) - 1)]:.3f}")
            if values
            else (".", ".")
        )
    return out


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
