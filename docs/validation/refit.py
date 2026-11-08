#!/usr/bin/env python3
"""Refit the quality score of `genotype::score` on labelled PASS calls.

The score is a naive Bayes model: a prior and, for each of seven pieces of
evidence, the log10 likelihood ratios (true against false insertions) when it
is present and when it is absent, with one call added to each count. The
output lists the values to put in `src/genotype/score.rs` and the quality of
the ranking: AUC leaving one sample out, AUC in sample, precision and recall
at several QUAL thresholds, and the AUC on a held-out sample.

Labels come from `evaluate.py` run on the same calls: the PASS calls of its
`--fp-out` file are false; those of its `--tp-out` file are true (true
positives, and the few calls near an MEI of another sample or of a sample not
genotyped there). A PASS call in neither file (outside the padded targets) is
left out. The test samples and the held-out sample (HG002) are optional.

Usage:
    refit.py --fp fp.tsv --tp tp.tsv SAMPLE=calls.tsv ...
             [--test SAMPLE=calls.tsv ... --test-fp fp.tsv --test-tp tp.tsv]
             [--held-out SAMPLE=calls.tsv --held-out-fp fp.tsv --held-out-tp tp.tsv]

Only the Python standard library is used.
"""

import argparse
import csv
import math


def read_calls(sample, path, false_calls, true_calls):
    """PASS calls of a `scan` TSV, labelled 1 (true) or 0 (false); others are left out."""
    rows, columns = [], None
    with open(path) as handle:
        for fields in csv.reader(handle, delimiter="\t"):
            if not fields:
                continue
            if fields[0] == "#id":
                columns = fields
                continue
            if fields[0].startswith("#"):
                continue
            row = dict(zip(columns, fields))
            if row["filter"] != "PASS":
                continue
            row["sample"] = sample
            key = (sample, row["contig"], row["position"])
            if key in false_calls:
                row["label"] = 0
            elif key in true_calls:
                row["label"] = 1
            else:
                continue
            rows.append(row)
    return rows


def read_keys(path, default_sample=None):
    """`(sample, contig, position)` of the calls of an `--fp-out` or `--tp-out` file."""
    with open(path) as handle:
        return {
            (r.get("sample", default_sample), r["contig"], r["position"])
            for r in csv.DictReader(handle, delimiter="\t")
        }


def number(text, default=0.0):
    """A float, or `default` for `.`."""
    try:
        return float(text)
    except ValueError:
        return default


def features(row):
    """Presence of each piece of evidence (the thresholds of `genotype::score`)."""
    support = max(int(row["support"]), 1)
    checked = int(row["checked_mates"] or 0)
    mei = int(row["mei_mates"] or 0)
    identity = number(row["identity"]) if row["family"] != "." else 0
    return {
        "both": row["left_junction"] != "." and row["right_junction"] != ".",
        "tail": int(row["polya_clips"]) + int(row["polyt_clips"]) > 0,
        "tsd": row["tsd_length"] != ".",
        "ident": identity >= 0.98,
        "starts": int(row["distinct_starts"]) / support >= 0.7,
        "mei_mates": checked > 0 and mei / max(checked, 1) >= 0.6,
        "ratio": row["signal_ratio"] != "." and number(row["signal_ratio"]) >= 0.2,
        # Evidence of the mate-only calls, tried and left out (#67).
        "by_mates": row["typed_by"] == "mates",
        "bp_clips": int(row["breakpoint_clips"] or 0) > 0,
        "imprecise": row["cipos"] != ".",
        "junction_clips": int(row["left_junction_clips"]) + int(row["right_junction_clips"]) >= 2,
        "mei_many": mei >= 8,
    }


SCORE = ["both", "tail", "tsd", "ident", "starts", "mei_mates", "ratio"]
EXTRA = ["by_mates", "bp_clips", "imprecise", "junction_clips", "mei_many"]


def fit(rows, names):
    """Prior and likelihood ratios, one call added to each count."""
    true = [r for r in rows if r["label"] == 1]
    false = [r for r in rows if r["label"] == 0]
    prior = math.log10((len(true) + 1) / (len(false) + 1))
    ratios = {}
    for name in names:
        p_true = (sum(r["f"][name] for r in true) + 1) / (len(true) + 2)
        p_false = (sum(r["f"][name] for r in false) + 1) / (len(false) + 2)
        ratios[name] = (
            math.log10(p_true / p_false),
            math.log10((1 - p_true) / (1 - p_false)),
        )
    return prior, ratios


def score(model, row, names):
    """log10 odds of a true insertion."""
    prior, ratios = model
    return prior + sum(ratios[n][0] if row["f"][n] else ratios[n][1] for n in names)


def auc(rows, scores):
    """Probability that a true call scores above a false one."""
    pos = [s for r, s in zip(rows, scores) if r["label"] == 1]
    neg = [s for r, s in zip(rows, scores) if r["label"] == 0]
    if not pos or not neg:
        return float("nan")
    return sum((p > n) + 0.5 * (p == n) for p in pos for n in neg) / (len(pos) * len(neg))


def phred(log_odds):
    """QUAL of a score, as `genotype::quality` (capped at 60)."""
    return min(round(-10 * math.log10(max(1 / (1 + 10**log_odds), 1e-12))), 60)


def leave_one_out(rows, samples, names):
    """Scores of each sample from a model fitted on the others."""
    scored = []
    for sample in samples:
        model = fit([r for r in rows if r["sample"] != sample], names)
        scored += [(r, score(model, r, names)) for r in rows if r["sample"] == sample]
    return [r for r, _ in scored], [s for _, s in scored]


def table(label, rows, scores):
    """True / false PASS calls kept at each QUAL threshold."""
    cells = []
    for q in (0, 2, 5, 8, 10):
        tp = sum(1 for r, s in zip(rows, scores) if r["label"] == 1 and phred(s) >= q)
        fp = sum(1 for r, s in zip(rows, scores) if r["label"] == 0 and phred(s) >= q)
        cells.append(f"q>={q}: {tp}/{fp} ({tp / max(tp + fp, 1):.3f})")
    print(f"{label:22s}", " | ".join(cells))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--fp", required=True, help="--fp-out of evaluate.py on the samples")
    parser.add_argument("--tp", required=True, help="--tp-out of evaluate.py on the samples")
    parser.add_argument("--test", action="append", default=[], metavar="SAMPLE=FILE",
                        help="calls of a test sample, not used to fit (repeat)")
    parser.add_argument("--test-fp", help="--fp-out of evaluate.py on the test samples")
    parser.add_argument("--test-tp", help="--tp-out of evaluate.py on the test samples")
    parser.add_argument("--held-out", metavar="SAMPLE=FILE", help="calls of a sample not used to fit")
    parser.add_argument("--held-out-fp", help="--fp-out of evaluate.py on the held-out sample")
    parser.add_argument("--held-out-tp", help="--tp-out of evaluate.py on the held-out sample")
    parser.add_argument("runs", nargs="+", metavar="SAMPLE=FILE")
    args = parser.parse_args()

    runs = [r.split("=", 1) for r in args.runs]
    samples = [s for s, _ in runs]
    false_calls, true_calls = read_keys(args.fp), read_keys(args.tp)
    rows = [r for s, p in runs for r in read_calls(s, p, false_calls, true_calls)]
    for r in rows:
        r["f"] = features(r)
    print(f"PASS {len(rows)} true {sum(r['label'] for r in rows)} "
          f"false {sum(1 - r['label'] for r in rows)} ({len(samples)} samples)")

    print("\n# AUC of the score and of variants with more evidence")
    variants = [("score", SCORE)] + [(f"+{e}", SCORE + [e]) for e in EXTRA] + [("+all", SCORE + EXTRA)]
    for label, names in variants:
        model = fit(rows, names)
        in_sample = auc(rows, [score(model, r, names) for r in rows])
        loo_rows, loo_scores = leave_one_out(rows, samples, names)
        print(f"{label:16s} leave-one-sample-out {auc(loo_rows, loo_scores):.3f}  in sample {in_sample:.3f}")

    print("\n# true / false PASS calls kept at each QUAL")
    model = fit(rows, SCORE)
    table("refit, in sample", rows, [score(model, r, SCORE) for r in rows])
    loo_rows, loo_scores = leave_one_out(rows, samples, SCORE)
    table("refit, leave-one-out", loo_rows, loo_scores)

    if args.test:
        test_fp, test_tp = read_keys(args.test_fp), read_keys(args.test_tp)
        test = [r for s, p in (x.split("=", 1) for x in args.test) for r in read_calls(s, p, test_fp, test_tp)]
        for r in test:
            r["f"] = features(r)
        test_scores = [score(model, r, SCORE) for r in test]
        current = [number(r["qual"]) for r in test]
        print(f"\n# test ({len(args.test)} samples): PASS {len(test)} true {sum(r['label'] for r in test)}, "
              f"AUC {auc(test, test_scores):.3f} (QUAL of the scans: {auc(test, current):.3f})")
        table("test", test, test_scores)

    if args.held_out:
        sample, path = args.held_out.split("=", 1)
        held = read_calls(
            sample, path, read_keys(args.held_out_fp, sample), read_keys(args.held_out_tp, sample)
        )
        for r in held:
            r["f"] = features(r)
        held_scores = [score(model, r, SCORE) for r in held]
        print(f"\n# held out {sample}: PASS {len(held)} true {sum(r['label'] for r in held)}, "
              f"AUC {auc(held, held_scores):.3f}")
        table(sample, held, held_scores)

    true = [r for r in rows if r["label"] == 1]
    false = [r for r in rows if r["label"] == 0]
    print("\n# values for src/genotype/score.rs")
    print(f"PRIOR = {model[0]:.3f}")
    for name in SCORE:
        present, absent = model[1][name]
        print(f"{name:10s} ({sum(r['f'][name] for r in true)}/{len(true)}, "
              f"{sum(r['f'][name] for r in false)}/{len(false)})  log_ratios: ({present:.3f}, {absent:.3f})")


if __name__ == "__main__":
    main()
