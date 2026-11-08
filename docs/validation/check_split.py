#!/usr/bin/env python3
"""Check that no family of `split.tsv` spans the tuning and the test sets.

`split.tsv` gives each sample a `set` (`tuning`, `test`, or `related` for a
sample left out of both because a relative is in the tuning set) and the
`family_id` of the 1000 Genomes pedigree (`20130606_g1k.ped`, see MANIFEST.md;
an unrelated sample is its own family). Relatives share inherited insertions,
so a family across sets lets the tuning set see part of the test set.

With `--pedigree`, the `family_id` column is also compared with the pedigree.
Exit code 1 and one line per problem when the check fails.

Only the Python standard library is used.
"""

import argparse
import collections
import csv
import sys

SETS = ("tuning", "test", "related")


def read_split(path):
    """Rows of `split.tsv` as dicts."""
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_pedigree(path):
    """Family of each individual of a 1000 Genomes `.ped` file (tab-separated, with a header)."""
    with open(path) as handle:
        return {row["Individual ID"]: row["Family ID"] for row in csv.DictReader(handle, delimiter="\t")}


def split_leaks(rows):
    """Problems of the split: unknown sets, families in tuning and test, related samples without kin in tuning."""
    problems = []
    sets = collections.defaultdict(set)
    for row in rows:
        if row["set"] not in SETS:
            problems.append(f"{row['sample']}: unknown set {row['set']!r}")
        sets[row["family_id"]].add(row["set"])
    for family, used in sorted(sets.items()):
        if {"tuning", "test"} <= used:
            members = ", ".join(f"{r['sample']} ({r['set']})" for r in rows if r["family_id"] == family)
            problems.append(f"family {family} in the tuning and the test sets: {members}")
        if "related" in used and "tuning" not in used:
            problems.append(f"family {family}: a `related` sample without a relative in the tuning set")
    return problems


def pedigree_mismatches(rows, pedigree):
    """List the samples whose `family_id` differs from the pedigree."""
    return [
        f"{row['sample']}: family_id {row['family_id']}, pedigree {pedigree.get(row['sample'], 'absent')}"
        for row in rows
        if pedigree.get(row["sample"]) != row["family_id"]
    ]


def main():
    """Print the problems of the split and exit with 1 when there is any."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("split", help="split.tsv")
    parser.add_argument("--pedigree", help="20130606_g1k.ped of the 1000 Genomes project")
    args = parser.parse_args()
    rows = read_split(args.split)
    problems = split_leaks(rows)
    if args.pedigree:
        problems += pedigree_mismatches(rows, read_pedigree(args.pedigree))
    for problem in problems:
        print(problem, file=sys.stderr)
    counts = collections.Counter(row["set"] for row in rows)
    print(f"split: {', '.join(f'{counts[s]} {s}' for s in SETS)}; {len(problems)} problems")
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
