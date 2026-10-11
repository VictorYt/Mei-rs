#!/usr/bin/env python3
"""Download the default MEI consensus sequences of mei-rs from Dfam.

Writes mei_consensus.fa (see README.md): the human Alu (AluY, AluS and AluJ
subfamilies), young LINE-1 (L1HS 5' and 3' ends, L1PA2-L1PA8 3' ends, and
the L1P1 ORF2 model for the middle of the element) and SVA (SVA_A-SVA_F)
consensus sequences of Dfam, released under CC0 1.0.

Headers: `>name#class accession`, class `SINE/Alu`, `LINE/L1` or
`Retroposon/SVA`. Only the Python standard library is used.
"""

import json
import sys
import urllib.request

API = "https://dfam.org/api/families"
GROUPS = [
    ("AluY", "SINE/Alu"),
    ("AluS", "SINE/Alu"),
    ("AluJ", "SINE/Alu"),
    ("L1HS", "LINE/L1"),
    ("L1PA", "LINE/L1"),
    ("L1P1", "LINE/L1"),
    ("SVA_", "Retroposon/SVA"),
]
YOUNG_L1PA = {f"L1PA{n}_3end" for n in range(2, 9)}
# Dfam models L1HS and L1PA as 5' and 3' ends only: the youngest ORF2 model
# (L1P1, 96-99% identical to the ORF2 of L1HS insertions) covers the middle.
ORF2 = {"L1P1_orf2"}


def get(url):
    with urllib.request.urlopen(url, timeout=60) as response:
        return response.read().decode()


def families(prefix):
    data = json.loads(get(f"{API}?name={prefix}&clade=9606&clade_relatives=both&limit=200"))
    return [r for r in data["results"] if r["name"].startswith(prefix)]


def main(path):
    records = []
    for prefix, kind in GROUPS:
        for family in families(prefix):
            name = family["name"]
            if prefix == "L1PA" and name not in YOUNG_L1PA:
                continue
            if prefix == "L1P1" and name not in ORF2:
                continue
            fasta = get(f"{API}/{family['accession']}/sequence?format=fasta")
            header, *lines = fasta.strip().split("\n")
            accession = header[1:].split()[0]
            records.append((name, kind, accession, "".join(lines).upper()))
    records.sort(key=lambda r: ([g for g, _ in GROUPS].index(next(g for g, _ in GROUPS if r[0].startswith(g))), r[0]))
    with open(path, "w") as out:
        for name, kind, accession, seq in records:
            out.write(f">{name}#{kind} {accession}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i : i + 60] + "\n")
    print(f"{len(records)} sequences written to {path}", file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "mei_consensus.fa")
