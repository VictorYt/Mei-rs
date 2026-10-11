# Synthetic test fixtures

Small, deterministic data for the unit and integration tests, simulated by
`simulate.py` and converted by `generate.sh`:

```bash
tests/data/generate.sh    # requires python3, samtools and perl
```

The output is identical from one run to the next for a given samtools version
(generated with samtools 1.24). Do not edit the files by hand: change the
simulation and regenerate.

## Files

| File | Content |
|---|---|
| `reference.fa`, `.fai` | Two random contigs: `chr1` (20 kb) and `chr2` (10 kb). `chr2:7001-7281` holds a reference copy of the Alu; `chr2:1961-2012` the last 40 bp of the Alu and a 12 bp poly(A) tract. |
| `targets.bed` | Five capture targets (see below). |
| `targets_tiled.bed` | The same, with `target_insertion` tiled with two overlapping probes (see below). |
| `targets_tract.bed` | One target around a reference poly(A) tract on `chr2` (see below). |
| `positive.bam`, `.bai`, `.cram`, `.crai` | Heterozygous Alu, LINE-1 and SVA insertions, probe-edge artefact, PCR duplicates; expanded reference poly(A) tract; chimeric junction. |
| `negative.bam`, `.bai`, `.cram`, `.crai` | Same targets without the insertions; same probe-edge artefact. |
| `truth.tsv` | The simulated insertions, machine-readable. |

BAM and CRAM hold the same records (samtools adds `MD`/`NM` when decoding the
CRAM). The CRAM needs `reference.fa`.

## Targets (`targets.bed`, 0-based half-open)

| Name | Region (1-based) | Content in `positive` |
|---|---|---|
| `target_insertion` | `chr1:9851-10150` | Alu insertion |
| `target_probe_edge` | `chr1:3001-3200` | Probe-edge artefact |
| `target_negative` | `chr2:4001-4300` | Chimeric junction (not an insertion) |
| `target_line1` | `chr1:14851-15150` | LINE-1 insertion |
| `target_sva` | `chr1:5851-6150` | SVA insertion |

## Simulated insertions (`positive`)

Three heterozygous (`0/1`) insertions on `chr1`, each followed by a 30 bp
poly(A) tail and flanked by a target site duplication (TSD). The elements are
Dfam consensus sequences (`data/consensus/mei_consensus.fa`, CC0) without
their own poly(A):

| Element | Strand | TSD (1-based) | Inserted |
|---|---|---|---|
| AluY (281 bp) | `+` | `chr1:10001-10012` (`CGCGGTCCCCCG`) | AluY + A(30) |
| L1HS_3end, last 500 bp (5'-truncated LINE-1) | `-` | `chr1:15001-15015` | reverse complement of L1 + A(30) |
| SVA_F, last 600 bp | `+` | `chr1:6001-6010` | SVA + A(30) |

```
alt = chr1[1..tsd_end] + inserted + chr1[tsd_start..]
```

For the Alu (all positions 1-based; the others are alike):

- **Left junction:** reads anchored on the left end at `chr1:10012` and are
  soft-clipped on their right side; the clipped part starts with the Alu.
- **Right junction:** reads anchored on the right start at `chr1:10001` and are
  soft-clipped on their left side; the clipped part ends with the poly(A) tail.
  For the minus-strand LINE-1, the poly(T) is on the left junction instead.
- **Discordant pairs:** reads entirely inside the Alu are placed on the `chr2`
  copy with MAPQ 0; other reads inside an insertion (poly(A), short anchors,
  the LINE-1 and the SVA, which have no reference copy) are unmapped and
  placed at their mate.
- **PCR duplicates:** 3 junction fragments of the Alu are copied and flagged
  `0x400`.

Reads are paired-end 2x100 bp, fragments of 300 ± 30 bp, about 100x over the
targets (50x per haplotype), with 0.2% substitution errors and base qualities
Q30-Q40. Reads are aligned from the simulation truth: a read crossing a
junction is soft-clipped, however short the clip, as long as at least 20 bp are
anchored.

## Probe-edge artefact (both samples)

10 reads start exactly at the start of `target_probe_edge` (`chr1:3001`) with a
25 bp random soft-clip on their left side (`25S75M`), and no poly(A). Every
soft-clip of that window stops at the probe edge: the probe-edge filter must
flag it.

## Insertion at probe edges (`targets_tiled.bed`)

`target_insertion` is replaced by two overlapping probes,
`target_insertion_left` (`chr1:9851-10012`) and `target_insertion_right`
(`chr1:10001-10150`): the first ends at the left junction, the second starts
at the right one. Every soft-clip of the insertion stops at a probe edge, but
its poly(A) tail, typed clips and TSD show a real insertion: it must stay
`PASS`, while the probe-edge artefact stays flagged.

## Expanded reference poly(A) tract (`positive`, `targets_tract.bed`)

`chr2:1961-2000` holds the last 40 bp of the Alu, followed by a 12 bp poly(A)
tract (`chr2:2001-2012`). One haplotype of `positive` carries 14 more A
(`chr2[1..2012] + A(14) + chr2[2013..]`), about 50x, simulated after the other
reads with their own random generator (the other reads do not change). A read
crossing the expansion is anchored on its longer side, as an aligner does:

- anchored upstream, it matches through the reference tract and is clipped
  after `chr2:2012` (the extra A, then the downstream flank);
- anchored downstream, it matches from `chr2:2001` and is clipped before it
  (the upstream flank, then the extra A).

The clips have a poly(A) tail and the two junctions overlap by 12 bp, like a
TSD, and those before the tract hold the end of the Alu (typed): without the
reference this is a `PASS` `high` candidate. With
`--reference`, the far ends of the clips are the reference flanks of the
tract: `tandem_repeat`.

## Chimeric junction (`positive`, `target_negative`)

Simulated after the other reads, with their own random generator:

- 3 reads anchored on `chr2` up to 4150 and soft-clipped on their right side
  with the first 25 bases of the Alu (`75M25S`);
- 6 fragments anchored just before `chr2:4150` whose mate is unmapped and
  random (chimeric or adapter sequence).

The junction is typed (Alu clips) and has enough clips for its mates to be
checked against the MEI consensus sequences, and none is a mobile element:
`non_mei_mates`, 6 checked mates, 0 MEI mates. A chimeric read with a piece
of Alu looks like this.
