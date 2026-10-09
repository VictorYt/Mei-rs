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
| `reference.fa`, `.fai` | Two random contigs: `chr1` (20 kb) and `chr2` (10 kb). `chr2:7001-7281` holds a reference copy of the Alu. |
| `targets.bed` | Three capture targets (see below). |
| `positive.bam`, `.bai`, `.cram`, `.crai` | Heterozygous Alu insertion, probe-edge artefact, PCR duplicates. |
| `negative.bam`, `.bai`, `.cram`, `.crai` | Same targets without the insertion; same probe-edge artefact. |
| `truth.tsv` | The simulated insertion, machine-readable. |

BAM and CRAM hold the same records (samtools adds `MD`/`NM` when decoding the
CRAM). The CRAM needs `reference.fa`.

## Targets (`targets.bed`, 0-based half-open)

| Name | Region (1-based) | Content in `positive` |
|---|---|---|
| `target_insertion` | `chr1:9851-10150` | Alu insertion |
| `target_probe_edge` | `chr1:3001-3200` | Probe-edge artefact |
| `target_negative` | `chr2:4001-4300` | Nothing |

## Simulated insertion (`positive`)

An Alu (281 bp, approximate AluY consensus, for test purposes only) followed by
a 30 bp poly(A) tail, inserted on the plus strand on `chr1`, heterozygous
(`0/1`). The 12 bp target site `chr1:10001-10012` (`CGCGGTCCCCCG`) is
duplicated (TSD):

```
alt = chr1[1..10012] + Alu + A(30) + chr1[10001..]
```

All positions are 1-based:

- **Left junction:** reads anchored on the left end at `chr1:10012` and are
  soft-clipped on their right side; the clipped part starts with the Alu.
- **Right junction:** reads anchored on the right start at `chr1:10001` and are
  soft-clipped on their left side; the clipped part ends with the poly(A) tail.
- **Discordant pairs:** reads entirely inside the Alu are placed on the `chr2`
  copy with MAPQ 0; other reads inside the insertion (poly(A), short anchors)
  are unmapped and placed at their mate. Their mates are anchored in
  `target_insertion`.
- **PCR duplicates:** 3 junction fragments are copied and flagged `0x400`.

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
