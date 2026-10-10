# Validation against truth sets

Precision and recall of `mei-rs scan` on real exomes, measured against mobile
element insertions called from long-read assemblies (issue #26).

## Truth set

HGSVC3 MEI calls against GRCh38 (Logsdon et al., *Complex genetic variation in
nearly complete human genomes*, Nature 2025): 12,642 non-reference insertions
(Alu 10,270, L1 1,604, SVA 764, other 4) called with L1ME-AID and PALMER on
phased assemblies of 65 samples, with a genotype per haplotype.

```bash
B=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Mobile_Elements/1.0
curl -O $B/MEI_Callset_GRCh38.ALL.20241211.csv.gz   # md5 0131ecb529c9225dab677aced8d72b0a
```

Seven 1000 Genomes exomes have HGSVC3 genotypes: HG00512, HG00513, HG00731,
HG00733, NA19238, NA19239, NA19240 (GRCh38 alignments of the 1000 Genomes
GRCh38 exome collection, 76 to 101 bp reads). HG002 (`NA24385`) is in the
truth set too, but its GIAB exome is aligned to GRCh37 (not evaluated yet).

Capture targets: the 1000 Genomes exome pull-down targets
(`20130108.exome.targets.bed`, GRCh37), lifted to GRCh38 with CrossMap and the
UCSC `hg19ToHg38` chain, keeping the targets mapped in one piece (99.4%).

## Method

```bash
mei-rs scan --bam SAMPLE.exome.bam --bed exome_targets.hg38.bed \
  --keep-filtered -o SAMPLE.tsv
python3 docs/validation/evaluate.py \
  --truth MEI_Callset_GRCh38.ALL.20241211.csv.gz --bed exome_targets.hg38.bed \
  HG00512=HG00512.tsv HG00513=HG00513.tsv ...
```

`evaluate.py` (Python standard library only; `--help` for the options):

* keeps the truth MEIs whose insertion point lies in a target padded by
  `--padding` (300 bp, as `scan`); a truth MEI belongs to a sample when one
  haplotype carries it;
* matches a candidate to a truth MEI of the same sample when its `position`,
  `left_junction` or `right_junction` is within `--tolerance` (20 bp) of the
  truth position (both 1-based, insertion right after the position);
* precision: PASS candidates matching a truth MEI of the sample, over PASS
  candidates matching it or no MEI of any HGSVC3 sample (candidates near an
  MEI of another sample, or of a sample with a missing haplotype, are counted
  apart);
* recall: truth MEIs of the sample matched by a PASS candidate (`recall_pass`)
  or by any candidate, filtered or not (`recall_any`), also by distance to the
  unpadded target and by family;
* `--fp-out` and `--fn-out` list the PASS false positives and the missed MEIs.

## Results

`develop` after #24 (unreleased v0.1.1), default options, 12 cores.

| Sample | Truth MEIs | PASS | True PASS | Precision | Recall (PASS) | Recall (any) | Time, max RSS |
|---|---|---|---|---|---|---|---|
| HG00512 | 49 | 80 | 9 | 0.114 | 0.184 | 0.388 | 25 s, 0.7 GB |
| HG00513 | 45 | 94 | 10 | 0.106 | 0.222 | 0.444 | 30 s, 0.8 GB |
| HG00731 | 45 | 527 | 1 | 0.002 | 0.022 | 0.222 | 106 s, 2.3 GB |
| HG00733 | 45 | 382 | 2 | 0.005 | 0.044 | 0.244 | 41 s, 1.4 GB |
| NA19238 | 59 | 1,054 | 24 | 0.023 | 0.407 | 0.610 | 61 s, 1.2 GB |
| NA19239 | 60 | 71 | 7 | 0.100 | 0.117 | 0.300 | 25 s, 0.7 GB |
| NA19240 | 68 | 184 | 7 | 0.038 | 0.103 | 0.368 | 41 s, 1.0 GB |
| **All** | **371** | **2,392** | **60** | **0.025** | **0.162** | **0.375** | |

By confidence: `high` 5 true of 84, `medium` 54 of 2,259, `low` 1 of 49.

Recall by distance of the truth MEI to the unpadded target:

| Distance | Any candidate | PASS |
|---|---|---|
| inside the target | 30 / 36 | 16 / 36 |
| 1-100 bp | 67 / 95 | 24 / 95 |
| 101-300 bp | 42 / 240 | 20 / 240 |

Recall (PASS) by family: Alu 50 / 325, L1 6 / 32, SVA 4 / 14.

### What the false positives are

* **`high`: poly(A) length variants of reference Alus.** 79 of the 84 `high`
  candidates are false; the soft-clips of 53 of them are explained by the
  nearby reference (see below). They sit at the 3' end of an Alu of the
  reference (`...GACTCCATCTCAAAAAAAAA`): the reads carry a longer A tract (or
  `(CAAAA)n` tract) than the reference, so the aligner clips them on both sides
  of the tract, and the overlap of the two junctions looks like a TSD. The same
  sites recur across samples (e.g. chr2:201350686, a top `high` of NA12878).
* **`medium` in deep samples: background mate signals.** HG00731, HG00733 and
  NA19238 (local depth up to ~1,000x) carry most of the false positives: the
  background of reads with an unmapped mate reaches `--min-support`, and a
  single soft-clip with a poly(A/T) tail passes `no_tail` (461 of the 527 PASS
  of HG00731 have one tailed clip).
* True insertions are mostly supported by mates mapped elsewhere
  (`discordant_mates`, typically in another copy of the element) with a
  `signal_ratio` of at least 0.2 (133 of the 139 matched candidates), and
  often fail `no_tail`: with 76 bp reads the tailed junction is frequently
  not covered by clips of 20 bp or more.

### Rules tried (not applied)

Post-filters on the same outputs:

| PASS rule | PASS | True | Precision | Recall (PASS) |
|---|---|---|---|---|
| current (#22) | 2,392 | 60 | 0.025 | 0.162 |
| + soft-clips explained by the nearby reference (homopolymer-compressed clip found within 60 bp) | 2,115 | 55 | 0.026 | 0.148 |
| tailed clips ≥ 2 and `signal_ratio` ≥ 0.3 | 219 | 38 | 0.175 | 0.102 |
| both | 138 | 34 | 0.250 | 0.092 |

A grid of read-count thresholds (discordant mates, ratio, junction clips,
tailed clips) does not exceed about 0.3 precision at any useful recall. Read
counts alone do not separate insertions from capture and alignment noise;
the sequence content of the clips and mates must be checked against MEI
consensus sequences (#30, #31).

## Not done yet

* HG002: GIAB exome on GRCh37, truth set on GRCh38 (lift one of them over).
* Comparison with Scramble, MELT and xTea on the same exomes (precision,
  recall, time and memory).
