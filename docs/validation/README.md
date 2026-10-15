# Validation against truth sets

Precision and recall of `mei-rs scan` on real exomes, measured against mobile
element insertions called from long-read assemblies, and comparison with
Scramble on the same exomes (issue #26).

## Truth set

HGSVC3 MEI calls against GRCh38 (Logsdon et al., *Complex genetic variation in
nearly complete human genomes*, Nature 2025): 12,642 non-reference insertions
(Alu 10,270, L1 1,604, SVA 764, other 4) called with L1ME-AID and PALMER on
phased assemblies of 65 samples, with a genotype per haplotype.

```bash
B=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Mobile_Elements/1.0
curl -O $B/MEI_Callset_GRCh38.ALL.20241211.csv.gz   # md5 0131ecb529c9225dab677aced8d72b0a
```

Exomes:

* seven 1000 Genomes exomes with HGSVC3 genotypes: HG00512, HG00513, HG00731,
  HG00733, NA19238, NA19239, NA19240 (GRCh38 alignments of the 1000 Genomes
  GRCh38 exome collection, 76 to 101 bp reads; 371 truth MEIs in the padded
  targets);
* HG002 (`NA24385` in the truth set): GIAB exome aligned to GRCh37 (hs37d5);
  the truth positions are lifted to GRCh37 with the UCSC `hg38ToHg19` chain
  (12,620 of 12,642 lifted; 49 truth MEIs in the padded targets). HG002 was
  not used to fit anything (the quality score of #34 was fitted on the seven
  1000 Genomes exomes), so it is a held-out sample.

Capture targets: the 1000 Genomes exome pull-down targets
(`20130108.exome.targets.bed`, GRCh37) for HG002, and the same targets lifted
to GRCh38 with CrossMap and the UCSC `hg19ToHg38` chain for the 1000 Genomes
exomes, keeping the targets mapped in one piece (99.4%).

## Method

The whole evaluation is scripted and versioned here; `MANIFEST.md` lists every
public input with its checksum.

```bash
docs/validation/fetch-data.sh data/       # truth, references, 7 CRAMs + HG002, targets (~70 GB)
cargo build --release
docs/validation/run.sh data/ out/         # scan, evaluate, refit
```

`run.sh` runs `mei-rs scan --keep-filtered --reference` on each exome (the
GRCh38 targets for the 1000 Genomes exomes, the GRCh37 targets and `hs37d5` for
HG002), then `evaluate.py` and `refit.py`. By hand, for one sample:

```bash
mei-rs scan --bam SAMPLE.exome.bam --bed exome_targets.hg38.bed \
  --reference GRCh38_full_analysis_set_plus_decoy_hla.fa --keep-filtered -o SAMPLE.tsv
python3 docs/validation/evaluate.py \
  --truth MEI_Callset_GRCh38.ALL.20241211.csv.gz --bed exome_targets.hg38.bed \
  --min-qual 0,2,5,10,15 \
  HG00512=HG00512.tsv HG00513=HG00513.tsv ...

# HG002 (GRCh37): lift the truth set with a UCSC chain file
python3 docs/validation/evaluate.py \
  --truth MEI_Callset_GRCh38.ALL.20241211.csv.gz --bed exome_targets.GRCh37.bed \
  --chain hg38ToHg19.over.chain.gz NA24385=HG002.tsv
```

`refit.py` refits the quality score of `genotype::score` (prior and likelihood
ratios of the seven pieces of evidence, one call added to each count) on the
labelled PASS calls: `evaluate.py --fp-out` and `--tp-out` give the labels. It
prints the values for `src/genotype/score.rs`, the AUC leaving one sample out
and in sample, the true and false calls kept at each QUAL, and the AUC on a
held-out sample (HG002). On the v0.4.0 calls it gives the constants of the
release (prior 0.501, AUC 0.928 and 0.787).

`evaluate.py` (Python standard library only; `--help` for the options):

* keeps the truth MEIs whose insertion point lies in a target padded by
  `--padding` (300 bp, as `scan`), after the liftover with `--chain`; a truth
  MEI belongs to a sample when one haplotype carries it; contig names are
  compared without a leading `chr`;
* reads the calls from a `mei-rs scan` TSV (`--keep-filtered`, for the recall
  of filtered candidates) or from a VCF of `mei-rs` or of another caller
  (`FILTER` `PASS` or `.` passes); calls outside the padded targets are left
  out (other callers scan the whole alignment file);
* matches a call to a truth MEI of the same sample when its position, left
  junction or right junction is within `--tolerance` (20 bp) of the truth
  position (both 1-based, insertion right after the position);
* precision: PASS calls matching a truth MEI of the sample, over PASS calls
  matching it or no MEI of any HGSVC3 sample (calls near an MEI of another
  sample, or of a sample with a missing haplotype, are counted apart);
* recall: truth MEIs of the sample matched by a PASS call (`recall_pass`) or
  by any candidate, filtered or not (`recall_any`), also by distance to the
  unpadded target and by family;
* genotype concordance of the true PASS calls, and precision and recall of the
  PASS calls above each `--min-qual` threshold;
* 95% bootstrap intervals of precision and recall (PASS and each `--min-qual`)
  resampling the samples (`--bootstrap`, 2000 resamples, fixed seed): the
  samples are the unit, since the calls of one sample share depth, library and
  ancestry; HG002 alone has no interval;
* `--fp-out`, `--tp-out` (true positives, with calls near an MEI of another
  sample or of a sample not genotyped there under their own status) and
  `--fn-out` list the PASS false positives, the PASS calls that are not, and
  the missed MEIs.

## Results

v0.4.0 (`develop` after #67), default options with `--reference`: seven
exomes 126 true / 40 false PASS (precision 0.759, recall 0.337, recall of any
candidate 0.437); HG002 (held out) 17 true / 24 false (0.415 / 0.347). Per
sample:

| Sample | Truth MEIs | PASS | True PASS | Precision | Recall (PASS) |
|---|---|---|---|---|---|
| HG00512 | 49 | 27 | 19 | 0.731 | 0.388 |
| HG00513 | 45 | 28 | 22 | 0.786 | 0.467 |
| HG00731 | 45 | 11 | 9 | 0.818 | 0.200 |
| HG00733 | 45 | 15 | 11 | 0.786 | 0.244 |
| NA19238 | 59 | 38 | 31 | 0.816 | 0.525 |
| NA19239 | 60 | 26 | 18 | 0.720 | 0.300 |
| NA19240 | 68 | 24 | 16 | 0.667 | 0.235 |

By confidence: `high` 13 true of 14, `medium` 35 of 38, `low` 78 of 117. Recall
(PASS) by distance to the target: inside 27 / 36, 1-100 bp 62 / 95, 101-300 bp
36 / 240; by family: Alu 106 / 325, L1 16 / 32, SVA 3 / 14. Genotypes of the
true PASS calls: 91 of 123 right (107 of 123 once a PASS call is no longer
genotyped `0/0`, #91; HG002 9 of 17 then 11 of 17). `--min-qual 2`: 125 / 15 (0.893; HG002
17 / 13), `--min-qual 5`: 114 / 6 (0.950; HG002 15 / 11). The quality score
was refit on these calls (#67): AUC 0.928 leaving one sample out, 0.787 on
HG002. NA12878 (GRCh38, 12 cores): 23 PASS, BAM 120 s and 1.4 GB; HG00731:
170 s and 2.5 GB. The numbers below are those of v0.3.0 and v0.2.0.

95% bootstrap intervals over the seven exomes (`evaluate.py`; rerun with the
scripts above, #83): precision 0.717 to 0.799, recall 0.256 to 0.423;
`--min-qual 2` precision 0.871 to 0.919, `--min-qual 5` precision 0.920 to
0.984 and recall 0.216 to 0.394. The intervals say how much the seven samples
vary, not how a new cohort would: HG002 alone (precision 0.415) lies far
below the interval.

v0.3.0 (`develop` after #70), default options with `--reference`: seven
exomes 121 true / 39 false PASS (precision 0.756, recall 0.326, recall of any
candidate 0.437; 377 truth MEIs in the padded targets of any sample, 371 in
these seven); HG002 (held out) 16 true / 23 false (0.410 / 0.327). Per sample:

| Sample | Truth MEIs | PASS | True PASS | Precision | Recall (PASS) |
|---|---|---|---|---|---|
| HG00512 | 49 | 28 | 20 | 0.741 | 0.408 |
| HG00513 | 45 | 27 | 20 | 0.741 | 0.444 |
| HG00731 | 45 | 11 | 9 | 0.818 | 0.200 |
| HG00733 | 45 | 15 | 11 | 0.786 | 0.244 |
| NA19238 | 59 | 36 | 30 | 0.833 | 0.508 |
| NA19239 | 60 | 23 | 16 | 0.727 | 0.267 |
| NA19240 | 68 | 23 | 15 | 0.652 | 0.221 |

By confidence: `high` 13 of 14, `medium` 34 of 37, `low` 74 of 112. Recall
(PASS) by distance to the target: inside 28 / 36, 1-100 bp 61 / 95, 101-300 bp
32 / 240; by family: Alu 103 / 325, L1 15 / 32, SVA 3 / 14. Genotypes of the
true PASS calls: 90 of 119 right. `--min-qual 2`: 119 / 14 (0.895; HG002
15 / 12), `--min-qual 5`: 103 / 5 (0.954; HG002 13 / 8). The numbers below
are those of v0.2.0.

v0.2.0 (`develop` after #56), default options with `--reference`.

| Sample | Truth MEIs | PASS | True PASS | Precision | Recall (PASS) | Recall (any) |
|---|---|---|---|---|---|---|
| HG00512 | 49 | 25 | 16 | 0.667 | 0.327 | 0.429 |
| HG00513 | 45 | 25 | 17 | 0.680 | 0.378 | 0.467 |
| HG00731 | 45 | 11 | 7 | 0.636 | 0.156 | 0.244 |
| HG00733 | 45 | 16 | 10 | 0.667 | 0.222 | 0.311 |
| NA19238 | 59 | 30 | 23 | 0.767 | 0.390 | 0.678 |
| NA19239 | 60 | 20 | 12 | 0.632 | 0.200 | 0.300 |
| NA19240 | 68 | 23 | 14 | 0.609 | 0.206 | 0.412 |
| **All seven** | **371** | **150** | **99** | **0.673** | **0.267** | **0.412** |
| HG002 (held out) | 49 | 34 | 13 | 0.382 | 0.265 | 0.408 |

(Precision is true PASS over true and false PASS; 3 PASS calls of the seven
exomes are near an MEI of another sample or of a missing haplotype.)

By confidence (seven exomes): `high` 12 true of 13, `medium` 33 of 38, `low`
54 of 99. HG002: `high` 1 of 2, `medium` 3 of 6, `low` 9 of 26.

Recall by distance of the truth MEI to the unpadded target (seven exomes):

| Distance | Any candidate | PASS |
|---|---|---|
| inside the target | 32 / 36 | 20 / 36 |
| 1-100 bp | 71 / 95 | 53 / 95 |
| 101-300 bp | 50 / 240 | 26 / 240 |

Recall (PASS) by family: Alu 82 / 325, L1 14 / 32, SVA 3 / 14.

Genotypes of the true PASS calls: 81 of 99 right (heterozygous called `0/1`
57 of 74, homozygous called `1/1` 24 of 25); HG002 8 of 13.

### Quality threshold

`--min-qual` is 0 (off) by default. PASS calls with `qual` at or above a
threshold:

| `--min-qual` | Seven exomes: true / false | Precision | Recall | HG002: true / false | Precision | Recall |
|---|---|---|---|---|---|---|
| 0 (default) | 99 / 48 | 0.673 | 0.267 | 13 / 21 | 0.382 | 0.265 |
| 2 | 97 / 18 | 0.843 | 0.261 | 13 / 11 | 0.542 | 0.265 |
| 5 | 82 / 9 | 0.901 | 0.221 | 10 / 7 | 0.588 | 0.204 |
| 10 | 48 / 4 | 0.923 | 0.129 | 6 / 6 | 0.500 | 0.122 |
| 15 | 19 / 0 | 1.000 | 0.051 | 1 / 1 | 0.500 | 0.020 |

The score was fitted on the seven exomes, so their rows are in sample; HG002
is not. On both, a threshold of 2 removes about half of the false positives
for at most two true calls.

### Remaining errors

* False positives are mostly `low` (no poly(A/T) tail): 44 of the 48 of the
  seven exomes, 17 of the 21 of HG002 (14 of these 17 are called `1/1`).
* Recall is limited by coverage in the padding: 20 of 36 truth MEIs inside
  the targets are called PASS, 26 of 240 at 101-300 bp from them. The
  missed MEIs within 100 bp of the targets are analysed in the v0.3.0 issues
  (#57, #58, #59): candidates left untyped or without enough junction clips,
  and mate-only candidates placed tens of bases off the insertion.
* LINE-1: Dfam models L1HS and L1PA as a 5' end and a 3' end, without the
  middle of ORF2, where the 5' junction of many 5'-truncated LINE-1 falls.
  The L1P1 ORF2 model of Dfam, added in #56, types them (identity 0.96-1.00 on
  the clips): LINE-1 recall 6 / 32 before, 14 / 32 after, for 3 more false
  positives (`qual` 0 to 4).

## Comparison with Scramble

[Scramble](https://github.com/GeneDx/scramble) 1.0.2 (GeneDx, commit
`14806b4`), the MEI caller designed for exomes, on the same BAM files:
`cluster_identifier` then `SCRAMble.R --eval-meis --ref`, default options,
PASS records of its VCF.

| | PASS | True | False | Precision | Recall (PASS) | HG002: true / false | HG002 recall |
|---|---|---|---|---|---|---|---|
| Scramble | 84 | 57 | 26 | 0.687 | 0.154 | 6 / 5 | 0.122 |
| `mei-rs` (default) | 150 | 99 | 48 | 0.673 | 0.267 | 13 / 21 | 0.265 |
| `mei-rs --min-qual 5` | 91 | 82 | 9 | 0.901 | 0.221 | 10 / 7 | 0.204 |

Of the 371 truth MEIs of the seven exomes, 57 are called PASS by both tools
and 42 by `mei-rs` only; every MEI called by Scramble is called by `mei-rs`.
Before #56, 8 of them (all LINE-1, at two sites whose 5' junction falls in
ORF2) were left untyped by `mei-rs`. Recall by family: Scramble Alu 45 / 325,
L1 10 / 32, SVA 2 / 14. Scramble does not genotype.

Run time and memory, both tools one after the other in the same Linux arm64
container (Docker on an Apple M5 Pro, 18 cores; alignment files on an
external USB SSD; `mei-rs` before #56, whose one more consensus sequence does
not change the time). `mei-rs` writes the VCF (PASS only); its wall time uses all
cores, Scramble is single-threaded (`cluster_identifier`, then R):

| Sample | `mei-rs` wall | CPU | Max RSS | Scramble wall | CPU | Max RSS |
|---|---|---|---|---|---|---|
| HG00512 | 38 s | 314 s | 0.6 GB | 139 s | 137 s | 3.6 GB |
| HG00513 | 44 s | 396 s | 0.7 GB | 129 s | 127 s | 3.5 GB |
| HG00731 | 159 s | 1,283 s | 1.7 GB | 849 s | 843 s | 4.1 GB |
| HG00733 | 60 s | 520 s | 1.2 GB | 59 s | 58 s | 3.5 GB |
| NA19238 | 87 s | 743 s | 1.0 GB | 136 s | 132 s | 3.6 GB |
| NA19239 | 30 s | 279 s | 0.5 GB | 55 s | 53 s | 3.5 GB |
| NA19240 | 52 s | 449 s | 0.8 GB | 74 s | 72 s | 3.5 GB |
| HG002 | 93 s | 655 s | 1.5 GB | 118 s | 115 s | 3.5 GB |

`mei-rs` is 1 to 5 times faster in wall time on 18 cores, but uses 1.5 to 9
times more CPU time; Scramble's memory is mostly the reference genome loaded
by R. Scramble's time grows with its clusters (HG00731: 261,677 clusters,
most of its time in R).

Scramble has no Linux arm64 package: it was built from source (Ubuntu 22.04,
htslib 1.13, R 4.1, Biostrings 2.62) without `rBLAST`, which only deletion
calling needs; `cluster_identifier` gives the same clusters as the bioconda
amd64 build on a test region.

## Comparison with MELT and xTea

Same seven exomes and truth set (371 MEIs in the padded targets), PASS calls
of each tool (`evaluate.py` reads their VCF; `mei-rs` is v0.4.0, Scramble is
the v0.2.0 run above):

| | PASS | True | False | Precision | Recall (PASS) | Alu | L1 | SVA |
|---|---|---|---|---|---|---|---|---|
| `mei-rs` v0.4.0 | 169 | 126 | 40 | 0.759 | 0.337 | 106 / 325 | 16 / 32 | 3 / 14 |
| MELT 2.2.2 | 110 | 96 | 11 | 0.897 | 0.259 | 79 / 325 | 14 / 32 | 3 / 14 |
| xTea 0.1.9, default | 27 | 23 | 4 | 0.852 | 0.062 | 21 / 325 | 2 / 32 | 0 / 14 |
| xTea 0.1.9, exome cutoffs | 43 | 33 | 10 | 0.767 | 0.089 | 30 / 325 | 3 / 32 | 0 / 14 |
| Scramble 1.0.2 | 84 | 57 | 26 | 0.687 | 0.154 | 45 / 325 | 10 / 32 | 2 / 14 |

HG002 (49 MEIs, GRCh37): MELT 5 true / 1 false (precision 0.833, recall
0.102), xTea default 2 / 0 (recall 0.041), xTea exome cutoffs 3 / 0 (0.061),
against 17 / 24 for `mei-rs` v0.4.0. For xTea on GRCh37 the `xtea_hg19` entry
point, its hg19 library and GENCODE v19 were used, with the `chr` prefix
removed from the annotation files (the hs37d5 contigs are `1`..`22`).

* MELT 2.2.2 (`Single -exome`, one run per element: ALU, LINE1, SVA; `-c 100`,
  `-r` the read length, human gene annotation of its distribution; the BAM is
  linked into the working directory because MELT writes next to it). MELT is
  the most precise caller on these exomes, `mei-rs` finds more of the MEIs.
* xTea 0.1.9 (`-f 5907 -y 7`, one pipeline per element, repeat library and
  GENCODE v46 added). Its automatic cutoffs (clip 3, discordant 6 from a
  coverage of about 40x) leave 3 to 5 calls per exome. "Exome cutoffs" are
  `-U --nclip 2 --cr 1 --nd 2 --nfclip 2 --nfdisc 2`: more candidates survive
  the first steps (573 against 137 on HG00512, Alu) but its later filters
  remove almost all of them, so recall barely changes.
* Run time and memory are not compared: both tools only run as linux/amd64
  images, emulated on this arm64 machine (and while other jobs were running),
  so their wall time says nothing about the tools.

## Genotypes (#91)

Of the 123 true PASS calls with reads at the site (v0.4.0), 14 heterozygous
insertions were called `0/0` and 15 `1/1`. A candidate that passed is an
insertion: `0/0` is no longer a choice for it, and the heterozygous genotype
has a prior of 0.77 (the share of heterozygous genotypes among the HGSVC3
insertions of the seven exomes). Result on the same calls (unchanged): 107 of
123 right on the seven exomes, 11 of 17 on HG002 (9 of 17 before). Six
heterozygous insertions of HG002 and 13 of the seven exomes are still called
`1/1`: their sites have no reference read (a probe edge, where the reads
upstream are missing; the alternative reads start at the target start), which
no count of the alternative allele corrects. Fitting the heterozygous fraction
by candidate type, or a beta-binomial model, gave the same 105 of 123 in a
leave-one-sample-out check (before the prior): the remaining errors are in the
counts, not in the model.

## VCF with bcftools, tabix and VEP (#85)

* Fixtures: `tests/vcf_tools.sh` (CI job `vcf`) writes the VCF of every
  fixture (BAM and CRAM, plain and bgzipped) and reads it with bcftools 1.21,
  tabix and `vcf-validator` (vcftools): no warning, the same records plain and
  bgzipped, the REF alleles equal to the reference (`bcftools norm
  --check-ref e`), and the `.tbi` we write returns the same records as the one
  tabix builds, on every record and contig. The validator reports an
  undeclared INFO field, so it does check the header.
* NA12878 (23 `PASS` records, `--format vcf`): Ensembl VEP 116.3 (database
  mode, GRCh38) reads all 23 records without warning and annotates them with
  the alleles `Alu_insertion` and `LINE1_insertion` (the symbolic ALTs
  `<INS:ME:ALU>` and `<INS:ME:LINE1>` are recognised): 559 lines of
  transcript consequences, mostly intron variants. VEP spans each insertion
  over `SVLEN` reference bases (`chr1:28571787-28572087` for an insertion of
  301 bp, with or without an `END` INFO field), so a long element is annotated
  against every transcript of a window as long as the element, not at one
  point. AnnotSV was not tried.

## History

Same seven exomes, PASS calls:

| Version | PASS | True | Precision | Recall (PASS) | `high` true / all |
|---|---|---|---|---|---|
| after #24 (v0.1.1 baseline) | 2,392 | 60 | 0.025 | 0.162 | 5 / 84 |
| #42 `tandem_repeat` (v0.1.1) | 2,002 | 58 | 0.029 | 0.156 | 5 / 24 |
| #30 typing | 216 | 72 | 0.338 | 0.194 | 5 / 6 |
| #31 MEI mates | 110 | 70 | 0.654 | 0.189 | 5 / 6 |
| #32 short clips | 140 | 91 | 0.664 | 0.245 | 9 / 10 |
| #33 assembly, #34 score, #35 VCF | 139 | 91 | 0.669 | 0.245 | 12 / 13 |
| #56 LINE-1 ORF2 (v0.2.0) | 150 | 99 | 0.673 | 0.267 | 12 / 13 |
| #57 typing by the mates | 163 | 111 | 0.694 | 0.299 | 12 / 13 |
| #58, #59 insertions without junction clips | 170 | 118 | 0.707 | 0.318 | 12 / 13 |
| #69 mate checks | 159 | 117 | 0.750 | 0.315 | 13 / 14 |
| #70 single typed clips (v0.3.0) | 163 | 121 | 0.756 | 0.326 | 13 / 14 |
| #71 tail mates and reference copies, #66, #67 (v0.4.0) | 169 | 126 | 0.759 | 0.337 | 13 / 14 |

#57 types untyped candidates with a poly(A/T) tail clip by their mates
(junction clips that are the tail alone): 12 more true calls, 1 more false
(`--min-qual 2`: 109 true / 19 false, precision 0.852); within 100 bp of the
targets, 81 of 131 truth MEIs are PASS (73 before). HG002: 16 true / 25
false (13 / 21 before). Without the identity (96%) and tailed-clip (2)
conditions, the mates typed 35 more false candidates for 13 true: most
false ones had one tailed clip and mates matching older reference copies.

#58 places the candidates without junction clips at the inner edge of their
mate signals and, for those with MEI mates, among the reads soft-clipped by
3 bases or more (of 15 truth MEIs whose candidate was 43-431 bp off, 8 are
then within 20 bp and 4 within the interval). #59 checks their mates and
passes them when young MEI mates (98% identity, at least 60% and 4 of the
checked mates) type them and reads are clipped at their position: 7 more
true calls, no false one (`--min-qual 2`: 116 true / 19 false, precision
0.859); within 100 bp of the targets, 85 of 131 truth MEIs are PASS. HG002
unchanged. Looser thresholds (96% identity, without the clips) added 11
true and 14 false candidates.

#69 changes how mates are checked: mapped mates are read over the span of
their alignment (`MC` tag) and typed without their A/T runs, and mapped
mates found in the reference within 600 bp of their candidate are not MEI
mates (next to a reference copy of an element, the mates of the reads in
that copy map with a low MAPQ to other copies). 1 fewer true call, 10 fewer
false (`--min-qual 2`: 115 true / 14 false, precision 0.891; `--min-qual 5`:
100 / 5, 0.952); within 100 bp of the targets, 85 of 131 truth MEIs are PASS.
HG002: 16 true / 23 false (`--min-qual 2`: 15 / 12). It was meant to raise
recall: the mates of several missed MEIs (`non_mei_mates`) carry 30 to 50 A
of the poly(A) tail and 30 to 45 bases of the element end, too few to be
typed. Not counting untyped mates with a long A/T run raised recall (124
true) but added false calls next to reference Alus (53 false); set aside in
#71.

#70 checks the mates of candidates typed by a single clip without enough
junction clips and lets them pass as mate-only insertions (#59), lets one
typed clip without a tail pass with mostly young MEI mates, and lifts
`low_ratio` with 5 or more typed junction clips: 4 more true calls, no false
one (`--min-qual 2`: 119 true / 14 false, precision 0.895); within 100 bp of
the targets, 89 of 131 truth MEIs are PASS. HG002 unchanged.

Before typing, most false positives were poly(A) length variants of reference
Alus (`high`) and background unmapped mates of the deep samples (`medium`):
#42 and #31 address them. The last two rows leave out one PASS call outside the
padded targets, counted before the evaluation restricted the calls to them.
