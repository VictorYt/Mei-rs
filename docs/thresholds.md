# Thresholds of the pipeline

Every threshold that decides a call or a score, with its value and where it
comes from. The HGSVC3 truth set and the seven 1000 Genomes exomes used to fit
or choose them are the ones of [docs/validation](validation/README.md); HG002
was never used. A value marked **fitted** or **chosen** was set on those seven
exomes, so the precision and recall measured on them are optimistic for these
values (HG002, held out, is the only independent check; independent exomes
are planned for the 1.0 publication).

Origins:

* **design**: set from the data type or the molecular mechanism, before any
  truth set existed (v0.1.0 defaults) or without comparing values;
* **chosen**: a few values compared on the seven exomes, the best kept (the
  issue says what was compared);
* **fitted**: estimated from the labelled `PASS` calls of the seven exomes
  (`docs/validation/refit.py` regenerates them).

## Options (defaults)

| Name | Value | Origin |
|---|---|---|
| `--padding` | 300 bp | design (fragment length minus read length on exomes) |
| `--min-mapq` | 20 | design |
| `--min-clip-len`, `--min-short-clip-len` | 20, 15 bp | design; short clips #32 (chosen on the seven exomes) |
| `--min-clip-quality` | 20 | design; #25 measured the filter on older data |
| `--polya-min-len`, `--polya-min-frac` | 10 bp, 0.8 | design |
| `--window` | 100 bp | design |
| `--min-support` | 3 reads | design |
| `--min-clips` | 2 | design (#22); revisited once truth sets existed: left |
| `--min-typed-clips` | 2 | chosen (#30) |
| `--min-mei-mate-fraction` | 0.3 | chosen (#31) |
| `--min-signal-ratio` | 0.05 | design (v0.1.0) |
| `--probe-edge-tolerance` | 2 bp | design (#29) |
| `--min-qual` | 0 (off) | decided 2026-10-11 to catch every case; see the README for what 2 and 5 give |

## Internal constants

| Name | Value | Where | Origin |
|---|---|---|---|
| `MIN_IDENTITY`, `MIN_COVERAGE` (typing) | 0.9, 0.7 | `typer` | design (#30) |
| Seed `K`, `BAND`, `ALIGNED_STRANDS` (typing) | 11, 20 bp, 8 | `typer` | design |
| Smith-Waterman scores | +2 / -3, gap 5 + 2k | `typer::align` | design |
| `MIN_TYPED_LENGTH` | 15 bp | `typer` | design |
| `JUNCTION_TOLERANCE` | 2 bp | `cluster` | design |
| `MATE_SPAN`, `MATE_OVERLAP` | 300, 20 bp | `cluster` | design |
| `MIN_MATE_IDENTITY` | 0.96 | `cluster::filters` | chosen (#57; true candidates typed by mates at 0.96 or more, most false ones at 0.90-0.95) |
| `MATE_ONLY_MIN_IDENTITY`, `_MIN_FRACTION`, `_MIN_MEI_MATES` | 0.98, 0.6, 4 | `cluster::filters` | fitted on 7 true mate-only `PASS` calls (#59, #70); revisited in #82 |
| `RATIO_EXEMPT_TYPED_CLIPS` | 5 | `cluster::filters` | chosen (#70; one true site at depth 1,570, no false candidate with 5 or more) |
| `MATE_ONLY_MIN_MATES`, `MATE_ONLY_MIN_RATIO` | 3, 0.2 | `typer::mates` | chosen (59 of the 61 true candidates above the ratio, a seventh of the others) |
| `MAX_CHECKED_MATES` | 40 | `typer::mates` | design (cost) |
| `LOCAL_WINDOW`, `LOCAL_K`, `LOCAL_MIN_SHARED` | 600 bp, 20, 0.5 | `typer::mates` | design (#69) |
| `MIN_MATE_RUN`, `TAIL_MATE_RUN` | 10, 20 bases | `typer::mates` | chosen (#69, #71) |
| `COPY_WINDOW` | 600 bp | `typer::mates` | chosen (#71: 300 bp gave 128 true / 42 false, 600 bp 126 / 40) |
| `tandem_repeat` constants (`WINDOW` 250, `ANCHOR_LENGTH` 8, `ANCHOR_SLACK` 30, `MIN_EXPLAINED` 2, `MIN_HOMOPOLYMER` 10 ...) | see `cluster::tracts` | `cluster::tracts` | chosen (#42; variants compared on the baseline) |
| Assembly `K`, `MIN_COUNT`, `MAX_CONTIG`, `MAX_TSD` | 21, 2, 600, 50 | `assembly` | design (#33) |
| `MIN_INSERT`, `MIN_SEED_GAP`, `MAX_SEED_GAP`, `SNP_CHECK` | 10, 3, 40, 8 | `assembly` | design (#33, #66) |
| `ERROR` (genotype) | 0.02 | `genotype` | design |
| `HET_FRACTION` | 0.3 | `genotype` | fitted (#34; median alternative fraction of the true heterozygous calls 0.28; 76 of 91 genotypes right, 68 with 0.5) |
| `PRIOR` and the 14 log10 likelihood ratios of the score | see `genotype::score` | `genotype::score` | fitted (#34, refit in #67 on 129 true and 40 false `PASS` calls) |
| Score cutoffs: identity 0.98, distinct starts 0.7, MEI mates 0.6, signal ratio 0.2 | | `genotype::score` | chosen (#34) |

The quality score is a naive Bayes model on seven binary pieces of evidence:
its 15 parameters (prior and two ratios per piece) are the most exposed to the
seven exomes; see the intervals and the HG002 result in
[docs/validation](validation/README.md).
