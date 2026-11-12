# Thresholds of the pipeline

Every threshold that decides a call or a score, with its value and where it
comes from. The HGSVC3 truth set and the seven 1000 Genomes exomes used to fit
or choose them are the ones of [docs/validation](validation/README.md); HG002
was never used. A value marked **fitted** or **chosen** was set on those seven
exomes (from #90 on, on the 18 tuning exomes of `validation/split.tsv`), so
the precision and recall measured on them are optimistic for these values.
The 17 test exomes were never fitted on but were looked at (false-positive
analysis, model and threshold comparisons, #122): they are a development set;
HG002 is the only fully external check.

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
| `--profile`, `--min-qual` | `exome`: 10 (`discovery`: 0); `exome` requires `--reference` (#120) | chosen (#90) on the 18 tuning exomes of `docs/validation/split.tsv`: lowest quality with a precision of 0.95 there, an internal target to choose the value, not a precision to expect (test exomes: 0.908, 95% interval over families and loci 0.809 to 0.977); 0 until v0.4.1 |

## Internal constants

| Name | Value | Where | Origin |
|---|---|---|---|
| `MIN_IDENTITY`, `MIN_COVERAGE` (typing) | 0.9, 0.7 | `typer` | design (#30) |
| Seed `K`, `BAND`, `ALIGNED_STRANDS` (typing) | 11, 20 bp, 8 | `typer` | design |
| Smith-Waterman scores | +2 / -3, gap 5 + 2k | `typer::align` | design |
| `MIN_TYPED_LENGTH` | 15 bp | `typer` | design; false typing rate measured below (#133) |
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
| `HOM_MIN_ALT_FRACTION` | 0.2 | `genotype` | fitted (#91) on the 18 tuning exomes: without reference read at the site, a `1/1` call needs this fraction of alternative reads among the reads at the site (best of 0.15 to 0.3) |
| `HET_FRACTION` | 0.3 | `genotype` | fitted (#34; median alternative fraction of the true heterozygous calls 0.28; 76 of 91 genotypes right, 68 with 0.5) |
| `PRIOR` and the 14 log10 likelihood ratios of the score | see `genotype::score` | `genotype::score` | fitted (#34, refit in #67 on 129 true and 40 false `PASS` calls) |
| `MIN_SEX_READS`, `MALE_MAX_X_RATIO` | 200 reads, 0.95 | `genotype::ploidy` | chosen (#127: reads per targeted base of the X against the autosomes, 0.61 to 0.84 in the 22 male and 1.12 to 1.63 in the 14 female 1000 Genomes exomes of the evaluation) |
| Score cutoffs: identity 0.98, distinct starts 0.7, MEI mates 0.6, signal ratio 0.2 | | `genotype::score` | chosen (#34) |

The quality score is a naive Bayes model on seven binary pieces of evidence:
its 15 parameters (prior and two ratios per piece) are the most exposed to the
seven exomes; see the intervals and the HG002 result in
[docs/validation](validation/README.md).

## False typing of short clips

A short clip (15 to 19 bp, `--min-short-clip-len`) is kept when it has a
poly(A/T) tail or is typed to an MEI consensus, and counts towards
`no_junction` and `untyped` only when the candidate has an MEI mate (#32).
How often does a clip that is not an insertion get typed? `examples/clip_typing.rs`
(#133) takes the soft-clips of 40 bases or more without a tail from a scan
(`--profile discovery --keep-filtered --format json`) whose origin is known:

* **local**: found in the reference within 1 kb of the junction, on either
  strand and without gaps, with at most 10% of the bases different (a local
  rearrangement or an alignment artefact), and the matched reference not
  itself typed to an MEI (clips from a reference copy are left out: other
  filters handle them, #42, #71);
* **adapter**: holding the start of the Illumina adapter (read-through of a
  short fragment).

It types the bases next to the junction (what a shorter clip of that read
would hold) at each length, with the default consensus set, and random
sequences of the same lengths for comparison. Two 1000 Genomes exomes:
NA12878 (76 bp reads; 242 local and 745 adapter clips) and NA19331 (100 bp
reads; 2,141 local and 104 adapter clips).

| Length (bp) | Local, typed | Rate | Adapter, typed | Rate | Random (20,000) |
|---|---|---|---|---|---|
| 15 | 76 / 2,383 | 0.032 | 7 / 849 | 0.008 | 0.020 |
| 17 | 42 / 2,383 | 0.018 | 4 / 849 | 0.005 | 0.009 |
| 19 | 26 / 2,383 | 0.011 | 4 / 849 | 0.005 | 0.003 |
| 20 | 30 / 2,383 | 0.013 | 4 / 849 | 0.005 | 0.004 |
| 25 | 11 / 2,383 | 0.005 | 3 / 849 | 0.004 | 0.0001 |
| 30 | 10 / 2,383 | 0.004 | 0 / 849 | 0 | 0.0001 |
| 35, 40 | 1 / 2,383 | 0.0004 | 0 / 849 | 0 | 0 |

Real clips are typed more often than random sequences: about 3% at 15 bp
and 1% at 19-20 bp for local clips, against 2% and 0.3-0.4% at random. A likely
reason (not checked clip by clip): genomic sequence is closer to old repeats
than uniform bases are, and the piece next to the junction can hold a
fragment of an element too short or diverged for the whole matched reference
to be typed. The rate
falls below 0.5% only from 25 bp. A single typed short clip is therefore weak
evidence, which is why short clips count only with an MEI mate.
