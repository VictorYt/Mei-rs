# Changelog

All notable changes to Mei-rs are documented in this file. The format is
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
project follows [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.4.1] - 2026-10-15

Reproducibility and distribution: no change in the calls (HGSVC3 exomes 126
true / 40 false `PASS`, HG002 17 / 24, NA12878 23 `PASS`).

### Changed

- README (#86): the speed claim is the measured one (faster in wall time,
  1.5 to 9 times more CPU time than Scramble), the positioning against MELT,
  Scramble and xTea and the held-out precision (HG002) are stated, a usage
  guide gives the precision and recall of `--min-qual` 0, 2, 5 and 10, and the
  scope (research tool) and the AGPL implications are spelled out.

### Added

- Tests and benchmarks (#87): property tests (`tests/properties.rs`, proptest)
  of the BED parser, the CIGAR clip lengths, the poly(A/T) tail detection,
  k-mers, the local alignment and the typing; `cargo fuzz` targets for the same
  routines (`fuzz/`, `fuzz/run.sh`); `cargo bench` (criterion) of extraction,
  clustering and typing on the fixtures and of the routines they use; a
  fixture of a mate-only insertion (`mateonly`: short clips at the breakpoint,
  mates in the element) and its test. `reader::signals::clip_lengths` and
  `poly_tail` are public for them.

- Release binaries and container image (#84): a `Release` workflow, run when a
  release is published, builds Linux (x86_64, aarch64; glibc 2.17 or later,
  linked with cargo-zigbuild) and macOS (arm64, x86_64) binaries, attaches
  them to the release with their SHA-256 checksums, and publishes the image
  (`Dockerfile`) to `ghcr.io/victoryt/mei-rs` for linux/amd64 and linux/arm64.
  musl builds were tried and rejected: the same exome region took 40 s
  against 1 to 2 s, because of musl's allocator.

- VCF checks (#85): `tests/vcf_tools.sh` (CI job `vcf`) reads the VCF of every
  fixture with bcftools, tabix and the VCF validator; VEP reads the NA12878
  VCF (23 records) without warning. No change to the VCF was needed.

- `docs/thresholds.md` (#86): every threshold of the pipeline with its value
  and how it was set (design, chosen or fitted on the seven HGSVC3 exomes).
- Versioned evaluation (#83): `docs/validation/` has `fetch-data.sh`,
  `build-targets.sh` (exome targets lifted to GRCh38, checksummed),
  `run.sh` (scan, evaluate, refit), `refit.py` (the quality score of
  `genotype::score`, with the AUC leaving one sample out and on a held-out
  sample) and `MANIFEST.md` (URLs and checksums of the public inputs).
  `evaluate.py` gives 95% bootstrap intervals over the samples and writes the
  true PASS calls (`--tp-out`).

### Fixed

- Container image: the `Dockerfile` copies `benches/`, which `Cargo.toml` declares (#84).
- Release workflow: the image name is lowercase (`ghcr.io/victoryt/mei-rs`), as registries require (#84).

## [0.4.0] - 2026-10-15

### Changed

- Quality score (#67): the prior and the likelihood ratios of the seven pieces
  of evidence are refit on the `PASS` calls of the seven HGSVC3 exomes (129
  true, 40 false). Leaving one sample out, AUC 0.913 to 0.928; HG002 (held
  out) 0.728 to 0.787. `--min-qual 5`: 107 to 114 true, 6 false (precision
  0.950); `--min-qual 2`: 125 true, 15 false (0.893). Evidence of the
  candidates placed by their mates does not raise the AUC and is not used.
- Assembly of candidates placed by their mates (#66): those that pass as
  insertions without junction clips (#59, #70) have their breakpoint tried as
  the left and as the right junction, and the sides whose inserted bases hold
  a poly(A/T) tail or are an element of the family are kept; an assembled
  junction replaces the estimate and drops `IMPRECISE`. Needs `--reference`.
  HGSVC3 exomes: 2 of the 13 mate-only PASS (all true) get an assembled
  junction, one at the position of the MELT call; the calls are unchanged
  (126 true, 40 false). The others have too few reads clipped by 10 bases or
  more to extend a contig.
- Mates in the poly(A) tail and reference copies (#71): a checked mate that
  is untyped and holds a run of 20 or more A or T is a mate in the tail of the
  element, with too few element bases to be typed, and is left out of the
  checked mates, unless a reference copy of the family of the candidate lies
  within 600 bp of its site (the mates of the reads in the tail of that copy
  look like it). Candidates typed by their mates or without enough junction
  clips (`no_junction`) with such a copy do not use their mate evidence
  (`mate_only` pass and typing by mates). Needs `--reference`. HGSVC3
  exomes: 121 → 126 true, 39 → 40 false (precision 0.759; `--min-qual 2`:
  0.892, `--min-qual 5`: 0.947); HG002: 16 → 17 true, 23 → 24 false;
  NA12878: 20 → 23 PASS.

## [0.3.0] - 2026-10-13

### Added

- Typing by the mates (#57): untyped candidates with a poly(A/T) tail clip
  have their mates checked too (on short reads the clips at the tail junction
  are often the tail alone), and at least 2 mates of one family type them
  (new column `typed_by`, `clips` or `mates`; VCF `TYPED_BY`). Such a
  candidate passes `untyped` with at least `--min-typed-clips` tailed clips
  and MEI mates, the mates matching the consensus with 96% identity, at
  `low` confidence. New fixture `tailonly`. HGSVC3 exomes: 99 → 111 true,
  48 → 49 false (precision 0.694, recall 0.299; `--min-qual 2`: 0.852).
- Position of the candidates placed by their mates (#58): without enough
  junction soft-clips, at the inner edge of the larger group of mate signals
  (or between the two groups when they do not overlap), with the interval
  where the insertion lies (new column `cipos`, VCF `CIPOS`). For those with
  MEI mates, the interval is read again and the boundary where at least 2
  reads are soft-clipped by 3 bases or more becomes the position (new column
  `breakpoint_clips`, VCF `BREAKPOINT_CLIPS`).
- Insertions without junction soft-clips (#59): the mates of candidates
  without junction clips (at least 3 mate signals, `signal_ratio` 0.2) are
  checked and may type them (#57); such a candidate passes `no_junction` and
  `untyped` (confidence `low`, VCF `IMPRECISE`) when typed by its mates with
  98% identity, at least 60% and 4 of its checked mates are MEI mates, and
  reads are soft-clipped by a few bases at its position (#58).

### Changed

- Mate checks (#69): a mapped mate is read from the reference over the span
  of its alignment (mate CIGAR, `MC` tag; else a read length, as before),
  and mates are typed without their runs of 10 or more A or T (the poly(A)
  tail of the element or of the reference copy the mate maps to). A mapped
  mate found in the reference within 600 bp of its candidate is not an MEI
  mate: the pair fits the reference (next to a reference copy of an element,
  the mates of the reads in that copy map with a low MAPQ to other copies).
  HGSVC3 exomes: 118 → 117 true, 49 → 39 false (precision 0.750;
  `--min-qual 2`: 0.891, `--min-qual 5`: 0.952); HG002: 16 true, 25 → 23
  false.
- Candidates typed by a single clip (#70): the mates of candidates with
  too few junction clips are checked when typed too, and they pass as
  mate-only insertions (#59) with the mean identity of their MEI mates and
  their junction clip as breakpoint evidence; one typed clip without a
  poly(A/T) tail passes `untyped` with at least 60% and 4 MEI mates matching
  the consensus with 96% identity; `low_ratio` no longer applies with 5 or
  more typed junction clips (very deep sites). HGSVC3 exomes: 117 → 121
  true, 39 false (precision 0.756; `--min-qual 2`: 0.895); within 100 bp of
  the targets 85 → 89 of 131; HG002 unchanged.

## [0.2.0] - 2026-10-11

### Added

- Element typing (#30): the junction soft-clips of each candidate, without
  the pure run of their tail, are seeded against a k-mer index (k = 11, both
  strands) of MEI consensus sequences and aligned locally (Smith-Waterman,
  affine gaps); a clip is typed when 70% of it aligns with 90% identity.
  New columns `family` (`ALU`, `LINE1`, `SVA`), `subfamily`,
  `element_strand`, `element_start`, `element_end`, `identity` and
  `typed_clips` (TSV and JSON).
- Default consensus sequences built into the binary: 42 Alu, 10 young LINE-1
  (L1HS 5' and 3' ends, L1PA2-L1PA8 3' ends, L1P1 ORF2) and 6 SVA consensus
  sequences of Dfam 4.0 (CC0), with the script that downloads them
  (`data/consensus/`); `--consensus` replaces them. (#30) The L1P1 ORF2 model
  covers the middle of LINE-1, missing from the L1HS and L1PA models: a
  5'-truncated LINE-1 whose 5' junction falls there is now typed (HGSVC3
  exomes: LINE-1 recall 6 / 32 to 14 / 32). (#56)
- `untyped` filter: fewer than `--min-typed-clips` (2) junction soft-clips
  typed to an MEI consensus. (#30)
- Fixtures: the simulated Alu is the AluY consensus of Dfam, and `positive`
  also holds a minus-strand 5'-truncated LINE-1 and an SVA insertion. (#30)
- MEI mates (#31): up to 40 mate signals per candidate (evenly sampled, for
  typed candidates with enough junction soft-clips) are checked against the MEI
  consensus sequences: unmapped mates by their read, mapped mates (low MAPQ,
  another contig, large insert) by the reference at their position, with
  `--reference`. New columns `checked_mates` and `mei_mates`, and `mei_mate`
  for each checked mate signal in JSON.
- `non_mei_mates` filter: fewer than `--min-mei-mate-fraction` (0.3) of the
  checked mates are mobile elements. (#31)
- Fixture: a chimeric junction with Alu clips and random unmapped mates in
  `positive`. (#31)
- Short soft-clips (#32): clips of `--min-short-clip-len` (15) to
  `--min-clip-len` (20) bases are kept when they have a poly(A/T) tail or are
  typed to an MEI consensus, for short reads that see a short piece of the
  junction. They count for `no_junction` and `untyped` only when the
  candidate has an MEI mate. New columns `short_junction_clips` and
  `short_typed_clips`, and `short` for each soft-clip in JSON.
- `untyped` also passes with one typed clip, a poly(A/T) tail and at least
  `--min-typed-clips` MEI mates. (#32)
- Local assembly of the junctions (#33), with `--reference`, for typed
  candidates with junction soft-clips: the reads clipped at each junction and
  the unmapped reads placed near the candidate are assembled from the
  reference flank (de Bruijn graph, k = 21; long homopolymers measured on
  the reads); the contigs give the exact junctions, the TSD and the length of
  the poly(A) tail. The assembled junctions and TSD replace the estimates
  from the clips. New columns `assembled`, `tsd` and `polya_length`, and
  `assembly` with the contigs in JSON.
- Genotype (#34): at each junction of the written candidates, the reads
  clipped there (alternative) and the reads aligned over the whole site
  (reference) give the most likely of `0/0`, `0/1` and `1/1` (binomial
  likelihoods) and its quality. New columns `genotype`, `gq`, `alt_reads`
  and `ref_reads`.
- Fixture `homozygous` (BAM and CRAM): the Alu insertion on both haplotypes.
  (#34)
- VCF output (#35): `--format vcf`, VCF 4.2 with symbolic alleles
  (`<INS:ME:ALU>`, `<INS:ME:LINE1>`, `<INS:ME:SVA>`, `<INS:ME>`), `MEINFO`,
  `SVLEN`, `TSD`, `POLYA_LEN`, confidence and support counts in `INFO`,
  `GT:GQ:AD:DP`, the quality score as `QUAL` and the filters as `FILTER`;
  bgzipped and indexed (`.tbi`) when `--output` ends with `.gz`.
- Quality score (#34): column `qual`, the phred-scaled probability that a
  candidate is not an insertion, from its evidence (naive Bayes with
  likelihood ratios measured on the HGSVC3 exomes); `low_qual` filter below
  `--min-qual` (0, off, by default).
- Validation against the HGSVC3 MEI calls (#26): `docs/validation/evaluate.py`
  (Python standard library) computes precision, recall, genotype concordance
  and precision by `--min-qual` from `mei-rs` TSV or VCF files, or the VCF of
  another caller, and lifts the truth set to GRCh37 with a UCSC chain file.
  Seven 1000 Genomes exomes: precision 0.673, recall 0.267 (0.901 and 0.221
  with `--min-qual 5`); HG002: 0.382 and 0.265; Scramble 1.0.2 on the same
  exomes: 0.687 and 0.154. Results, method and run times in
  `docs/validation/README.md`.

### Changed

- The poly(A/T) tail is no longer required to pass (`--require-tail` turns
  `no_tail` back on; `--allow-no-tail` is removed): typed clips are the
  evidence of a retrotransposon. Confidence levels: `high` (tail and TSD),
  `medium` (tail), `low` (no tail). A typed clip also keeps a candidate from
  `probe_edge`. On the NA12878 exome (`--reference`): 40 `PASS` (1 `high`, 4
  `medium`, 35 `low`) instead of 317 (2 `high`, 315 `medium`). (#30)

## [0.1.1] - 2026-10-10

Specificity on real data: `PASS` candidates need soft-clip junctions and a
poly(A/T) tail, get a confidence level, and lose the length variants of
reference poly(A) tracts; lower memory and faster CRAM reading.

### Added

- Confidence level of `PASS` candidates (`confidence` column): `high` (poly(A/T)
  tail and TSD), `medium` (tail and mate signals), `low` (soft-clips only);
  `left_junction_clips` and `right_junction_clips` columns. (#22)
- `--min-clips` and `--allow-no-tail` options. (#22)
- `--probe-edge-tolerance` option (default 2 bp). (#29)
- `tests/data/targets_tiled.bed`: the fixture insertion with probe edges on
  both junctions. (#29)
- The log counts, among the soft-clips set aside for their base quality,
  those without any base at the threshold (masked read ends): 2,682,180 of
  2,870,963 on the NA12878 exome. The quality filter is unchanged: trimming
  the low-quality read end would rescue 0.05% of them. (#25)
- `tandem_repeat` filter, applied when `--reference` is given: the soft-clips
  of the candidate are a length variant of a reference homopolymer or short
  tandem repeat (typically a longer poly(A) tract at the 3' end of an Alu of
  the reference), not an insertion. Such clips have a poly(A/T) tail and the
  two junctions overlap like a TSD, so they made most `high` candidates.
  On the NA12878 exome: 317 `PASS` (2 `high`, 315 `medium`) instead of 358
  (6 `high`, 350 `medium`, 2 `low`); 4 of the 6 `high` were such variants.
  Without `--reference`, a warning says the filter is off. (#42)
- `tests/data/targets_tract.bed`: an expanded reference poly(A) tract in the
  `positive` fixture. (#42)

### Changed

- `PASS` now requires soft-clips at the junctions (new `no_junction` filter)
  and a soft-clip with a poly(A/T) tail (new `no_tail` filter); a TSD is only
  reported when each junction has `--min-clips` soft-clips. On the NA12878
  exome: 358 `PASS` candidates instead of 56,956. (#22)
- `probe_edge` now applies within `--probe-edge-tolerance` of a probe edge,
  and only to candidates without any poly(A/T) tail or TSD: a real insertion
  next to a probe edge stays `PASS`. (#29)
- Faster CRAM reading: neighbouring regions (and candidate breakpoints, for
  the local depth) are read with one query, so that each CRAM container is
  decoded once instead of once per region; regions are processed contig by
  contig, with the reference sequence loaded once per contig and released
  after it; the index is read once instead of once per reader; fewer readers
  are opened. Same output. On the NA12878 exome (12 cores): CRAM 507 s and
  6.5 GB -> 124 s (memory 1.4 GB with the macOS large-allocation cache
  disabled, 3.8 GB with it); BAM 72 s and 1.65 GB -> 66 s and 1.2 GB. (#24)
- Template lengths of proper pairs are counted in a bounded histogram instead
  of being kept in memory: same insert size estimate and output, peak memory
  on the NA12878 exome down from 3.2 GB to 1.7 GB. (#23)

## [0.1.0] - 2026-10-10

First release: discovery of candidate mobile element insertion sites in
targeted capture and exome data.

### Added

- `mei-rs scan`, with global `--threads`, `-v` and `-q` options. (#1, #2)
- Inputs: indexed BAM or CRAM (with `--reference`), capture BED (plain or
  gzipped) padded with `--padding` and merged, `--region`. Inconsistent
  inputs fail early with a fix: contig naming (`chr1` vs `1`), missing index,
  reference contigs or lengths, empty region. (#5, #6, #7)
- Read filtering: secondary, supplementary, QC-fail and duplicate reads
  (`--include-duplicates`) and reads below `--min-mapq` do not anchor
  signals. (#8)
- Signal extraction, in parallel and deterministic: soft-clips with their
  poly(A/T) tail (`--min-clip-len`, `--min-clip-quality`, `--polya-min-len`,
  `--polya-min-frac`), unmapped mates, low-MAPQ mates, mates on another
  contig, and large inserts beyond a threshold estimated from the sample.
  (#9, #10, #11)
- Clustering into candidates (`--window`, `--min-support`), with the
  breakpoint, the junctions of both flanks, the TSD length and the strand.
  (#12)
- Capture filters: `probe_edge` (soft-clips stopping at probe edges) and
  `low_ratio` (`--min-signal-ratio`, against the local depth);
  `--keep-filtered`. (#13)
- Outputs: TSV, BED and JSON (`--format`, `--output`), specified in
  `docs/output.md`. (#14)
- Synthetic fixtures (a heterozygous Alu insertion with a TSD and a poly(A)
  tail, a negative sample, a probe-edge artefact) and end-to-end tests. (#4,
  #15)
- GitHub Actions CI: lint, tests on Linux and macOS, MSRV 1.91, audit. (#3)

### Known limitations

- Not yet specific on real data: on the NA12878 exome (1000 Genomes,
  GRCh38), 56,956 `PASS` candidates, mostly mate-only clusters without a
  soft-clip junction; to be fixed in v0.1.1 (#22). See the README.
- CRAM is about 7 times slower than BAM (#24); the insert size sample is
  kept in memory (#23); most soft-clips of older data fall below the base
  quality threshold (#25).

[Unreleased]: https://github.com/VictorYt/Mei-rs/compare/v0.4.0...HEAD
[0.4.0]: https://github.com/VictorYt/Mei-rs/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/VictorYt/Mei-rs/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/VictorYt/Mei-rs/compare/v0.1.1...v0.2.0
[0.1.1]: https://github.com/VictorYt/Mei-rs/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/VictorYt/Mei-rs/releases/tag/v0.1.0
