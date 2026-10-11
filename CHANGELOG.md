# Changelog

All notable changes to Mei-rs are documented in this file. The format is
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
project follows [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- Element typing (#30): the junction soft-clips of each candidate, without
  the pure run of their tail, are seeded against a k-mer index (k = 11, both
  strands) of MEI consensus sequences and aligned locally (Smith-Waterman,
  affine gaps); a clip is typed when 70% of it aligns with 90% identity.
  New columns `family` (`ALU`, `LINE1`, `SVA`), `subfamily`,
  `element_strand`, `element_start`, `element_end`, `identity` and
  `typed_clips` (TSV and JSON).
- Default consensus sequences built into the binary: 42 Alu, 9 young LINE-1
  (L1HS 5' and 3' ends, L1PA2-L1PA8 3' ends) and 6 SVA consensus sequences of
  Dfam 4.0 (CC0), with the script that downloads them
  (`data/consensus/`); `--consensus` replaces them. (#30)
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

[Unreleased]: https://github.com/VictorYt/Mei-rs/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/VictorYt/Mei-rs/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/VictorYt/Mei-rs/releases/tag/v0.1.0
