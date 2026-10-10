# Changelog

All notable changes to Mei-rs are documented in this file. The format is
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
project follows [Semantic Versioning](https://semver.org/).

## [Unreleased]

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

### Changed

- `PASS` now requires soft-clips at the junctions (new `no_junction` filter)
  and a soft-clip with a poly(A/T) tail (new `no_tail` filter); a TSD is only
  reported when each junction has `--min-clips` soft-clips. On the NA12878
  exome: 358 `PASS` candidates instead of 56,956. (#22)
- `probe_edge` now applies within `--probe-edge-tolerance` of a probe edge,
  and only to candidates without any poly(A/T) tail or TSD: a real insertion
  next to a probe edge stays `PASS`. (#29)
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

[Unreleased]: https://github.com/VictorYt/Mei-rs/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/VictorYt/Mei-rs/releases/tag/v0.1.0
