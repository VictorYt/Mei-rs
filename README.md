# Mei-rs

An ultra-fast, memory-efficient Rust tool for detecting mobile element
insertions (Alu, LINE-1, SVA) in targeted capture and exome sequencing data.

Mei-rs reads only the capture targets of an indexed BAM or CRAM file, collects
the reads that betray an insertion (soft-clipped reads, discordant pairs,
unmapped mates), groups them into candidate insertion sites, and removes the
artefacts specific to capture data, such as soft-clips that stop at the probe
edges.

> **Status: v0.1, candidate discovery.** `mei-rs scan` finds candidate
> insertion sites, their junctions, target site duplication (TSD) and strand.
> Element typing (Alu / LINE-1 / SVA), local assembly, genotyping and VCF
> output come in the next releases (see [Limitations](#limitations)).

## Installation

From source, with Rust 1.91 or later:

```bash
cargo install --git https://github.com/VictorYt/Mei-rs --locked
# or, from a clone:
cargo build --release   # binary in target/release/mei-rs
```

The binary has no runtime dependency: BAM, CRAM, BED and FASTA files are read
in pure Rust.

## Usage

```bash
mei-rs scan \
  --bam sample.bam \
  --bed targets.bed \
  --output sample.candidates.tsv \
  --threads 8

# CRAM: the reference is required
mei-rs scan --bam sample.cram --reference genome.fa --bed targets.bed
```

The candidates are written to standard output unless `--output` is given; the
logs (inputs, read counts, insert size, signals, candidates) go to standard
error.

### Inputs

| Input | Requirements |
|---|---|
| `--bam` | BAM or CRAM, coordinate-sorted and indexed (`.bam.bai`, `.bam.csi`, `.bai`, `.cram.crai` or `.crai`). The sample name is the `SM` of the read groups, or else the file name. |
| `--bed` | capture targets: at least 3 columns (contig, 0-based start, end), tabs or spaces, plain or gzipped. `#`, `track` and `browser` lines are skipped. The contig names must match the alignment (`chr1` vs `1` is reported). |
| `--reference` | reference genome indexed with `samtools faidx` (`.fai`); required for CRAM. When given, its contigs must match the lengths of the alignment header. |

Inconsistent inputs fail early with a message that says how to fix them.
Targets on contigs absent from the alignment, or beyond their contig, are
skipped with a warning.

### Options

| Option | Default | Meaning |
|---|---|---|
| **Input** | | |
| `--padding` | 300 | bases added on each side of every target; padded targets that overlap are merged |
| `--region` | | only scan this region (`chr`, or `chr:start-end`, 1-based inclusive) |
| **Read filtering** | | |
| `--min-mapq` | 20 | minimum mapping quality of an anchored read; mates below it are low-MAPQ mates |
| `--include-duplicates` | off | keep reads flagged as duplicates (files without duplicate marking) |
| **Signal extraction** | | |
| `--min-clip-len` | 20 | minimum length of a soft-clip |
| `--min-clip-quality` | 20 | minimum median base quality of a soft-clip |
| `--polya-min-len` | 10 | minimum length of a poly(A/T) tail at the junction end of a soft-clip |
| `--polya-min-frac` | 0.8 | minimum fraction of A (or T) in a tail, to tolerate sequencing errors |
| **Clustering and filters** | | |
| `--window` | 100 | maximum distance between two consecutive signals of a candidate |
| `--min-support` | 3 | minimum number of distinct supporting reads of a candidate |
| `--min-signal-ratio` | 0.05 | minimum number of supporting reads per read of local depth |
| `--keep-filtered` | off | also write the candidates flagged `probe_edge` or `low_ratio` |
| **Output** | | |
| `-o`, `--output` | stdout | output file |
| `--format` | `tsv` | `tsv`, `bed` or `json` (with every signal of every candidate) |
| **Global** | | |
| `-t`, `--threads` | all cores | worker threads |
| `-v`, `-q` | | more (`-v`, `-vv`) or fewer (`-q`) logs; `RUST_LOG` takes precedence |

### How it works

1. **Signal extraction**, in parallel over the padded targets: reads that are
   primary, mapped, not duplicates and above `--min-mapq` can anchor a signal.
   - **Soft-clips** of at least `--min-clip-len` bases: the junction gives the
     breakpoint to the base, and a poly(A) or poly(T) tail at the junction
     gives the strand of the element.
   - **Mate signals**: the mate is unmapped, maps with a low MAPQ, maps to
     another contig, or maps unusually far away on the same contig (beyond the
     median insert size + 6 × 1.4826 × MAD of the sample).
2. **Clustering**: signals at most `--window` bases apart form a candidate.
   Its position is the most supported junction; the junctions of the two
   flanks give the TSD.
3. **Capture filters**: `probe_edge` when every soft-clip stops exactly at a
   probe start or end (digestion or ligation artefact of the capture);
   `low_ratio` when the supporting reads are too few for the local depth.

### Output

One line per candidate, sorted by position. Positions are 1-based: a position
`p` means the insertion is right after base `p`. Example:

```
#id    contig  position  left_junction  right_junction  tsd_length  strand  filter  support  ...
MEI_1  chr1    10012     10012          10000           12          +       PASS    160      ...
```

Every column of the TSV, BED and JSON formats is described in
[docs/output.md](docs/output.md).

## Limitations

- No element typing (Alu / LINE-1 / SVA), local assembly, genotype or VCF yet.
- Candidates closer than `--window` are merged into one.
- The TSD length comes from the soft-clip junctions only, and the length of a
  poly(A/T) tail is approximate (it can extend over a few bases after it).
- Validated on synthetic data so far (see [tests/data](tests/data/README.md)).

Calls are research results: confirm them with an orthogonal method before any
clinical use.

## Development

```bash
cargo build --release
cargo test                  # unit and end-to-end tests on synthetic fixtures
tests/data/generate.sh      # regenerate the fixtures (python3, samtools, perl)
UPDATE_EXPECTED=1 cargo test --test scan   # update the expected outputs
```

See [CHANGELOG.md](CHANGELOG.md) for the release history.

## License

[GNU Affero General Public License v3.0](LICENSE) or later.
