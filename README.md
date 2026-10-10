# Mei-rs

An ultra-fast, memory-efficient Rust tool for detecting mobile element
insertions (Alu, LINE-1, SVA) in targeted capture and exome sequencing data.

Mei-rs reads only the capture targets of an indexed BAM or CRAM file, collects
the reads that betray an insertion (soft-clipped reads, discordant pairs,
unmapped mates), groups them into candidate insertion sites, types the
inserted element against Alu, LINE-1 and SVA consensus sequences, and removes
the artefacts specific to capture data, such as soft-clips that stop at the
probe edges.

> **Status: v0.2 in progress.** `mei-rs scan` finds candidate insertion
> sites, their junctions, target site duplication (TSD), strand and element
> (family, subfamily). Local assembly, genotyping and VCF output come in the
> next releases (see [Limitations](#limitations)).

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
  --reference genome.fa \
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
| `--reference` | reference genome indexed with `samtools faidx` (`.fai`); required for CRAM, and recommended for BAM: without it, the `tandem_repeat` filter is off. When given, its contigs must match the lengths of the alignment header. |
| `--consensus` | MEI consensus sequences (FASTA, headers `>name#SINE/Alu`, `>name#LINE/L1`, `>name#Retroposon/SVA`). Default: the Alu, young LINE-1 and SVA consensus sequences of Dfam (CC0) built into mei-rs, see [data/consensus](data/consensus/README.md). |

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
| `--min-clips` | 2 | minimum number of soft-clips at the junctions (`no_junction` below); each junction needs as many to give a TSD |
| `--min-typed-clips` | 2 | minimum number of junction soft-clips typed to an MEI consensus (`untyped` below; 0 turns the filter off) |
| `--require-tail` | off | also require a soft-clip with a poly(A/T) tail (`no_tail`) |
| `--min-signal-ratio` | 0.05 | minimum number of supporting reads per read of local depth |
| `--probe-edge-tolerance` | 2 | distance in bp to a probe start or end within which a soft-clip stops at the probe edge |
| `--keep-filtered` | off | also write the candidates flagged `no_junction`, `no_tail`, `untyped`, `tandem_repeat`, `probe_edge` or `low_ratio` |
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
3. **Typing**: each junction soft-clip, without the pure run of its tail, is
   seeded against a k-mer index of the consensus sequences (both strands) and
   aligned locally (Smith-Waterman) to the best ones; a clip is typed when at
   least 70% of it aligns with 90% identity or more. The candidate takes the
   family and subfamily with the highest summed score, the strand of the
   element and the part of the consensus covered.
4. **Filters**: `no_junction` when fewer than `--min-clips` soft-clips mark
   the junctions (mate signals alone); `untyped` when fewer than
   `--min-typed-clips` of them are typed; `no_tail` (with `--require-tail`)
   without any poly(A/T) tail;
   `tandem_repeat` (with `--reference`) when the soft-clips are the reference
   with a longer homopolymer or tandem repeat, such as the poly(A) tract of
   an Alu of the reference, rather than an insertion;
   `probe_edge` when every soft-clip stops within `--probe-edge-tolerance`
   of a probe start or end, without any poly(A/T) tail, typed clip or TSD
   (digestion or ligation artefact of the capture); `low_ratio` when the
   supporting reads are too few for the local depth.
5. **Confidence** of the `PASS` candidates: `high` with a poly(A/T) tail and a
   TSD (complete evidence of a retrotransposition), `medium` with a tail,
   `low` without a tail (typed clips only).

### Output

One line per candidate, sorted by position. Positions are 1-based: a position
`p` means the insertion is right after base `p`. Example:

```
#id    contig  position  left_junction  right_junction  tsd_length  strand  filter  confidence  support  ...
MEI_1  chr1    10012     10012          10000           12          +       PASS    high        160      ...
```

Every column of the TSV, BED and JSON formats is described in
[docs/output.md](docs/output.md).

## Limitations

v0.1.0 found candidate sites but was not specific on real data. On the
NA12878 exome of the 1000 Genomes Project (GRCh38, 172 million reads,
1000 Genomes exome targets), it wrote 56,956 `PASS` candidates, mostly
clusters of unmapped mates or of mates on another contig without any soft-clip
junction, while a few dozen insertions are expected in the targets.

Since #22, a `PASS` candidate needs soft-clips at its junctions
(`no_junction`), and since #30, at least two of them typed to an Alu, LINE-1
or SVA consensus (`untyped`): the poly(A/T) tail, often missed by short reads,
only sets the confidence level. With `--reference`, the `tandem_repeat` filter
(#42) rejects the longer poly(A) tracts of reference Alus, whose clips have a
tail and junctions that overlap like a TSD. On the same exome: v0.1.1 wrote
317 `PASS` (2 `high`, 315 `medium`), this version **40** (1 `high`, 4
`medium`, 35 `low`). `high` and `medium` candidates (typed, with a tail) are
the most reliable; `low` ones (typed clips without a tail) still include
artefacts until the mates are checked against MEI sequences (#31).

Other limitations:

- No local assembly, genotype or VCF yet.
- Typing uses the junction soft-clips only: a candidate whose clips are
  shorter than 15 bases (without their tail) stays untyped. Close subfamilies
  (e.g. AluY and AluYa5) are often not told apart, and the middle of LINE-1
  (most of ORF2) is not in the default Dfam set.
- Candidates closer than `--window` are merged into one.
- The TSD length comes from the soft-clip junctions only, and the length of a
  poly(A/T) tail is approximate (it can extend over a few bases after it).
- On older data, most soft-clips fall below `--min-clip-quality`: 2.9 million
  set aside against 105,000 kept on NA12878. They are not lost junctions:
  93% have no base at the threshold (read ends masked by the base caller,
  soft-clipped because they mismatch the reference), and trimming the
  low-quality read end, as BWA `-q` does, would rescue 0.05% of them (chr1).
  The log reports both counts (#25).
- CRAM is slower than BAM: decoding its blocks costs more than reading BAM
  (NA12878 on 12 cores, same output: BAM 66 s and 1.2 GB, CRAM 124 s; CRAM
  was 507 s and 6.5 GB before #24). On macOS, the allocator keeps the freed
  reference sequences in memory (3.8 GB, 1.4 GB with `MallocLargeCache=0`).

### Next steps (v0.2.0)

- Mates checked against the MEI consensus sequences (#31), so that the `low`
  level gets specific.
- Short soft-clips of 2x75 bp reads (#32), local assembly of the junctions
  and TSD (#33), scoring and genotype (#34), VCF output (#35).
- Precision and recall against truth sets (HGSVC, GIAB HG002) and comparison
  with Scramble, MELT and xTea (#26).

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
