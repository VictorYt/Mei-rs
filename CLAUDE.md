# CLAUDE.md - Mei-rs (Mobile Element Insertion detection in Rust)

## Executive Overview
**Mei-rs** is a high-performance detector, written in Rust, of mobile element insertions (MEI: **Alu**, **LINE-1**, **SVA**) in targeted sequencing data (**capture panels / exomes**).
It uses abnormal read signals (soft-clips, discordant pairs, unmapped mates) restricted to the capture regions, types the inserted element, resolves the insertion site to the nucleotide (TSD, poly(A) tail) and writes a **standard VCF**.

Clinical goal: accuracy comparable to **xTea / Scramble / MELT** on capture data, where those tools fit poorly (uneven coverage, probe-edge effects), with a single fast binary.

---

## 1. Analysis pipeline

```
BAM/CRAM + BED (capture) + MEI consensus (FASTA) + reference genome (FASTA)
                                │
                                ▼
      Phase 1: Targeted signal extraction (soft-clips / discordant)    → reader
                                │
                                ▼
      Phase 2: Spatial clustering & capture artefact filtering         → cluster
                                │
                                ▼
      Phase 3: Element typing (Alu / LINE-1 / SVA) & orientation       → typer
                                │
                                ▼
      Phase 4: Local assembly & TSD / poly(A) characterisation         → assembly
                                │
                                ▼
      Phase 5: Genotyping, scoring, filtering & VCF export             → genotype / writer
```

### Target command (UX)
A single calling command; debug outputs are optional.

```bash
mei-rs call \
  --bam sample.bam \               # or .cram (+ .bai / .crai)
  --bed targets.bed \
  --padding 300 \
  --consensus mei_consensus.fa \   # AluY, AluSx, L1Hs, SVA-E/F...
  --reference genome.fa \          # required for CRAM and TSD detection
  --output sample.mei.vcf.gz \
  --debug-out sample.mei.json \    # optional: JSON or BED of the candidates
  --threads 16
```

The v0.1.0 milestone ships the `mei-rs scan` subcommand first (phases 1 and 2, candidates as TSV/BED/JSON); `call` (typing, assembly, VCF) follows in v0.2.

### Inputs / outputs
* **Inputs:** `ALIGNMENT.bam/.cram` + `.bai/.crai` index; `TARGETS.bed` (capture regions, ~300 bp padding); `MEI_CONSENSUS.fasta` (Alu, L1Hs, SVA-E/F consensus sequences); reference genome FASTA (+ `.fai`).
* **Outputs:** `OUTPUT.vcf.gz` (annotated VCF 4.2+); optional `OUTPUT.json` or `OUTPUT.bed` for debugging.

---

## 2. Module specification

### Phase 1 — `reader`: selective ingestion and signal extraction
Goal: discard ~95% of irrelevant reads and keep only structural anomalies.
* **BED-driven processing:** BED loaded into an interval tree (`coitrees` or `bio::data_structures::interval_tree`); indexed BAM queries restricted to the targets + padding (no intron scanning).
* **Soft-clipped reads (split reads):** `S` operation in the CIGAR (e.g. `35S115M`). Inspect the clipped part: poly(A)/poly(T) tail or repetitive sequence?
* **Discordant pairs:** read 1 anchored in the target, read 2 unmapped, mapped with near-zero MAPQ, or mapped elsewhere in the genome.
* **Unmapped reads with an anchored mate:** keep unmapped reads whose mate is anchored in a target.

### Phase 2 — `cluster`: spatial clustering and noise reduction
* **Aggregation:** group split reads and discordant pairs by proximity (50–100 bp sliding window or DBSCAN); the candidate breakpoint is the densest soft-clip junction.
* **Probe-edge effect:** if 100% of the soft-clips in a window stop exactly at the start or end of a BED probe, flag the region as a capture artefact (digestion/ligation bias).
* **Dynamic coverage threshold:** require a minimum ratio of MEI signals to local depth.

### Phase 3 — `typer`: retrotransposon typing
* **Classification:** align the clipped sequences (or the assembled contig) to the consensus sequences (AluY, AluSx, L1Hs, SVA-E...) with a k-mer pre-filter then vectorised local alignment (Smith-Waterman / Myers bit-parallel) → family and subfamily.
* **Strand / orientation:** poly(A) vs poly(T) tail → insertion strand (`+` or `-`).

### Phase 4 — `assembly`: fine resolution of the insertion site
* **Local de novo assembly:** gather anchored and clipped reads around the candidate site; small de Bruijn graph (k = 21 or 31, `petgraph` if useful) to assemble the junction contig.
* **TSD (Target Site Duplication):** compare the contig with the reference genome to find the target site duplication (usually 7–20 bp) on both sides of the insertion. A valid TSD confirms the retrotransposition mechanism (TPRT) and rules out complex rearrangements.

### Phase 5 — `genotype` / `writer`: genotyping, scoring and VCF
* **Genotype:** reads supporting the insertion (alt) vs reads spanning the reference site without the insertion (ref) → `0/1`, `1/1` or `0/0` (noise) with a simple Bayesian model or a likelihood test.
* **QUAL / FILTER**, score based on:
  * supporting reads on both sides (5' and 3' junctions);
  * poly(A/T) tail detection;
  * a clear TSD;
  * diversity of fragment start positions (PCR duplicate removal).
* **VCF:** spec-compliant (symbolic ALTs `<INS:ME:ALU>`, `<INS:ME:LINE1>`, `<INS:ME:SVA>`; INFO `SVTYPE=INS`, `MEINFO=AluYa5,1,280,+`, `SVLEN`, `TSD=AAGGACT...`); declare every INFO/FORMAT/FILTER field in the header.

---

## 3. Project structure & tech stack

```
mei-rs/
├── Cargo.toml
├── src/
│   ├── main.rs              # CLI entry point (clap)
│   ├── cli.rs               # Subcommands and arguments
│   ├── inputs.rs            # Loading and consistency checks of the inputs
│   ├── scan.rs              # `scan` subcommand (phases 1-2)
│   ├── regions/             # --region, BED loading, padding, merging, interval tree
│   ├── reader/              # Phase 1: indexed BAM/CRAM reading, signal extraction
│   ├── cluster/             # Phase 2: windowing, breakpoints, capture filters
│   ├── typer/               # Phase 3: k-mers + alignment to consensus, strand
│   ├── assembly/            # Phase 4: local de Bruijn, TSD, poly(A)
│   ├── genotype/            # Phase 5: genotype, QUAL/FILTER
│   └── writer/              # VCF 4.2+, debug JSON/BED
└── tests/
```

### Key dependencies (Cargo.toml)
* **Genomic I/O:** `noodles` + `noodles-util` (pure Rust, static binary; one indexed reader for BAM and CRAM), `flate2` (gzipped BED).
* **Parallelism:** `rayon` (parallel iterators over chromosomes / BED regions).
* **Algorithms:** `coitrees` (probe interval trees), `bio` (Rust-Bio: Smith-Waterman, string algorithms), `petgraph` (light de Bruijn graphs if needed).
* **Concurrency / cache:** `dashmap` or `lru` for shared in-memory structures.
* **CLI, errors & logging:** `clap`, `thiserror` (module errors), `anyhow` (`main`), `tracing`.

---

## 4. Roadmap (prototyping)

| Sprint | Duration | Scope |
| :--- | :--- | :--- |
| **1 — PoC** | 2 weeks | BAM ingestion + extraction of soft-clipped reads with poly(A/T) on a targeted BED region. |
| **2 — Clustering & alignment** | 2 weeks | Spatial windowing + alignment of clipped sequences to the Alu consensus (`bio::alignment`). |
| **3 — Capture bias & assembly** | 3 weeks | Probe-edge filtering + small assembly for TSD detection. |
| **4 — VCF & benchmark** | 2 weeks | VCF generation, run-time benchmarks, validation against **Scramble** and **MELT** on control data (HG002 / GIAB). |

---

## 5. Development Guidelines & Commands

### Build & tests
```bash
# Optimised release build (LTO enabled)
cargo build --release

# Unit and integration tests
cargo test

# Benchmarks (BAM reading, clustering, alignment)
cargo bench
```

### Style & linting
* Everything tracked in the repository (code, comments, docs, `CLAUDE.md`) and everything on GitHub (issues, pull requests, milestones, labels, commit messages) is written in English.
* Follow standard Rust naming conventions (strict `clippy`).
* Format the code with `cargo fmt` before every commit.
* Document every public function (`///`).
* Coordinates: always state the system (BED 0-based half-open, VCF 1-based) at module boundaries.
