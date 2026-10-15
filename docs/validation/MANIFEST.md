# Data manifest

Every public input of the evaluation, with the checksum of the copy used for the
numbers of the README (v0.4.0). `fetch-data.sh DIR` downloads and verifies them
(about 70 GB, most of it the reference genomes and the HG002 alignment);
`run.sh DATA OUT` then reruns the evaluation and the refit.

| Input | Source | md5 |
| :--- | :--- | :--- |
| HGSVC3 MEI calls, GRCh38 | `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Mobile_Elements/1.0/MEI_Callset_GRCh38.ALL.20241211.csv.gz` | `0131ecb529c9225dab677aced8d72b0a` |
| UCSC chain `hg38ToHg19` (HG002 truth lift) | `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz` | `ff3031d93792f4cbb86af44055efd903` |
| UCSC chain `hg19ToHg38` (targets lift) | `https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz` | `35887f73fe5e2231656504d1f6430900` |
| 1000 Genomes exome pull-down targets, GRCh37 | `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/exome_pull_down_targets/20130108.exome.targets.bed` | `3285b0ba60e45550297d5e072d105483` |
| GRCh38 reference (full analysis set, decoys, HLA) | `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa` | `64b32de2fc934679c16e83a2bc072064` |
| GRCh37 reference hs37d5 | `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz` | `d10eebe06c0dbbcb04253e3294d63efc` |
| HG002 exome (Oslo University Hospital, GRCh37) | `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam` (+ `.bai`) | `c80f0cab24bfaa504393457b8f7191fa` |

Derived targets (`build-targets.sh`, CrossMap 0.7.3):

| File | Lines | md5 |
| :--- | ---: | :--- |
| `exome_targets.GRCh37.bed` (pull-down targets, contigs `1`..`22`) | 193,637 | `14805f7965dc8b44772b410d1141393c` |
| `exome_targets.hg38.bed` (lifted, mapped in one piece, sorted) | 194,644 | `f4f589406344a74b624e1af7d8ff8e33` |

The UCSC chain files are served as they are today: if UCSC replaces one, the
script stops on the checksum instead of changing the targets silently.

## 1000 Genomes exome alignments (GRCh38)

HG00512, HG00513, HG00731, HG00733, NA19238, NA19239 and NA19240: the
`*.exome.cram` and `.crai` listed, with their md5, in
`https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/1000genomes.exome.GRCh38DH.alignment.index`
(`fetch-data.sh` reads the checksums from this index). The numbers of the
README were measured on the BAM conversions of these CRAMs (same records); the
calls of CRAM and BAM are identical.

## Not public

NA12878 (timing and memory check of each release) is a 1000 Genomes exome that
is not genotyped in HGSVC3: it is run, never scored.

## Order of the steps

```bash
docs/validation/fetch-data.sh data/
cargo build --release
docs/validation/run.sh data/ out/
```

`run.sh` writes in `out/`: the TSV of every sample, `evaluation.txt` and
`evaluation_hg002.txt` (`evaluate.py`, with the 95% bootstrap intervals over the
samples), the true and false PASS calls (`tp.tsv`, `fp.tsv`), the missed
insertions (`fn.tsv`) and `refit.txt` (`refit.py`: the values of
`src/genotype/score.rs`, the AUC leaving one sample out and on HG002).
