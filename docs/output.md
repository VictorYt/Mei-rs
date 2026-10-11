# Output of `mei-rs scan`

`mei-rs scan` writes one line (or JSON object) per candidate insertion site,
sorted by contig (alignment header order) and position. By default only the
candidates that pass every filter are written; `--keep-filtered` also writes
those flagged `no_junction`, `no_tail`, `untyped`, `non_mei_mates`,
`tandem_repeat`, `probe_edge` or `low_ratio`. Clusters with fewer than
`--min-support` distinct reads are never written.

## Coordinates

TSV and JSON positions are **1-based**: a position `p` means the insertion is
right after base `p`. For an insertion with a target site duplication (TSD),
the duplicated bases are `right_junction + 1` to `left_junction`.

The BED output is 0-based half-open, like every BED file.

## TSV (`--format tsv`, default)

Header lines start with `##` (version, command line, sample, coordinate
convention), then one `#` line names the columns:

| Column | Description |
|---|---|
| `id` | Candidate identifier (`MEI_1`, `MEI_2`, ... in output order). |
| `contig` | Contig. |
| `position` | Estimated insertion point: the most supported soft-clip junction, or, without soft-clips, the midpoint between the reads of the two flanks. |
| `left_junction` | End of the upstream flank: the assembled junction (`assembled` is `both` or `left`), or else the most supported junction of the reads clipped on their right end (`.` if none). |
| `right_junction` | Start of the downstream flank: the assembled junction (`assembled` is `both` or `right`), or else the most supported junction of the reads clipped on their left end (`.` if none). |
| `tsd_length` | Length of the target site duplication: with both junctions assembled, the length of `tsd` (`.` when the assembled flanks do not overlap by 1 to 50 bp); otherwise estimated from the clips, `left_junction - right_junction` when each junction has at least `--min-clips` soft-clips and the flanks overlap by 1 to 50 bp (`.` otherwise). |
| `strand` | Strand of the inserted element: `+` when poly(A) tails outnumber poly(T) tails in the soft-clips, `-` for the reverse, `.` on ties. |
| `filter` | `PASS`, or the failed filters separated by `;`: `no_junction` (fewer than `--min-clips` soft-clips at the junctions, short ones counting only when the candidate has an MEI mate: mate signals alone do not resolve the insertion to the base), `no_tail` (no soft-clip with a poly(A/T) tail; only with `--require-tail`), `untyped` (fewer than `--min-typed-clips` junction soft-clips typed to an MEI consensus, short ones counting only when the candidate has an MEI mate, and not one typed clip with a poly(A/T) tail and at least as many MEI mates: nothing shows that the inserted sequence is a mobile element), `non_mei_mates` (mates were checked against the consensus sequences and fewer than `--min-mei-mate-fraction` of them are mobile elements: the mate signals are mostly noise), `tandem_repeat` (with `--reference` only: the soft-clips are a length variant of a reference homopolymer or tandem repeat, such as a longer poly(A) tract of an Alu of the reference: at least 2 junction soft-clips, and half of the informative ones, end in the reference where a longer tract would put them, or the tailed soft-clips hold only their tail and the reference has a homopolymer of at least 10 bp of the tail base at their junction), `probe_edge` (every soft-clip stops within `--probe-edge-tolerance` of a probe start or end, none has a poly(A/T) tail, none is typed and there is no TSD: capture artefact), `low_ratio` (`signal_ratio` below `--min-signal-ratio`), `low_qual` (`qual` below `--min-qual`; off by default). |
| `confidence` | For `PASS` candidates (typed, unless `--min-typed-clips 0`): `high` (a poly(A/T) tail and a TSD: complete evidence of a retrotransposition), `medium` (a poly(A/T) tail, without a TSD), `low` (no poly(A/T) tail: often the tailed junction is not covered by clips long enough); `.` for filtered candidates. |
| `support` | Distinct supporting reads (by name). |
| `left_clips` | Reads soft-clipped on their left end (the insertion lies before them). |
| `right_clips` | Reads soft-clipped on their right end (the insertion lies after them). |
| `left_junction_clips` | Soft-clips within 2 bp of `left_junction`. |
| `right_junction_clips` | Soft-clips within 2 bp of `right_junction`. |
| `short_junction_clips` | Of the soft-clips of both junctions, short ones (`--min-short-clip-len` to `--min-clip-len` - 1 bp, kept with a poly(A/T) tail or MEI sequence); they count for `no_junction` and `untyped` only when `mei_mates` is at least 1. |
| `unmapped_mates` | Reads whose mate is unmapped. |
| `discordant_mates` | Reads whose mate maps with a low MAPQ, to another contig, or far away on the same contig. |
| `polya_clips` | Soft-clips with a poly(A) tail at the junction. |
| `polyt_clips` | Soft-clips with a poly(T) tail at the junction. |
| `both_sides` | `yes` if signals come from both flanks of the insertion. |
| `distinct_starts` | Distinct fragment start positions (strand-aware) of the supporting reads; far below `support`, it points to PCR duplicates. |
| `depth` | Reads (after filtering) covering one of the two bases around `position`. |
| `signal_ratio` | `support / depth` (3 decimals). Mate signals come from reads that do not cover `position`, so it can exceed 1. |
| `start`, `end` | Positions of the first and last signal of the cluster. |
| `family` | Family of the inserted element, typed from the junction soft-clips: `ALU`, `LINE1` or `SVA` (`.` if no clip is typed). |
| `subfamily` | Best-scoring consensus of that family, e.g. `AluYa5`, `L1HS_3end`, `SVA_F`. Close subfamilies differ by a few bases: short clips often cannot tell them apart. |
| `element_strand` | Strand of the element on the reference: `+` when the clips match the consensus, `-` when they match its reverse complement. It should agree with `strand` (from the tails) when both are known. |
| `element_start`, `element_end` | Part of the `subfamily` consensus covered by its clips, 1-based inclusive (e.g. `1` for a clip of the 5' end of the element). |
| `identity` | Mean identity of the clips aligned to the `subfamily` consensus (3 decimals). |
| `typed_clips` | Junction soft-clips typed to `family`. |
| `short_typed_clips` | Of these, short soft-clips. |
| `checked_mates` | Mate signals whose mate sequence was checked against the consensus sequences: unmapped mates (their read), and, with `--reference`, mapped mates (the reference at the mate position). At most 40 per candidate, evenly sampled; only for typed candidates with enough junction soft-clips (not `no_junction`), the others cannot pass. |
| `mei_mates` | Of these, mates typed to `family` (to any family if the candidate is untyped). |
| `assembled` | Junctions assembled (with `--reference`, for typed candidates with junction soft-clips): `both`, `left`, `right` or `none`; `.` when not assembled. |
| `tsd` | With both junctions assembled, the TSD: the reference bases between the junctions (`.` otherwise, or when the flanks do not overlap). When the inserted bases next to a junction match the reference beyond it (microhomology), they are counted in the TSD, as aligners and PALMER do. |
| `polya_length` | Length of the poly(A) (or poly(T)) tail at an assembled junction: the run of A (T) from the junction, when the contig holds its end (`.` otherwise). |
| `qual` | Quality score: phred-scaled probability that the candidate is not a mobile element insertion, from its evidence (naive Bayes: a poly(A/T) tail, mostly MEI mates (60%), typed clips of 98% identity, distinct fragment starts (70% of the supporting reads), `signal_ratio` of 0.2, signals on both sides, a TSD), with likelihood ratios measured on the `PASS` candidates of 7 HGSVC3 exomes. It ranks the candidates better than it estimates the probability. |
| `genotype` | `0/0`, `0/1` or `1/1`: the most likely genotype from `alt_reads` and `ref_reads` at each junction (binomial: 2% of alternative reads for `0/0`, 30% for `0/1`, as alternative reads are under-counted, 98% for `1/1`; with two junctions, their log-likelihoods are averaged). `.` without junction or without any counted read. |
| `gq` | Genotype quality: phred-scaled probability that `genotype` is wrong, equal priors, at most 99. |
| `alt_reads` | Reads clipped at a junction (within 5 bp) by at least `--min-short-clip-len` bases of median quality at least `--min-clip-quality`, summed over the junctions. |
| `ref_reads` | Reads aligned over the whole site (both junctions) with at least `--min-short-clip-len` bases on each side, counted at each junction (twice with two junctions). |

Missing values are written `.`.

## BED (`--format bed`)

| Column | Description |
|---|---|
| 1 | Contig. |
| 2, 3 | The TSD when `tsd_length` is known, otherwise the base right before the insertion. |
| 4 | Candidate identifier. |
| 5 | `support`, capped at 1000. |
| 6 | `strand`. |
| 7 | `confidence` for a `PASS` candidate, otherwise the failed filters. |

## JSON (`--format json`)

An object with `mei_rs` (version), `command`, `sample`, `positions` (the
coordinate convention) and `candidates`. Each candidate has the TSV columns
(`filter` as a list, `both_sides` as a boolean, missing values, such as
the `confidence` of a filtered candidate or the `genotype` of a candidate
without junction, as `null`), except `assembled`,
`tsd` and `polya_length`, which are in `assembly` (for assembled candidates
only, with `tail_base`, and `left` and `right`, each with its 1-based
`junction`, its `contig` in reference orientation and the `inserted` bases
of the contig), and a `signals` list. Each signal has:

| Field | Description |
|---|---|
| `kind` | `soft_clip`, `mate_unmapped`, `mate_low_mapq`, `mate_elsewhere` or `large_insert`. |
| `position` | Junction of a soft-clip, or end of the read facing the insertion for a mate signal. |
| `side` | Side of the read where the inserted sequence lies: `left` or `right`. |
| `read`, `segment`, `strand` | Read name, segment (1 or 2) and strand of the read. |
| `read_start`, `read_end` | Aligned bases of the read, 1-based inclusive. |
| `mapq` | Mapping quality of the read (`null` if unavailable). |
| `clipped_sequence`, `median_quality`, `tail`, `supplementary`, `short` | Soft-clips only: the clipped bases (reference orientation), their median base quality, the poly(A/T) tail at the junction (e.g. `A20`), whether the read has supplementary alignments (`SA` tag), and whether the clip is short (below `--min-clip-len`). |
| `mate_mapq` | `mate_low_mapq` only: mapping quality of the mate. |
| `mate_contig` | `mate_elsewhere` only: contig of the mate. |
| `insert_length` | `large_insert` only: absolute template length. |
| `mei_mate` | Mate signals checked against the consensus sequences only: whether the mate is a mobile element of the candidate family. |
