# Output of `mei-rs scan`

`mei-rs scan` writes one line (or JSON object) per candidate insertion site,
sorted by contig (alignment header order) and position. By default only the
candidates that pass every filter are written; `--keep-filtered` also writes
those flagged `no_junction`, `no_tail`, `probe_edge` or `low_ratio`. Clusters with fewer than
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
| `left_junction` | End of the upstream flank: most supported junction of the reads clipped on their right end (`.` if none). |
| `right_junction` | Start of the downstream flank: most supported junction of the reads clipped on their left end (`.` if none). |
| `tsd_length` | `left_junction - right_junction` when each junction has at least `--min-clips` soft-clips and the flanks overlap by 1 to 50 bp: the length of the target site duplication (`.` otherwise). |
| `strand` | Strand of the inserted element: `+` when poly(A) tails outnumber poly(T) tails in the soft-clips, `-` for the reverse, `.` on ties. |
| `filter` | `PASS`, or the failed filters separated by `;`: `no_junction` (fewer than `--min-clips` soft-clips at the junctions: mate signals alone do not resolve the insertion to the base), `no_tail` (no soft-clip with a poly(A/T) tail; disabled by `--allow-no-tail`), `probe_edge` (every soft-clip stops within `--probe-edge-tolerance` of a probe start or end, none has a poly(A/T) tail and there is no TSD: capture artefact), `low_ratio` (`signal_ratio` below `--min-signal-ratio`). |
| `confidence` | For `PASS` candidates: `high` (a poly(A/T) tail and a TSD: complete evidence of a retrotransposition), `medium` (mate signals, without a TSD or, with `--allow-no-tail`, without a tail), `low` (soft-clips only, without a TSD); `.` for filtered candidates. |
| `support` | Distinct supporting reads (by name). |
| `left_clips` | Reads soft-clipped on their left end (the insertion lies before them). |
| `right_clips` | Reads soft-clipped on their right end (the insertion lies after them). |
| `left_junction_clips` | Soft-clips within 2 bp of `left_junction`. |
| `right_junction_clips` | Soft-clips within 2 bp of `right_junction`. |
| `unmapped_mates` | Reads whose mate is unmapped. |
| `discordant_mates` | Reads whose mate maps with a low MAPQ, to another contig, or far away on the same contig. |
| `polya_clips` | Soft-clips with a poly(A) tail at the junction. |
| `polyt_clips` | Soft-clips with a poly(T) tail at the junction. |
| `both_sides` | `yes` if signals come from both flanks of the insertion. |
| `distinct_starts` | Distinct fragment start positions (strand-aware) of the supporting reads; far below `support`, it points to PCR duplicates. |
| `depth` | Reads (after filtering) covering one of the two bases around `position`. |
| `signal_ratio` | `support / depth` (3 decimals). Mate signals come from reads that do not cover `position`, so it can exceed 1. |
| `start`, `end` | Positions of the first and last signal of the cluster. |

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
the `confidence` of a filtered candidate, as `null`) and
a `signals` list. Each signal has:

| Field | Description |
|---|---|
| `kind` | `soft_clip`, `mate_unmapped`, `mate_low_mapq`, `mate_elsewhere` or `large_insert`. |
| `position` | Junction of a soft-clip, or end of the read facing the insertion for a mate signal. |
| `side` | Side of the read where the inserted sequence lies: `left` or `right`. |
| `read`, `segment`, `strand` | Read name, segment (1 or 2) and strand of the read. |
| `read_start`, `read_end` | Aligned bases of the read, 1-based inclusive. |
| `mapq` | Mapping quality of the read (`null` if unavailable). |
| `clipped_sequence`, `median_quality`, `tail`, `supplementary` | Soft-clips only: the clipped bases (reference orientation), their median base quality, the poly(A/T) tail at the junction (e.g. `A20`), and whether the read has supplementary alignments (`SA` tag). |
| `mate_mapq` | `mate_low_mapq` only: mapping quality of the mate. |
| `mate_contig` | `mate_elsewhere` only: contig of the mate. |
| `insert_length` | `large_insert` only: absolute template length. |
