#!/usr/bin/env bash
# Regenerate the synthetic test fixtures (see README.md).
# Requires python3 and samtools; the output is deterministic for a given
# samtools version.
set -euo pipefail

cd "$(dirname "$0")"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

python3 -I simulate.py "$tmp"
mv "$tmp/reference.fa" "$tmp/targets.bed" "$tmp/targets_tiled.bed" "$tmp/truth.tsv" .
samtools faidx reference.fa

for sample in positive negative; do
    samtools sort --no-PG -O bam -o "$sample.bam" "$tmp/$sample.sam"
    samtools index "$sample.bam"
    samtools view --no-PG -C -T reference.fa -o "$tmp/$sample.cram" "$sample.bam"
    # Drop the UR tags: they hold the absolute path of the reference.
    samtools reheader --no-PG -c 'perl -pe "s/\tUR:[^\t\n]*//g"' "$tmp/$sample.cram" > "$sample.cram"
    samtools index "$sample.cram"
done
