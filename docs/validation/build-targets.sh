#!/bin/bash
# Rebuild the exome capture targets of the evaluation (see MANIFEST.md).
#
#   build-targets.sh DIR
#
# Downloads the 1000 Genomes exome pull-down targets (GRCh37) and the UCSC
# `hg19ToHg38` chain into DIR, lifts the targets to GRCh38 with CrossMap (the
# targets mapped in one piece are kept; the others are listed in
# `targets.hg38.raw.bed.unmap`) and writes `exome_targets.GRCh37.bed` (contigs
# `1`..`22`, as the GRCh37 exomes) and `exome_targets.hg38.bed` (`chr1`..).
# Needs curl, python3 (venv, pip) and md5sum or md5.
set -euo pipefail
dir=${1:?usage: build-targets.sh DIR}
mkdir -p "$dir"
cd "$dir"

md5() { if command -v md5sum >/dev/null; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi; }
check() {
  [ "$(md5 "$1")" = "$2" ] || { echo "$1: unexpected md5 $(md5 "$1") (expected $2)" >&2; exit 1; }
}

curl -sfL -o targets.GRCh37.chr.bed \
  https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/exome_pull_down_targets/20130108.exome.targets.bed
curl -sfL -o hg19ToHg38.over.chain.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
check targets.GRCh37.chr.bed 3285b0ba60e45550297d5e072d105483
check hg19ToHg38.over.chain.gz 35887f73fe5e2231656504d1f6430900

python3 -m venv venv
./venv/bin/pip install -q CrossMap==0.7.3
./venv/bin/CrossMap bed hg19ToHg38.over.chain.gz targets.GRCh37.chr.bed targets.hg38.raw.bed

awk -F'\t' '{print $1, $2, $3}' OFS='\t' targets.hg38.raw.bed | sort -k1,1 -k2,2n > exome_targets.hg38.bed
sed 's/^chr//' targets.GRCh37.chr.bed > exome_targets.GRCh37.bed
check exome_targets.hg38.bed f4f589406344a74b624e1af7d8ff8e33
check exome_targets.GRCh37.bed 14805f7965dc8b44772b410d1141393c
echo "targets: $(wc -l < exome_targets.GRCh37.bed) GRCh37, $(wc -l < exome_targets.hg38.bed) GRCh38"
