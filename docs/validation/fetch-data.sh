#!/bin/bash
# Download the public inputs of the evaluation (see MANIFEST.md): about 160 GB (36 exomes).
#
#   fetch-data.sh DIR
#
# Layout written in DIR, as read by `run.sh`:
#   truth/   HGSVC3 MEI calls, UCSC hg38ToHg19 chain
#   ref/     GRCh38 and hs37d5 reference genomes
#   1000g/   SAMPLE.exome.cram(.crai) of the 1000 Genomes exomes of split.tsv
#   giab/    HG002.exome.bam(.bai)
#   targets/ exome targets (build-targets.sh)
# Files already present with the right md5 are not downloaded again.
#
# SET selects the samples of split.tsv: `all` (default, 36 exomes, about 140 GB),
# `tuning` (18) or `test` (18). SAMPLES (a list) overrides it.
set -euo pipefail
dir=${1:?usage: fetch-data.sh DIR}
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$dir"/{truth,ref,1000g,giab,targets}

md5() { if command -v md5sum >/dev/null; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi; }
# get URL FILE [MD5]
get() {
  if [ -f "$2" ] && { [ -z "${3:-}" ] || [ "$(md5 "$2")" = "$3" ]; }; then return; fi
  curl -fL --retry 3 -o "$2" "$1"
  if [ -n "${3:-}" ] && [ "$(md5 "$2")" != "$3" ]; then echo "$2: unexpected md5" >&2; exit 1; fi
}

HGSVC3=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Mobile_Elements/1.0
get $HGSVC3/MEI_Callset_GRCh38.ALL.20241211.csv.gz "$dir/truth/MEI_Callset_GRCh38.ALL.20241211.csv.gz" 0131ecb529c9225dab677aced8d72b0a
get https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz "$dir/truth/hg38ToHg19.over.chain.gz" ff3031d93792f4cbb86af44055efd903

REF=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference
get $REF/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa "$dir/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa" 64b32de2fc934679c16e83a2bc072064
get $REF/phase2_reference_assembly_sequence/hs37d5.fa.gz "$dir/ref/hs37d5.fa.gz" d10eebe06c0dbbcb04253e3294d63efc

# The 1000 Genomes exome index lists the path and the md5 of every file.
INDEX=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/1000genomes.exome.GRCh38DH.alignment.index
get $INDEX "$dir/1000g/alignment.index"
set=${SET:-all}
samples=${SAMPLES:-$(awk -F'\t' -v set="$set" 'NR > 1 && (set == "all" || $2 == set) {print $1}' "$here/split.tsv")}
for s in $samples; do
  grep -F "/$s/" "$dir/1000g/alignment.index" | tr '\t' '\n' | paste - - | grep -E "\.exome\.cram(\.crai)?\s" |
  while read -r path sum; do
    get "https://${path#ftp:/}" "$dir/1000g/$s.exome.${path##*.exome.}" "$sum"
  done
done

# HG002 (NA24385): GIAB exome of Oslo University Hospital, aligned to GRCh37.
GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome
name=151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup
get $GIAB/$name.bam "$dir/giab/HG002.exome.bam" c80f0cab24bfaa504393457b8f7191fa
get $GIAB/$name.bai "$dir/giab/HG002.exome.bam.bai"

"$here/build-targets.sh" "$dir/targets"
