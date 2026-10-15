#!/bin/bash
# Run `mei-rs scan` on the evaluation exomes, evaluate the calls against the
# HGSVC3 MEI calls and refit the quality score (see MANIFEST.md).
#
#   run.sh DATA OUT
#
# DATA is the directory written by `fetch-data.sh` (or any directory with the
# same layout); OUT receives the calls, the evaluation and the refit. Nothing
# is written into the repository.
#
# Environment: MEI_RS (binary, default target/release/mei-rs), THREADS,
# SAMPLES (default: the seven 1000 Genomes exomes), SKIP_HG002=1 to skip HG002,
# ALIGNMENT_EXT (`cram`, or `bam` if you converted the exomes).
set -euo pipefail
data=${1:?usage: run.sh DATA OUT}
out=${2:?usage: run.sh DATA OUT}
here=$(cd "$(dirname "$0")" && pwd)
mei_rs=${MEI_RS:-"$here/../../target/release/mei-rs"}
ext=${ALIGNMENT_EXT:-cram}
samples=${SAMPLES:-"HG00512 HG00513 HG00731 HG00733 NA19238 NA19239 NA19240"}
threads=${THREADS:-$(getconf _NPROCESSORS_ONLN)}
truth=$data/truth/MEI_Callset_GRCh38.ALL.20241211.csv.gz
mkdir -p "$out"
"$mei_rs" --version | tee "$out/version.txt"

for s in $samples; do
  start=$(date +%s)
  "$mei_rs" scan --bam "$data/1000g/$s.exome.$ext" \
    --bed "$data/targets/exome_targets.hg38.bed" \
    --reference "$data/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
    --threads "$threads" --keep-filtered -o "$out/$s.tsv"
  echo "$s wall_seconds $(($(date +%s) - start))" >> "$out/times.txt"
done
if [ -z "${SKIP_HG002:-}" ]; then
  "$mei_rs" scan --bam "$data/giab/HG002.exome.bam" \
    --bed "$data/targets/exome_targets.GRCh37.bed" \
    --reference "$data/ref/hs37d5.fa.gz" \
    --threads "$threads" --keep-filtered -o "$out/HG002.tsv"
fi

runs=(); for s in $samples; do runs+=("$s=$out/$s.tsv"); done
python3 -I "$here/evaluate.py" --truth "$truth" --bed "$data/targets/exome_targets.hg38.bed" \
  --min-qual 0,2,5,10,15 --fp-out "$out/fp.tsv" --tp-out "$out/tp.tsv" --fn-out "$out/fn.tsv" \
  "${runs[@]}" | tee "$out/evaluation.txt"

held_out=()
if [ -z "${SKIP_HG002:-}" ]; then
  python3 -I "$here/evaluate.py" --truth "$truth" --bed "$data/targets/exome_targets.GRCh37.bed" \
    --chain "$data/truth/hg38ToHg19.over.chain.gz" --min-qual 0,2,5 \
    --fp-out "$out/fp_hg002.tsv" --tp-out "$out/tp_hg002.tsv" \
    NA24385="$out/HG002.tsv" | tee "$out/evaluation_hg002.txt"
  held_out=(--held-out "NA24385=$out/HG002.tsv" --held-out-fp "$out/fp_hg002.tsv" --held-out-tp "$out/tp_hg002.tsv")
fi

python3 -I "$here/refit.py" --fp "$out/fp.tsv" --tp "$out/tp.tsv" "${held_out[@]}" "${runs[@]}" \
  | tee "$out/refit.txt"
