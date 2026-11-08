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
# The samples are those of split.tsv: the tuning samples are evaluated and
# used to refit the quality score, the test samples are only evaluated (with
# the weights fitted on the tuning samples), and the `related` samples (left
# out of both sets because a relative is a tuning sample, #111) are evaluated
# apart. `check_split.py` stops the run when a family spans the tuning and the
# test sets. Every sample is scanned with
# `--profile discovery` (no quality threshold) so that the evaluation sees all
# the `PASS` calls and applies the `--min-qual` thresholds itself.
#
# Environment: MEI_RS (binary, default target/release/mei-rs), THREADS
# (default 4; the scans are run one after the other), NICE (default 19),
# SET (`all`, `tuning`, `test` or `related`; default `all`), SAMPLES, TEST_SAMPLES and
# RELATED_SAMPLES (lists replacing the samples of each set of split.tsv, e.g. SAMPLES="HG00512
# HG00513 HG00731 HG00733 NA19238 NA19239 NA19240" for the seven exomes of the
# score weights), SKIP_HG002=1 to skip HG002, ALIGNMENT_EXT (`cram`, or `bam` if you
# converted the exomes).
set -euo pipefail
data=${1:?usage: run.sh DATA OUT}
out=${2:?usage: run.sh DATA OUT}
here=$(cd "$(dirname "$0")" && pwd)
mei_rs=${MEI_RS:-"$here/../../target/release/mei-rs"}
ext=${ALIGNMENT_EXT:-cram}
threads=${THREADS:-4}
nice_level=${NICE:-19}
set=${SET:-all}
split() { awk -F'\t' -v want="$1" 'NR > 1 && $2 == want {print $1}' "$here/split.tsv"; }
python3 -I "$here/check_split.py" "$here/split.tsv"
tuning=""; test_samples=""; related=""
case $set in all | tuning) tuning=${SAMPLES:-$(split tuning)} ;; esac
case $set in all | test) test_samples=${TEST_SAMPLES:-$(split test)} ;; esac
case $set in all | related) related=${RELATED_SAMPLES:-$(split related)} ;; esac
truth=$data/truth/MEI_Callset_GRCh38.ALL.20241211.csv.gz
mkdir -p "$out"
"$mei_rs" --version | tee "$out/version.txt"

for s in $tuning $test_samples $related; do
  [ -s "$out/$s.tsv" ] && continue
  start=$(date +%s)
  nice -n "$nice_level" "$mei_rs" scan --bam "$data/1000g/$s.exome.$ext" \
    --bed "$data/targets/exome_targets.hg38.bed" \
    --reference "$data/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
    --threads "$threads" --profile discovery --keep-filtered -o "$out/$s.tsv.part"
  mv "$out/$s.tsv.part" "$out/$s.tsv"
  echo "$s wall_seconds $(($(date +%s) - start))" >> "$out/times.txt"
done
if [ -z "${SKIP_HG002:-}" ]; then
  [ -s "$out/HG002.tsv" ] || {
    nice -n "$nice_level" "$mei_rs" scan --bam "$data/giab/HG002.exome.bam" \
      --bed "$data/targets/exome_targets.GRCh37.bed" \
      --reference "$data/ref/hs37d5.fa.gz" \
      --threads "$threads" --profile discovery --keep-filtered -o "$out/HG002.tsv.part"
    mv "$out/HG002.tsv.part" "$out/HG002.tsv"
  }
fi

evaluate() { # NAME SAMPLES...
  local name=$1; shift
  local runs=(); for s in "$@"; do runs+=("$s=$out/$s.tsv"); done
  python3 -I "$here/evaluate.py" --truth "$truth" --bed "$data/targets/exome_targets.hg38.bed" \
    --min-qual 0,2,5,10,15 --fp-out "$out/fp_$name.tsv" --tp-out "$out/tp_$name.tsv" \
    --fn-out "$out/fn_$name.tsv" "${runs[@]}" | tee "$out/evaluation_$name.txt"
}
[ -n "$tuning" ] && evaluate tuning $tuning
[ -n "$test_samples" ] && evaluate test $test_samples
[ -n "$related" ] && evaluate related $related

held_out=()
if [ -z "${SKIP_HG002:-}" ]; then
  python3 -I "$here/evaluate.py" --truth "$truth" --bed "$data/targets/exome_targets.GRCh37.bed" \
    --chain "$data/truth/hg38ToHg19.over.chain.gz" --min-qual 0,2,5 \
    --fp-out "$out/fp_hg002.tsv" --tp-out "$out/tp_hg002.tsv" \
    NA24385="$out/HG002.tsv" | tee "$out/evaluation_hg002.txt"
  held_out=(--held-out "NA24385=$out/HG002.tsv" --held-out-fp "$out/fp_hg002.tsv" --held-out-tp "$out/tp_hg002.tsv")
fi

# Refit on the tuning samples; the test samples are scored with those weights.
tests=(); for s in $test_samples; do tests+=(--test "$s=$out/$s.tsv"); done
test_labels=(); [ -n "$test_samples" ] && test_labels=(--test-fp "$out/fp_test.tsv" --test-tp "$out/tp_test.tsv")
runs=(); for s in $tuning; do runs+=("$s=$out/$s.tsv"); done
python3 -I "$here/refit.py" --fp "$out/fp_tuning.tsv" --tp "$out/tp_tuning.tsv" \
  "${tests[@]}" "${test_labels[@]}" "${held_out[@]}" "${runs[@]}" | tee "$out/refit.txt"
