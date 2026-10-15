#!/usr/bin/env bash
# Run every fuzz target for a set time (default 30 s each), starting from a
# few seeds taken from tests/data. Needs a nightly toolchain and cargo-fuzz.
#
#   fuzz/run.sh [SECONDS]
#
# The corpus and any crash go to a temporary directory next to the target
# directory, never into the repository: a crash is printed with its path and
# the input as hex.
set -euo pipefail

cd "$(dirname "$0")/.."
seconds=${1:-30}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

for target in bed_parse consensus_parse type_sequence align_local poly_tail clip_lengths; do
    corpus="$work/corpus/$target"
    mkdir -p "$corpus"
    case $target in
        bed_parse) cp tests/data/targets.bed tests/data/targets_tiled.bed "$corpus"/ ;;
        consensus_parse) head -c 4000 data/consensus/mei_consensus.fa > "$corpus/seed" ;;
        type_sequence) head -c 300 data/consensus/mei_consensus.fa | tail -c 200 > "$corpus/seed" ;;
    esac
    echo "== $target ($seconds s)"
    cargo +nightly fuzz run "$target" "$corpus" -- -max_total_time="$seconds" -artifact_prefix="$work/"
done
