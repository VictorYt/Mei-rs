#!/usr/bin/env bash
# Check the VCF written by `mei-rs scan --format vcf` with the tools users run
# on it (see tests/data/README.md): bcftools, tabix and the VCF validator.
#
# For each fixture (BAM and CRAM), the VCF is written plain and bgzipped
# (`.vcf.gz` + `.tbi`), then:
#   * bcftools reads both without a warning and finds the same records;
#   * the records are well-formed for the VCF validator (vcf-validator, from
#     vcftools);
#   * the REF alleles match the reference (`bcftools norm --check-ref e`);
#   * the `.tbi` we write returns the same records as the one tabix builds,
#     for every record position and for each contig.
#
# Requires bcftools, tabix and vcf-validator. The binary is MEI_RS if set
# (default: `cargo run --release`).
set -euo pipefail

cd "$(dirname "$0")/.."
data=tests/data
if [[ -n "${MEI_RS:-}" ]]; then
    mei=("$MEI_RS")
else
    cargo build --release --locked -q
    mei=("${CARGO_TARGET_DIR:-target}/release/mei-rs")
fi
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
fail=0
bad() { echo "FAIL: $*" >&2; fail=1; }

for tool in bcftools tabix vcf-validator; do
    command -v "$tool" >/dev/null || { echo "missing: $tool" >&2; exit 2; }
done

check() {
    local name=$1 alignment=$2 targets=$3 base="$tmp/$1.$4"
    local args=(scan --bam "$alignment" --bed "$targets" --reference $data/reference.fa --keep-filtered --format vcf)
    "${mei[@]}" "${args[@]}" -o "$base.vcf" 2> "$tmp/log"
    "${mei[@]}" "${args[@]}" -o "$base.vcf.gz" 2> "$tmp/log"
    local n; n=$(grep -vc '^#' "$base.vcf" || true)
    echo "== $name.$4: $n records"

    # bcftools: no warning, same records plain and bgzipped.
    bcftools view "$base.vcf" > "$base.view.vcf" 2> "$tmp/err"
    [[ -s "$tmp/err" ]] && bad "$name.$4: bcftools view (plain) warns: $(head -3 "$tmp/err")"
    bcftools view "$base.vcf.gz" > "$base.view.gz.vcf" 2> "$tmp/err"
    [[ -s "$tmp/err" ]] && bad "$name.$4: bcftools view (bgzipped) warns: $(head -3 "$tmp/err")"
    diff -q <(grep -v '^##' "$base.view.vcf") <(grep -v '^##' "$base.view.gz.vcf") >/dev/null \
        || bad "$name.$4: plain and bgzipped records differ"
    bcftools index -n "$base.vcf.gz" > "$tmp/n" 2> "$tmp/err" \
        || bad "$name.$4: bcftools cannot use our .tbi: $(head -3 "$tmp/err")"
    [[ "$(cat "$tmp/n" 2>/dev/null)" == "$n" ]] || bad "$name.$4: .tbi counts $(cat "$tmp/n") records, file has $n"

    # VCF validator (header declarations, field types, ordering).
    if ! vcf-validator "$base.vcf" > "$tmp/validator" 2>&1; then
        bad "$name.$4: vcf-validator: $(head -5 "$tmp/validator")"
    fi

    # REF alleles against the reference.
    bcftools norm --check-ref e -f $data/reference.fa "$base.vcf" -o /dev/null 2> "$tmp/err" \
        || bad "$name.$4: REF differs from the reference: $(head -3 "$tmp/err")"

    # Our .tbi against the one tabix builds, on every record and contig.
    cp "$base.vcf.gz" "$base.ref.vcf.gz"
    tabix -p vcf -f "$base.ref.vcf.gz"
    while read -r contig pos; do
        for region in "$contig:$pos-$pos" "$contig:$((pos > 5 ? pos - 5 : 1))-$((pos + 5))" "$contig"; do
            diff -q <(tabix "$base.vcf.gz" "$region") <(tabix "$base.ref.vcf.gz" "$region") >/dev/null \
                || bad "$name.$4: .tbi differs from tabix on $region"
        done
    done < <(grep -v '^#' "$base.vcf" | cut -f1,2)
}

for sample in positive negative homozygous tailonly; do
    for ext in bam cram; do
        check "$sample" "$data/$sample.$ext" $data/targets.bed "$ext"
    done
done
exit $fail
