#!/usr/bin/env bash
set -euo pipefail

WORKDIR="${WORKDIR:-$HOME/work/ROH/roh_gl}"
GFF="${GFF:-$HOME/work/genetic_load/parmeniaca/5-ReferenceData/Prunus_armeniaca_Marouch_n14/genes.gff}"

CDS_RAW="results/regions/cds_raw.bed"
CDS_REDUCED="results/regions/cds_reduced.bed"
QC_LENGTHS="results/qc/cds_lengths.tsv"

fail() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

bed_length() {
    awk 'BEGIN {sum = 0} NF >= 3 {sum += $3 - $2} END {print sum + 0}' "$1"
}

cd "$WORKDIR"
mkdir -p results/regions results/qc

command -v bedtools >/dev/null 2>&1 || fail "bedtools was not found in PATH. Load bedtools before running this step."
[[ -s "$GFF" ]] || fail "GFF file does not exist or is empty: $GFF"

TMP_CDS_RAW="results/regions/cds_raw.unsorted.tmp.bed"
trap 'rm -f "$TMP_CDS_RAW"' EXIT

printf 'Extracting CDS features from %s\n' "$GFF"
awk 'BEGIN {OFS = "\t"}
     $0 !~ /^#/ && $3 == "CDS" {
         if ($4 < 1 || $5 < $4) {
             printf("Invalid CDS coordinates at line %d\n", NR) > "/dev/stderr";
             exit 1;
         }
         print $1, $4 - 1, $5;
     }' "$GFF" > "$TMP_CDS_RAW"

[[ -s "$TMP_CDS_RAW" ]] || fail "No CDS features were extracted from the GFF."

sort -k1,1 -k2,2n "$TMP_CDS_RAW" > "$CDS_RAW"
bedtools merge -i "$CDS_RAW" > "$CDS_REDUCED"

[[ -s "$CDS_REDUCED" ]] || fail "Reduced CDS BED is empty after bedtools merge."

total_CDS_bp="$(bed_length "$CDS_REDUCED")"

{
    printf 'region_class\tcds_length_bp\n'
    printf 'all_CDS\t%s\n' "$total_CDS_bp"
} > "$QC_LENGTHS"

printf 'total_CDS_bp: %s\n' "$total_CDS_bp"
printf 'Wrote %s\n' "$CDS_RAW"
printf 'Wrote %s\n' "$CDS_REDUCED"
printf 'Wrote %s\n' "$QC_LENGTHS"
