#!/usr/bin/env bash
set -euo pipefail

WORKDIR="${WORKDIR:-$HOME/work/sweep_gl}"
GFF="${GFF:-$HOME/work/genetic_load/parmeniaca/5-ReferenceData/Prunus_armeniaca_Marouch_n14/genes.gff}"
SWEEP_BED="${SWEEP_BED:-results/regions/sweep_regions_merged.bed}"

CDS_RAW="results/regions/cds_raw.bed"
CDS_REDUCED="results/regions/cds_reduced.bed"
SWEEP_CDS="results/regions/sweep_cds.bed"
CONTROL_CDS="results/regions/control_cds.bed"
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
[[ -s "$SWEEP_BED" ]] || fail "Merged sweep BED does not exist or is empty: $SWEEP_BED"

TMP_CDS_RAW="results/regions/cds_raw.unsorted.tmp.bed"
TMP_SWEEP_CDS="results/regions/sweep_cds.unsorted.tmp.bed"
TMP_SWEEP_CDS_SORTED="results/regions/sweep_cds.sorted.tmp.bed"
TMP_CONTROL_CDS="results/regions/control_cds.unsorted.tmp.bed"
trap 'rm -f "$TMP_CDS_RAW" "$TMP_SWEEP_CDS" "$TMP_SWEEP_CDS_SORTED" "$TMP_CONTROL_CDS"' EXIT

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

printf 'Intersecting reduced CDS with merged sweep regions\n'
bedtools intersect -a "$CDS_REDUCED" -b "$SWEEP_BED" > "$TMP_SWEEP_CDS"
sort -k1,1 -k2,2n "$TMP_SWEEP_CDS" > "$TMP_SWEEP_CDS_SORTED"
bedtools merge -i "$TMP_SWEEP_CDS_SORTED" > "$SWEEP_CDS"

printf 'Creating control CDS as reduced CDS minus sweep CDS\n'
bedtools subtract -a "$CDS_REDUCED" -b "$SWEEP_CDS" > "$TMP_CONTROL_CDS"
sort -k1,1 -k2,2n "$TMP_CONTROL_CDS" > "$CONTROL_CDS"

total_reduced_CDS_bp="$(bed_length "$CDS_REDUCED")"
sweep_CDS_bp="$(bed_length "$SWEEP_CDS")"
control_CDS_bp="$(bed_length "$CONTROL_CDS")"

[[ "$sweep_CDS_bp" -gt 0 ]] || fail "sweep_CDS_bp is zero. Check sweep/CDS overlap."
[[ "$control_CDS_bp" -gt 0 ]] || fail "control_CDS_bp is zero. Check sweep regions and CDS annotation."

{
    printf 'region_class\tcds_length_bp\n'
    printf 'all_CDS\t%s\n' "$total_reduced_CDS_bp"
    printf 'sweep\t%s\n' "$sweep_CDS_bp"
    printf 'control\t%s\n' "$control_CDS_bp"
} > "$QC_LENGTHS"

printf 'total_reduced_CDS_bp: %s\n' "$total_reduced_CDS_bp"
printf 'sweep_CDS_bp: %s\n' "$sweep_CDS_bp"
printf 'control_CDS_bp: %s\n' "$control_CDS_bp"
printf 'Wrote %s\n' "$CDS_RAW"
printf 'Wrote %s\n' "$CDS_REDUCED"
printf 'Wrote %s\n' "$SWEEP_CDS"
printf 'Wrote %s\n' "$CONTROL_CDS"
printf 'Wrote %s\n' "$QC_LENGTHS"
