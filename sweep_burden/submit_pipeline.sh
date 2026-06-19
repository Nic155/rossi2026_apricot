#!/usr/bin/env bash
set -euo pipefail

WORKDIR="$HOME/work/sweep_gl"
cd "$WORKDIR"

mkdir -p scripts slurm logs results results/regions results/variants results/qc

# Load the cluster modules needed by the pipeline before preflight checks.
# shellcheck disable=SC1091
source scripts/load_cluster_modules.sh

command -v sbatch >/dev/null 2>&1 || {
    printf 'ERROR: sbatch was not found in PATH. Run this on a SLURM login node or load SLURM tools.\n' >&2
    exit 1
}

command -v python3 >/dev/null 2>&1 || {
    printf 'ERROR: python3 was not found in PATH. Load a Python module before submitting.\n' >&2
    exit 1
}

command -v bedtools >/dev/null 2>&1 || {
    printf 'ERROR: bedtools was not found in PATH after loading bioinfo/bedtools/2.31.1.\n' >&2
    exit 1
}

job1="$(sbatch --parsable slurm/01_make_sweep_regions.sbatch)"
job2="$(sbatch --parsable --dependency=afterok:"$job1" slurm/02_make_cds_regions.sbatch)"
job3="$(sbatch --parsable --dependency=afterok:"$job2" slurm/03_filter_and_assign_load_variants.sbatch)"
job4="$(sbatch --parsable --dependency=afterok:"$job3" slurm/04_calculate_genotype_burden.sbatch)"

printf 'Submitted SLURM jobs:\n'
printf '  01_make_sweep_regions: %s\n' "$job1"
printf '  02_make_cds_regions: %s\n' "$job2"
printf '  03_filter_and_assign_load_variants: %s\n' "$job3"
printf '  04_calculate_genotype_burden: %s\n' "$job4"
printf '\nFinal expected output:\n'
printf '  results/final_genotype_sweep_control_burden.tsv\n'
