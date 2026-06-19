#!/bin/bash
#SBATCH --job-name=apricot_filter
#SBATCH --output=/home/nrossi/work/ROH/logs/%x_%j.out
#SBATCH --error=/home/nrossi/work/ROH/logs/%x_%j.err
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00

set -euo pipefail

WORKDIR=/home/nrossi/work/ROH
OUT=${WORKDIR}/plink/04_filtered

module load bioinfo/PLINK/1.90b7

plink \
    --bfile "${WORKDIR}/plink/03_subset" \
    --allow-extra-chr \
    --geno 0.05 \
    --make-bed \
    --threads 8 \
    --out "${OUT}"

echo "Step 4 complete: filtered dataset at ${OUT}"
grep -E "^[0-9]+ variants" "${OUT}.log" || true
