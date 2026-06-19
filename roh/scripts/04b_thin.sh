#!/bin/bash
#SBATCH --job-name=apricot_thin
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
OUT=${WORKDIR}/plink/05_thinned

module load bioinfo/PLINK/1.90b7

# Thin to 1 SNP per 10 kb and apply MAF filter.
# --bp-space 10000: retain at most one SNP per 10 000 bp window.
# --maf 0.05: safe to apply post-thinning because ~18 k SNPs remain;
#   at this density MAF filtering no longer disrupts ROH signal.
# --geno 0.05: redundant after step 4 but re-applied for safety.
plink \
    --bfile "${WORKDIR}/plink/04_filtered" \
    --allow-extra-chr \
    --maf 0.05 \
    --geno 0.05 \
    --bp-space 10000 \
    --make-bed \
    --threads 8 \
    --out "${OUT}"

echo "Step 4b complete: thinned dataset at ${OUT}"
grep -E "^[0-9]+ variants" "${OUT}.log" || true
