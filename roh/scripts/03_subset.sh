#!/bin/bash
#SBATCH --job-name=apricot_subset
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
KEEP=${WORKDIR}/plink/keep.txt
OUT=${WORKDIR}/plink/03_subset

module load bioinfo/PLINK/1.90b7

# Build keep file: skip header, FID=IID=sample_name (no _marouch_v3 suffix)
awk 'NR>1 {print $1, $1}' "${WORKDIR}/list.txt" > "${KEEP}"
echo "Keep file created with $(wc -l < "${KEEP}") samples."

plink \
    --bfile "${WORKDIR}/plink/01_raw" \
    --allow-extra-chr \
    --keep "${KEEP}" \
    --make-bed \
    --threads 8 \
    --out "${OUT}"

echo "Step 3 complete: subsetted to $(wc -l < "${OUT}.fam") samples at ${OUT}"
