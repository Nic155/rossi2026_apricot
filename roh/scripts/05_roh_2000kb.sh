#!/bin/bash
#SBATCH --job-name=apricot_roh_2000kb
#SBATCH --output=/home/nrossi/work/ROH/logs/%x_%j.out
#SBATCH --error=/home/nrossi/work/ROH/logs/%x_%j.err
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00

set -euo pipefail

WORKDIR=/home/nrossi/work/ROH
KB=2000
OUT=${WORKDIR}/roh/roh_${KB}kb

module load bioinfo/PLINK/1.90b7

plink \
    --bfile "${WORKDIR}/plink/05_thinned" \
    --homozyg \
    --homozyg-kb ${KB} \
    --homozyg-window-snp 50 \
    --homozyg-snp 50 \
    --homozyg-gap 100 \
    --homozyg-density 20 \
    --homozyg-window-missing 5 \
    --homozyg-window-het 2 \
    --homozyg-window-threshold 0.05 \
    --allow-extra-chr \
    --threads 8 \
    --out "${OUT}"

echo "Step 5 complete: ROH at ${KB}kb threshold written to ${OUT}"
