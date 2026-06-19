#!/bin/bash
#SBATCH --job-name=apricot_merge_froh
#SBATCH --output=/home/nrossi/work/ROH/logs/%x_%j.out
#SBATCH --error=/home/nrossi/work/ROH/logs/%x_%j.err
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0:15:00

set -euo pipefail

WORKDIR=/home/nrossi/work/ROH

python3 "${WORKDIR}/scripts/merge_froh.py"

echo "Step 6 complete: F_ROH table at ${WORKDIR}/results/FROH_all_thresholds.tsv"
head -3 "${WORKDIR}/results/FROH_all_thresholds.tsv"
