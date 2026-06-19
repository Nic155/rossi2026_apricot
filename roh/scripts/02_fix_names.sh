#!/bin/bash
#SBATCH --job-name=apricot_fixnames
#SBATCH --output=/home/nrossi/work/ROH/logs/%x_%j.out
#SBATCH --error=/home/nrossi/work/ROH/logs/%x_%j.err
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0:15:00

set -euo pipefail

FAM=/home/nrossi/work/ROH/plink/01_raw.fam

# Strip _marouch_v3 suffix from FID (col 1) and IID (col 2)
awk '{
    sub(/_marouch_v3$/, "", $1)
    sub(/_marouch_v3$/, "", $2)
    print $1, $2, $3, $4, $5, $6
}' "${FAM}" > "${FAM}.tmp" && mv "${FAM}.tmp" "${FAM}"

echo "Step 2 complete: _marouch_v3 suffix stripped from ${FAM}"
echo "First 5 lines:"
head -5 "${FAM}"
