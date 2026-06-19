#!/usr/bin/env bash
#SBATCH --job-name=04_calc_genotype_burden
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

set -euo pipefail

cd "$HOME/work/ROH/roh_gl"
# shellcheck disable=SC1091
source scripts/load_cluster_modules.sh
python3 scripts/04_calculate_genotype_burden.py
