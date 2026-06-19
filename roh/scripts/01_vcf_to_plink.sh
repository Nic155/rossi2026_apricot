#!/bin/bash
#SBATCH --job-name=apricot_vcf2plink
#SBATCH --output=/home/nrossi/work/ROH/logs/%x_%j.out
#SBATCH --error=/home/nrossi/work/ROH/logs/%x_%j.err
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00

set -euo pipefail

WORKDIR=/home/nrossi/work/ROH
VCF=/home/nrossi/work/genetic_load/parmeniaca/1-RawVCF/apricot_collection_2019_marouch_v3.1.snps.vcf.gz
OUT=${WORKDIR}/plink/01_raw

module load bioinfo/Bcftools/1.21
module load bioinfo/PLINK/1.90b7

# Extract biallelic SNPs only, then pipe directly into PLINK.
# bcftools -m2 -M2 -v snps: keep sites with exactly 2 alleles and type=SNP.
# PLINK reads VCF from stdin via /dev/stdin (single-pass, compatible with pipes).
bcftools view -m2 -M2 -v snps "${VCF}" | \
    plink \
        --vcf /dev/stdin \
        --allow-extra-chr \
        --chr-set 8 \
        --double-id \
        --make-bed \
        --threads 8 \
        --out "${OUT}"

echo "Step 1 complete: biallelic SNPs converted to PLINK binary at ${OUT}"
