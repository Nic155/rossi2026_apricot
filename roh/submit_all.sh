#!/bin/bash
set -euo pipefail

WORKDIR=/home/nrossi/work/ROH
SCRIPTS=${WORKDIR}/scripts

JOB1=$(sbatch --parsable "${SCRIPTS}/01_vcf_to_plink.sh")
JOB2=$(sbatch --parsable --dependency=afterok:${JOB1} "${SCRIPTS}/02_fix_names.sh")
JOB3=$(sbatch --parsable --dependency=afterok:${JOB2} "${SCRIPTS}/03_subset.sh")
JOB4=$(sbatch --parsable --dependency=afterok:${JOB3} "${SCRIPTS}/04_filter.sh")
JOB4B=$(sbatch --parsable --dependency=afterok:${JOB4} "${SCRIPTS}/04b_thin.sh")
JOB5A=$(sbatch --parsable --dependency=afterok:${JOB4B} "${SCRIPTS}/05_roh_100kb.sh")
JOB5B=$(sbatch --parsable --dependency=afterok:${JOB4B} "${SCRIPTS}/05_roh_500kb.sh")
JOB5C=$(sbatch --parsable --dependency=afterok:${JOB4B} "${SCRIPTS}/05_roh_800kb.sh")
JOB5D=$(sbatch --parsable --dependency=afterok:${JOB4B} "${SCRIPTS}/05_roh_1000kb.sh")
JOB5E=$(sbatch --parsable --dependency=afterok:${JOB4B} "${SCRIPTS}/05_roh_2000kb.sh")
ROH_DEP="${JOB5A}:${JOB5B}:${JOB5C}:${JOB5D}:${JOB5E}"
JOB6=$(sbatch --parsable --dependency=afterok:${ROH_DEP} "${SCRIPTS}/06_merge_froh.sh")
