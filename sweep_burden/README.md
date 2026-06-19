# Sweep-Control Deleterious Burden Pipeline

## Purpose

This SLURM-based command-line pipeline calculates deleterious genetic-load burden for every genotype by comparing candidate C1 sweep CDS regions against control CDS regions. The final burden is reported in the soybean-paper style:

```text
deleterious genotype dosage per 100 kb CDS
```

The final output has one row per genotype and exactly three columns:

```text
genotype
sweep_deleterious_dosage_per_100kb_CDS
control_deleterious_dosage_per_100kb_CDS
```

## Input Files

Expected paths:

```text
~/work/sweep_gl/pi_ratio_W2_C1.csv
~/work/sweep_gl/result3-provean_forGL+score+class.C1_W2.subset.tsv
~/work/genetic_load/parmeniaca/5-ReferenceData/Prunus_armeniaca_Marouch_n14/genes.gff
```

## Sweep-Region Definition

Sweep windows are read from `pi_ratio_W2_C1.csv`, a semicolon-separated table of 10 kb pi-ratio windows.

The default threshold-indicator column is:

```text
95% threshold
```

A window is selected when this threshold column is non-empty and not `NA`. Chromosomes are converted from values like `1` to BED chromosomes like `chr1`.

Selected windows are converted from 1-based inclusive coordinates to 0-based half-open BED coordinates:

```text
BED start = pos_start - 1
BED end   = pos_end
```

Selected 10 kb windows are sorted and merged with:

```text
bedtools merge -d 10000
```

This allows one missing 10 kb bin between selected windows. The merged sweep regions are written to:

```text
results/regions/sweep_regions_merged.bed
```

To change the threshold column, edit the argument in:

```text
slurm/01_make_sweep_regions.sbatch
```

or run:

```text
python3 scripts/01_make_sweep_regions.py --threshold-col "Top 1% threshold"
```

## Control-Region Definition

Control regions are reduced CDS bases outside sweep CDS, not the whole genome.

The pipeline creates:

```text
sweep_CDS   = reduced_CDS intersect merged sweep regions
control_CDS = reduced_CDS minus sweep_CDS
```

## Deleterious/Load Variant Definition

Load variants are defined only by `IMPACT`:

```text
IMPACT == "PROVEAN" OR IMPACT == "HIGH"
```

`PROVEAN_SCORE` is not used for filtering.

`DEL_CLASS` is not used for filtering.

## CDS Length Definition

CDS features are extracted from `genes.gff` using only rows where column 3 is `CDS`.

GFF CDS intervals are converted to BED:

```text
BED start = GFF start - 1
BED end   = GFF end
```

CDS intervals are sorted and merged with `bedtools merge` before length calculation. This prevents double-counting overlapping CDS intervals from multiple transcripts.

The CDS length table is:

```text
results/qc/cds_lengths.tsv
```

with rows for:

```text
all_CDS
sweep
control
```

## Formula

For every genotype:

```text
sweep_deleterious_dosage_per_100kb_CDS =
    sum genotype dosage for PROVEAN/HIGH variants in sweep CDS / sweep_CDS_bp * 100000

control_deleterious_dosage_per_100kb_CDS =
    sum genotype dosage for PROVEAN/HIGH variants in control CDS / control_CDS_bp * 100000
```

NA genotype calls are ignored in dosage sums. NA is not treated as zero. The denominator is CDS length, not the number of called genotypes.

The same `sweep_CDS_bp` and `control_CDS_bp` denominators are used for every genotype.

## Coordinate Conventions

Input coordinates:

```text
GFF CDS coordinates: 1-based inclusive
pi-ratio windows:    1-based inclusive
SNP POS values:      1-based
```

BED coordinates:

```text
0-based half-open
```

SNP variants are converted to BED as:

```text
BED start = POS - 1
BED end   = POS
```

## How To Run

From the cluster login node:

```bash
cd ~/work/sweep_gl
bash submit_pipeline.sh
```

The submit script creates required output directories and submits four SLURM jobs with `afterok` dependencies:

```text
01_make_sweep_regions
02_make_cds_regions
03_filter_assign_load
04_calc_genotype_burden
```

The Python scripts use only the Python standard library. No pandas or conda environment is required.

The pipeline loads these cluster modules from `scripts/load_cluster_modules.sh`:

```text
devel/python/Python-3.11.1
bioinfo/bedtools/2.31.1
```

`bedtools` is required for interval operations.

The pipeline does not assume conda is available.

## Final Output

```text
results/final_genotype_sweep_control_burden.tsv
```

This file has exactly three columns:

```text
genotype
sweep_deleterious_dosage_per_100kb_CDS
control_deleterious_dosage_per_100kb_CDS
```

## QC Outputs

```text
results/qc/sweep_region_qc.tsv
results/qc/cds_lengths.tsv
results/qc/load_variant_counts.tsv
results/qc/genotype_burden_extended.tsv
```

Additional region and variant outputs:

```text
results/regions/sweep_windows_selected.bed
results/regions/sweep_regions_merged.bed
results/regions/cds_raw.bed
results/regions/cds_reduced.bed
results/regions/sweep_cds.bed
results/regions/control_cds.bed
results/variants/load_variants_sweep.tsv
results/variants/load_variants_control.tsv
```

## Pipeline Files

```text
scripts/01_make_sweep_regions.py
scripts/02_make_cds_regions.sh
scripts/03_filter_and_assign_load_variants.py
scripts/04_calculate_genotype_burden.py
slurm/01_make_sweep_regions.sbatch
slurm/02_make_cds_regions.sbatch
slurm/03_filter_and_assign_load_variants.sbatch
slurm/04_calculate_genotype_burden.sbatch
submit_pipeline.sh
```
