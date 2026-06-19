# ROH Deleterious Burden Pipeline

## Purpose

This SLURM-based command-line pipeline calculates deleterious genetic-load burden for every genotype by comparing CDS regions inside each individual's runs of homozygosity (ROH) against CDS regions outside ROH. The final burden is reported in the soybean-paper style:

```text
deleterious genotype dosage per 100 kb CDS
```

Because ROH regions are per-individual, the inside-ROH vs outside-ROH classification is done separately for each sample. The final output has one row per genotype and six columns:

```text
genotype
population
roh_deleterious_dosage_per_100kb_CDS
nonroh_deleterious_dosage_per_100kb_CDS
roh_cds_bp
nonroh_cds_bp
```

## Input Files

Expected paths:

```text
~/work/ROH/roh/roh_800kb.hom
~/work/ROH/result3-provean_forGL+score+class.C1_W2.subset.tsv
~/work/genetic_load/parmeniaca/5-ReferenceData/Prunus_armeniaca_Marouch_n14/genes.gff
~/work/ROH/list.txt
```

## ROH Region Definition

ROH intervals are read from the PLINK `.hom` file. This file is whitespace-delimited with a header. The relevant columns are:

```text
IID  (column 2) -- sample identifier
CHR  (column 4) -- chromosome as an integer
POS1 (column 7) -- ROH start, 1-based inclusive
POS2 (column 8) -- ROH end, 1-based inclusive
```

Chromosome values are converted from integers to BED chromosome format:

```text
1 -> chr1
2 -> chr2
...
```

Coordinates are converted from 1-based inclusive to 0-based half-open BED:

```text
BED start = POS1 - 1
BED end   = POS2
```

One sorted BED file is written per individual to:

```text
results/regions/roh_per_sample/{sample}.bed
```

## Per-Individual CDS Intersection

For each individual, the pipeline computes:

```text
ROH CDS     = CDS regions intersecting that individual's ROH intervals
non-ROH CDS = CDS regions minus that individual's ROH intervals
```

These are written to:

```text
results/regions/roh_per_sample_cds/{sample}.bed
results/regions/nonroh_per_sample_cds/{sample}.bed
```

Per-individual ROH CDS bp and non-ROH CDS bp are written to:

```text
results/qc/per_individual_cds_lengths.tsv
```

Individuals absent from the `.hom` file (no ROH intervals) have `roh_cds_bp = 0`
and their `nonroh_cds_bp` equals the total reduced CDS bp. Their
`roh_deleterious_dosage_per_100kb_CDS` is reported as `NA`.

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

## Formula

For every genotype:

```text
roh_deleterious_dosage_per_100kb_CDS =
    sum genotype dosage for PROVEAN/HIGH variants in that individual's ROH CDS
    / that individual's ROH CDS bp * 100000

nonroh_deleterious_dosage_per_100kb_CDS =
    sum genotype dosage for PROVEAN/HIGH variants in that individual's non-ROH CDS
    / that individual's non-ROH CDS bp * 100000
```

NA genotype calls are ignored in dosage sums. NA is not treated as zero. The denominator is per-individual CDS length, not the number of called genotypes.

## Coordinate Conventions

Input coordinates:

```text
GFF CDS coordinates: 1-based inclusive
PLINK .hom POS1/POS2: 1-based inclusive
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
cd ~/work/ROH/roh_gl
bash submit_pipeline.sh
```

The submit script creates required output directories and submits four SLURM jobs with `afterok` dependencies:

```text
01_make_cds_regions
02_make_roh_regions
03_filter_load_variants
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
results/final_genotype_roh_nonroh_burden.tsv
```

This file has exactly six columns:

```text
genotype
population
roh_deleterious_dosage_per_100kb_CDS
nonroh_deleterious_dosage_per_100kb_CDS
roh_cds_bp
nonroh_cds_bp
```

## QC Outputs

```text
results/qc/cds_lengths.tsv
results/qc/per_individual_cds_lengths.tsv
results/qc/load_variant_counts.tsv
results/qc/genotype_burden_extended.tsv
```

Additional region and variant outputs:

```text
results/regions/cds_raw.bed
results/regions/cds_reduced.bed
results/regions/roh_per_sample/{sample}.bed
results/regions/roh_per_sample_cds/{sample}.bed
results/regions/nonroh_per_sample_cds/{sample}.bed
results/variants/load_variants_cds.tsv
results/variants/load_variants_cds.bed
```

## Pipeline Files

```text
scripts/01_make_cds_regions.sh
scripts/02_make_roh_regions.py
scripts/03_filter_load_variants.py
scripts/04_calculate_genotype_burden.py
slurm/01_make_cds_regions.sbatch
slurm/02_make_roh_regions.sbatch
slurm/03_filter_load_variants.sbatch
slurm/04_calculate_genotype_burden.sbatch
submit_pipeline.sh
```
