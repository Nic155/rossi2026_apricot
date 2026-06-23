# rossi2026_apricot

Analysis code for **Rossi et al. (2026), "Reduced deleterious load and tandem
gene-family signals in domesticated European apricot."**

This repository covers the downstream analyses: runs of homozygosity (ROH),
ROH-vs-load burden, selective-sweep burden, and the R analyses producing
Figures 1-2 and the GO/KEGG enrichment. It tests whether European cultivated
apricot (C1, n = 50) carries higher or lower deleterious load than its
Central Asian wild progenitor (W2, n = 43).

## Upstream dependency

The per-variant deleterious classification is produced by a separate pipeline:
**WP5_genetic_load** (https://github.com/Nic155/WP5_genetic_load). Its output,
`result3-provean_W2_C1_subset.tsv`, is the main input consumed by the stages
here. That file is not stored in git (1.3 GB); see Data below.

## Pipeline stages

Stages run in order. ROH, ROH-load and sweep-burden run on a SLURM HPC; the
R stage runs locally.

| Stage | Folder | What it does |
|-------|--------|--------------|
| ROH | `roh/` | VCF to PLINK, filter, thin, call ROH at several length thresholds, compute FROH |
| ROH-load | `roh_load/` | Per-individual deleterious dosage inside vs outside ROH (per 100 kb CDS) |
| Sweep burden | `sweep_burden/` | Deleterious dosage in selective-sweep regions vs control CDS |
| R analysis | `R_analysis/` | Figures 1-2: load, inbreeding, PROVEAN density, sweep burden, GO/KEGG enrichment, chr4 gene cloud |

`roh_load/` and `sweep_burden/` each contain their own detailed README.

## Data

Code lives here; bulk data is deposited externally.

**Derived data (this study)** - Recherche Data Gouv,
DOI 10.57745/D1LNG3 (https://doi.org/10.57745/D1LNG3):
the deleterious-variant table (`result3-...subset.tsv`), the R workspaces
(`.RData`) loaded by the R scripts, the sample list, and `roh_800kb.hom`.

**Raw sequencing reads (Groppi et al. 2021)** - European Nucleotide Archive:
Illumina DNASeq PRJEB42181 and PRJEB40984 (plus PRJEB42606, PRJEB40668,
PRJEB42479 for the de novo and RNASeq data).

**Reference genome & annotation (Marouch #14)** - Genome Database for Rosaceae:
https://www.rosaceae.org/Analysis/9642068 (source of `genes.gff`).

Small inputs (`data/metadata/list.txt`, `data/sweep_input/pi_ratio_W2_C1.csv`)
are committed here. Large inputs and PLINK/ROH binaries are gitignored.

## Environments

- HPC: PLINK 1.90b7, bcftools 1.21, bedtools 2.31.1, Python 3.11.1 (stdlib only).
- R: clusterProfiler 4.20.0 and the packages loaded at the top of each script.

## Citation

Cite Rossi et al. (2026) and the data DOI above.
