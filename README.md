# Transcriptomic analysis of cerebral organoids and tissue (scRNA-seq & bulk RNA-seq) done in Tynianskaia et al. (2026)

This repository contains R scripts used to analyze single-cell and bulk RNA sequencing data from embryoid bodies, cerebral organoids, and tissue samples.

## Overview

This project investigates cell composition, developmental trajectories, and functional networks of marmoset and human cerebral organoid and marmoset tissue using scRNA-seq and bulk RNA-seq.

## Repository structure
RNAseq_analysis/
├── preprocessing/
│   ├── sc-RNAseq/
│   │   ├── 01_slurm_cellranger_mkref.sh
│   │   ├── 02_slurm_cellranger_count.sh
│   │   └── 03_slurm_cellranger_aggregate.sh
│   │
│   ├── bulk-RNAseq/
│       ├── 01_slurm_genome_index.sh
│       └── 02_slurm_trim_fastqc_map.sh
│   
├── sc-RNAseq/
│   ├── 01_seurat_pipeline_integrated_d30_orgs.R
│   ├── 02_seurat_pipeline_marmoset_d50_orgs.R
│   ├── 03_seurat_pipeline_marmoset_GD90.R
│   ├── 04_cell_composition_analysis.R
│   ├── 05_functional_profiling_bRGs_d30.R
│   └── 06_sc_data_integration_d30_d50_GD90.R
│
├── bulk-RNAseq/
│   ├── 01_EBs_d3_d6_human_marmoset.R
│   ├── 02_marmoset_org_FACS_aRGs.R
│   ├── 03_orgs_brains_aRGs_Neu.R
│   └── metadata/
│       └── pheno_data_EBs.csv
│       └── pheno_data_aRG_orgs.csv
│       └── pheno_data_marmoset_orgs&brains.csv
│
├── renv.lock
└── README.md

## Requirements

- R version: 4.5.1  
- Package versions are tracked using `renv`.

To restore the environment:
```r
install.packages("renv")
renv::restore()
