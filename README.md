# Transcriptomic analysis of cerebral organoids and tissue (scRNA-seq & bulk RNA-seq) done in Tynianskaia et al. (2026)

This repository contains R scripts used to analyze single-cell and bulk RNA sequencing data from embryoid bodies, cerebral organoids, and tissue samples.

## Overview

This project investigates cell composition, developmental trajectories, and functional networks of marmoset and human cerebral organoid and marmoset tissue using scRNA-seq and bulk RNA-seq.

## Repository structure
```text
├── preprocessing/
│   ├── sc-RNAseq/
│   │   ├── 01_slurm_cellranger_mkref.sh
│   │   ├── 02_slurm_cellranger_count.sh
│   │   └── 03_slurm_cellranger_aggregate.sh
│   │
│   ├── bulk-RNAseq/
│   │   ├── 01_slurm_genome_index.sh
│   │   └── 02_slurm_trim_fastqc_map.sh
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
│   └── 03_orgs_GD90_aRGs_Neu.R
│
├── metadata/
│   ├── pheno_data_EBs.csv
│   ├── pheno_data_aRG_orgs.csv
│   └── pheno_data_marmoset_orgs&GD90.csv
│
├── renv.lock
└── README.md
```
## System Requirements
## Software

R 4.5.1

renv (for package version control)

Seurat

DESeq2

ggplot2

dplyr

Package versions are recorded in renv.lock.

Operating System: Linux (HPC cluster environment and local testing)

## Hardware
Downstream R analyses can run on a standard workstation (≥16 GB RAM recommended).

Preprocessing scripts require an HPC environment with a SLURM scheduler.

## Installation Guide

Clone the repository and restore the R environment:

install.packages("renv")

renv::restore()

Estimated installation time: 5–10 minutes on a standard desktop computer with internet access.

## Demo

A small demo dataset is provided in /demo_data.

To run the demo:

source("demo/run_sc_demo.R")

What the demo does: 
Loads a small Seurat object

Performs normalization and clustering

Generates a UMAP plot

## Expected runtime
~5 minutes on a standard laptop (16 GB RAM).

## Instructions for Use
- scRNA-seq analysis

Run preprocessing scripts in preprocessing/sc-RNAseq/ (HPC required)

Run Seurat scripts in sc-RNAseq/ in numerical order

- Bulk RNA-seq analysis

Run preprocessing scripts in preprocessing/bulk-RNAseq/

Run R scripts in bulk-RNAseq/ in numerical order

## Reproducing Manuscript Results
The following scripts correspond to key analyses presented in the manuscript:

Integrated scRNA-seq analysis (d30 organoids)	01_seurat_pipeline_integrated_d30_orgs.R

Cell composition analysis	04_cell_composition_analysis.R

Functional profiling of bRGs	05_functional_profiling_bRGs_d30.R

Multi-stage scRNA-seq integration	06_sc_data_integration_d30_d50_GD90.R

Bulk RNA-seq analyses	Scripts in bulk-RNAseq/

Running scripts in the indicated order on the full dataset (ENA accession PRJEB107058) reproduces the computational analyses presented in the study.
