#!/bin/bash
#SBATCH --job-name=cellranger_aggr
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --partition=medium
#SBATCH --output=cellranger_aggr_%j.out
#SBATCH --error=cellranger_aggr_%j.err

echo "Job: $SLURM_JOB_NAME"
echo "Job ID: $SLURM_JOB_ID"
echo "Host: $(hostname)"
echo "Start time: $(date)"
echo "Working directory: $(pwd)"

############################
# USER PARAMETERS          #
############################

# Base project directory (CHANGE THIS)
PROJECT_DIR=/scratch/path/to/project

# Aggregation CSV files
AGGR_HUMAN=${PROJECT_DIR}/aggregation/aggregation_sc_human.csv
AGGR_MARMOSET=${PROJECT_DIR}/aggregation/aggregation_sc_marmoset.csv

# Output directories
OUT_HUMAN=${PROJECT_DIR}/cellranger_output/human/aggr
OUT_MARMOSET=${PROJECT_DIR}/cellranger_output/marmoset/aggr

############################
# LOAD MODULES             #
############################

module load cellranger

############################
# HUMAN AGGREGATION        #
############################

mkdir -p ${OUT_HUMAN}
cd ${OUT_HUMAN}

cellranger aggr \
  --id=human_aggr \
  --csv=${AGGR_HUMAN}

############################
# MARMOSET AGGREGATION     #
############################

mkdir -p ${OUT_MARMOSET}
cd ${OUT_MARMOSET}

cellranger aggr \
  --id=marmoset_aggr \
  --csv=${AGGR_MARMOSET}

echo "Aggregation finished at: $(date)"
