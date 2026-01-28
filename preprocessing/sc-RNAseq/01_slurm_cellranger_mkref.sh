#!/bin/bash
#SBATCH --job-name=cellranger_mkref
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --partition=medium
#SBATCH --output=cellranger_mkref_%j.out
#SBATCH --error=cellranger_mkref_%j.err

echo "Job name: $SLURM_JOB_NAME"
echo "Job ID: $SLURM_JOB_ID"
echo "Started at: $(date)"
echo "Running on: $(hostname)"
echo "CPUs: $SLURM_NTASKS"

############################
# USER-DEFINED PARAMETERS #
############################

# Base project directory (CHANGE THIS)
PROJECT_DIR=/scratch/path/to/project

# Reference input files
REF_DIR=${PROJECT_DIR}/reference

# Output directory for Cell Ranger references
CR_INDEX_DIR=${PROJECT_DIR}/cellranger_index

# Memory for cellranger mkref
MEM_GB=100

############################
# REFERENCES               #
############################

# Human
HUMAN_FASTA_GZ=${REF_DIR}/human_genome.fa.gz
HUMAN_GTF_GZ=${REF_DIR}/human_annotation.gtf.gz

# Marmoset (standard / non-liftoff)
MARMOSET_FASTA_GZ=${REF_DIR}/marmoset_genome.fa.gz
MARMOSET_GTF_GZ=${REF_DIR}/marmoset_annotation.gtf.gz

############################
# LOAD MODULES             #
############################

module load cellranger

############################
# SETUP DIRECTORIES        #
############################

mkdir -p \
  ${CR_INDEX_DIR}/human \
  ${CR_INDEX_DIR}/marmoset \
  ${REF_DIR}/decompressed

############################
# DECOMPRESS REFERENCES    #
############################

gunzip -c ${HUMAN_FASTA_GZ} > ${REF_DIR}/decompressed/human.fa
gunzip -c ${HUMAN_GTF_GZ} > ${REF_DIR}/decompressed/human.gtf

gunzip -c ${MARMOSET_FASTA_GZ} > ${REF_DIR}/decompressed/marmoset.fa
gunzip -c ${MARMOSET_GTF_GZ} > ${REF_DIR}/decompressed/marmoset.gtf

############################
# CELLRANGER MKREF         #
############################

echo "Running cellranger mkref: Human"

cd ${CR_INDEX_DIR}/human

cellranger mkref \
  --genome=human \
  --fasta=${REF_DIR}/decompressed/human.fa \
  --genes=${REF_DIR}/decompressed/human.gtf \
  --memgb=${MEM_GB}

echo "Running cellranger mkref: Marmoset"

cd ${CR_INDEX_DIR}/marmoset

cellranger mkref \
  --genome=marmoset \
  --fasta=${REF_DIR}/decompressed/marmoset.fa \
  --genes=${REF_DIR}/decompressed/marmoset.gtf \
  --memgb=${MEM_GB}

############################
# CLEANUP (OPTIONAL)       #
############################

rm -f \
  ${REF_DIR}/decompressed/human.fa \
  ${REF_DIR}/decompressed/human.gtf \
  ${REF_DIR}/decompressed/marmoset.fa \
  ${REF_DIR}/decompressed/marmoset.gtf

echo "Cell Ranger reference generation completed at: $(date)"
