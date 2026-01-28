#!/bin/bash
#SBATCH --job-name=STAR_GenomeIndex
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=320G
#SBATCH --time=24:00:00
#SBATCH --partition=medium
#SBATCH --output=star_genome_index_%j.out
#SBATCH --error=star_genome_index_%j.err

echo "Job: $SLURM_JOB_NAME"
echo "Job ID: $SLURM_JOB_ID"
echo "Started at: $(date)"
echo "Running on: $(hostname)"
echo "CPUs: $SLURM_NTASKS"

############################
# USER-DEFINED PARAMETERS #
############################

# Base project directory
PROJECT_DIR=/scratch/path/to/project

# Directory containing genome fasta and GTF files
REF_DIR=${PROJECT_DIR}/reference

# Output directory for STAR indexes
INDEX_DIR=${PROJECT_DIR}/star_index

# STAR parameters
SJDB_OVERHANG=149
THREADS=${SLURM_NTASKS}
RAM_LIMIT=31000000000

############################
# REFERENCES               #
############################

# Human
HUMAN_FASTA=${REF_DIR}/human_genome.fa.gz
HUMAN_GTF=${REF_DIR}/human_annotation.gtf.gz

# Marmoset (standard)
MARMOSET_FASTA=${REF_DIR}/marmoset_genome.fa.gz
MARMOSET_GTF=${REF_DIR}/marmoset_annotation.gtf.gz

# Marmoset (liftoff / UCSC-style)
MARMOSET_UCSC_FASTA=${REF_DIR}/marmoset_ucsc_genome.fa
MARMOSET_LIFTOFF_GTF=${REF_DIR}/marmoset_liftoff_annotation.gtf.gz

############################
# LOAD MODULES             #
############################

module load star/2.7.3a

############################
# PREPARE FILES            #
############################

mkdir -p \
  ${INDEX_DIR}/human \
  ${INDEX_DIR}/marmoset \
  ${INDEX_DIR}/marmoset_liftoff \
  ${REF_DIR}/decompressed

# Decompress references (STAR prefers uncompressed)
gunzip -c ${HUMAN_FASTA} > ${REF_DIR}/decompressed/human.fa
gunzip -c ${HUMAN_GTF} > ${REF_DIR}/decompressed/human.gtf

gunzip -c ${MARMOSET_FASTA} > ${REF_DIR}/decompressed/marmoset.fa
gunzip -c ${MARMOSET_GTF} > ${REF_DIR}/decompressed/marmoset.gtf

gunzip -c ${MARMOSET_LIFTOFF_GTF} > ${REF_DIR}/decompressed/marmoset_liftoff.gtf

############################
# STAR GENOME INDEXING     #
############################

echo "Generating STAR index: Human"

STAR \
  --runMode genomeGenerate \
  --genomeDir ${INDEX_DIR}/human \
  --genomeFastaFiles ${REF_DIR}/decompressed/human.fa \
  --sjdbGTFfile ${REF_DIR}/decompressed/human.gtf \
  --sjdbOverhang ${SJDB_OVERHANG} \
  --runThreadN ${THREADS} \
  --limitGenomeGenerateRAM ${RAM_LIMIT}

echo "Generating STAR index: Marmoset (standard)"

STAR \
  --runMode genomeGenerate \
  --genomeDir ${INDEX_DIR}/marmoset \
  --genomeFastaFiles ${REF_DIR}/decompressed/marmoset.fa \
  --sjdbGTFfile ${REF_DIR}/decompressed/marmoset.gtf \
  --sjdbOverhang ${SJDB_OVERHANG} \
  --runThreadN ${THREADS} \
  --limitGenomeGenerateRAM ${RAM_LIMIT}

echo "Generating STAR index: Marmoset (liftoff / UCSC)"

STAR \
  --runMode genomeGenerate \
  --genomeDir ${INDEX_DIR}/marmoset_liftoff \
  --genomeFastaFiles ${MARMOSET_UCSC_FASTA} \
  --sjdbGTFfile ${REF_DIR}/decompressed/marmoset_liftoff.gtf \
  --sjdbOverhang ${SJDB_OVERHANG} \
  --runThreadN ${THREADS} \
  --limitGenomeGenerateRAM ${RAM_LIMIT}

echo "Genome indexing completed at: $(date)"
