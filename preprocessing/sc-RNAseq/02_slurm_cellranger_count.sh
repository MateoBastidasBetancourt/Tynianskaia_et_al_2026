#!/bin/bash
#SBATCH --job-name=cellranger_count_intronic
#SBATCH --mail-type=BEGIN,END
#SBATCH --ntasks=48
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --partition=medium
#SBATCH --output=cellranger_count_intronic_%j.out
#SBATCH --error=cellranger_count_intronic_%j.err

echo "Job name: $SLURM_JOB_NAME"
echo "Job ID: $SLURM_JOB_ID"
echo "Host: $(hostname)"
echo "Start time: $(date)"
echo "Working directory: $(pwd)"
echo "CPUs: $SLURM_NTASKS"

############################
# LOAD MODULES
############################

module load cellranger

############################
# BASE DIRECTORIES
############################

PROJECT_DIR=/scratch/path/to/project

REF_DIR=${PROJECT_DIR}/references
FASTQ_HUMAN=${PROJECT_DIR}/fastqs/human
FASTQ_MARMOSET=${PROJECT_DIR}/fastqs/marmoset

OUT_HUMAN=${PROJECT_DIR}/cellranger_output/human
OUT_MARMOSET=${PROJECT_DIR}/cellranger_output/marmoset

REF_HUMAN=${REF_DIR}/index_human
REF_MARMOSET=${REF_DIR}/index_marmoset

############################
# HUMAN DATA (INTRONIC)
############################

mkdir -p ${OUT_HUMAN}/out_intr
cd ${OUT_HUMAN}/out_intr

cellranger count \
  --id=sc \
  --transcriptome=${REF_HUMAN} \
  --fastqs=${FASTQ_HUMAN} \
  --sample=sc \
  --expect-cells=10000 \
  --localmem=120 \
  --no-bam \
  --include-introns

cellranger count \
  --id=IDP \
  --transcriptome=${REF_HUMAN} \
  --fastqs=${FASTQ_HUMAN} \
  --sample=IDP \
  --expect-cells=10000 \
  --localmem=120 \
  --no-bam \
  --include-introns

############################
# iLonza FASTQ MERGING (REQUIRED)
############################
# NOTE: This step is essential because reads are split across lanes

cd ${FASTQ_HUMAN}

cat iLonza36_S1_L144168_R1_001.fastq.gz \
    iLonza98_S1_L144168_R1_001.fastq.gz \
    > iLonza_S1_L144168_R1_001.fastq.gz

cat iLonza36_S1_L144168_R2_001.fastq.gz \
    iLonza98_S1_L144168_R2_001.fastq.gz \
    > iLonza_S1_L144168_R2_001.fastq.gz

cd ${OUT_HUMAN}/out_intr

cellranger count \
  --id=iLonza \
  --transcriptome=${REF_HUMAN} \
  --fastqs=${FASTQ_HUMAN} \
  --sample=iLonza \
  --expect-cells=10000 \
  --localmem=120 \
  --no-bam \
  --include-introns

############################
# MARMOSET DATA (INTRONIC)
############################

mkdir -p ${OUT_MARMOSET}/out_intr
cd ${OUT_MARMOSET}/out_intr

cellranger count \
  --id=cj4 \
  --transcriptome=${REF_MARMOSET} \
  --fastqs=${FASTQ_MARMOSET} \
  --sample=cj4 \
  --expect-cells=10000 \
  --localmem=200 \
  --no-bam \
  --include-introns

cellranger count \
  --id=cj5 \
  --transcriptome=${REF_MARMOSET} \
  --fastqs=${FASTQ_MARMOSET} \
  --sample=cj5 \
  --expect-cells=10000 \
  --localmem=200 \
  --no-bam \
  --include-introns

cellranger count \
  --id=cj6 \
  --transcriptome=${REF_MARMOSET} \
  --fastqs=${FASTQ_MARMOSET} \
  --sample=cj6 \
  --expect-cells=10000 \
  --localmem=200 \
  --no-bam \
  --include-introns

echo "Job finished at: $(date)"