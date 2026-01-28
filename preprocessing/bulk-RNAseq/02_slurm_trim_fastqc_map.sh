#!/bin/bash
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=star_trim_map
#SBATCH --mem=150G
#SBATCH --time=48:00:00
#SBATCH --output=star_trim_map_%j.out
#SBATCH --error=star_trim_map_%j.err

############################
# USER-DEFINED PARAMETERS #
############################

# Base project directory (change once per dataset)
PROJECT_DIR=/scratch/path/to/project

# Input FASTQ parent directory (contains subfolders per sample)
READS_DIR=${PROJECT_DIR}/reads

# Output directories
TRIMMED_DIR=${PROJECT_DIR}/trimmed_reads
FASTQC_DIR=${PROJECT_DIR}/fastqc
STAR_OUT_DIR=${PROJECT_DIR}/star_output
STAR_TMP_DIR=${PROJECT_DIR}/star_tmp

# STAR genome index
GENOME_INDEX=/scratch/path/to/STAR_genome_index

####################
# LOAD MODULES     #
####################

module load cutadapt/4.4
module load fastqc/0.12.1
module load star/2.7.3a

####################
# FUNCTIONS        #
####################

extract_sample_id() {
    basename "$1" | sed 's/_R1_001.fastq.gz//'
}

process_folder() {
    folder="$1"
    sample_group=$(basename "$folder")

    echo "Processing folder: ${sample_group}"

    mkdir -p \
        "${TRIMMED_DIR}/${sample_group}" \
        "${STAR_OUT_DIR}/${sample_group}" \
        "${FASTQC_DIR}/${sample_group}"

    cd "$folder" || exit 1

    for r1 in *_R1_001.fastq.gz; do
        [ -e "$r1" ] || continue
        sample=$(extract_sample_id "$r1")

        r2="${sample}_R2_001.fastq.gz"

        echo "Processing sample: ${sample}"

        # Trimming
        cutadapt \
            -o "${TRIMMED_DIR}/${sample_group}/${sample}_R1_trimmed.fastq.gz" \
            -p "${TRIMMED_DIR}/${sample_group}/${sample}_R2_trimmed.fastq.gz" \
            "$r1" "$r2"

        # FastQC
        fastqc \
            --outdir "${FASTQC_DIR}/${sample_group}" \
            "${TRIMMED_DIR}/${sample_group}/${sample}_R1_trimmed.fastq.gz" \
            "${TRIMMED_DIR}/${sample_group}/${sample}_R2_trimmed.fastq.gz"

        # Decompress for STAR
        gunzip -c "${TRIMMED_DIR}/${sample_group}/${sample}_R1_trimmed.fastq.gz" \
            > "${TRIMMED_DIR}/${sample_group}/${sample}_R1_trimmed.fastq"

        gunzip -c "${TRIMMED_DIR}/${sample_group}/${sample}_R2_trimmed.fastq.gz" \
            > "${TRIMMED_DIR}/${sample_group}/${sample}_R2_trimmed.fastq"

        rm -rf "${STAR_TMP_DIR}"

        # STAR alignment
        STAR \
            --genomeDir "${GENOME_INDEX}" \
            --readFilesIn \
                "${TRIMMED_DIR}/${sample_group}/${sample}_R1_trimmed.fastq" \
                "${TRIMMED_DIR}/${sample_group}/${sample}_R2_trimmed.fastq" \
            --runThreadN ${SLURM_CPUS_PER_TASK} \
            --quantMode GeneCounts \
            --outSAMtype None \
            --outTmpDir "${STAR_TMP_DIR}" \
            --outFileNamePrefix "${STAR_OUT_DIR}/${sample_group}/${sample}_" \
            --outFilterMatchNmin 20 \
            --outFilterMatchNminOverLread 0.3 \
            --outFilterMismatchNmax 10 \
            --outFilterMismatchNoverLmax 0.2 \
            --alignEndsType Local \
            --outFilterScoreMinOverLread 0.3

        # Cleanup
        rm \
            "${TRIMMED_DIR}/${sample_group}/${sample}_R1_trimmed.fastq" \
            "${TRIMMED_DIR}/${sample_group}/${sample}_R2_trimmed.fastq"
    done
}

####################
# MAIN LOOP        #
####################

for folder in "${READS_DIR}"/*/; do
    [ -d "$folder" ] && process_folder "$folder"
done
