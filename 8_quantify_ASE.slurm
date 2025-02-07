#!/bin/bash

#SBATCH --job-name=quant_ASE
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 08:00:00
#SBATCH --mem=32G
#SBATCH --array=1-246
#SBATCH --output=slurmout/quant_ASE/ASE_%A_%a.out
#SBATCH --error=slurmout/quant_ASE/ASE_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dpbudke@ucdavis.edu

# Set working directory
cd /90daydata/do2_projects/CCPups_redo/EMASE/

# Get the sample ID from the samples.txt file based on SLURM array task ID
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

# Create a temporary directory in work area instead of /tmp
TEMP_DIR="/90daydata/do2_projects/CCPups_redo/EMASE/temp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p ${TEMP_DIR}
cd ${TEMP_DIR}

# Define paths (using absolute paths since we're in temp directory)
WORK_DIR="/90daydata/do2_projects/CCPups_redo/EMASE"
GBRS_DATA="${WORK_DIR}/GBRS_DATA"
MERGED_COMPRESSED_EMASE_FILE="${WORK_DIR}/compressed_emase_output/${SAMPLE}_compressed.h5"
GENOTYPES_FILE="${WORK_DIR}/reconstruct_output/${SAMPLE}/${SAMPLE}.genotypes.tsv"

# Create final output directory structure
OUTPUT_DIR="${WORK_DIR}/quantify_ASE/${SAMPLE}"
mkdir -p ${OUTPUT_DIR}

# Run GBRS quantify
gbrs quantify \
    -i ${MERGED_COMPRESSED_EMASE_FILE} \
    -G ${GENOTYPES_FILE} \
    -g ${GBRS_DATA}/emase.gene2transcripts.tsv \
    -L ${GBRS_DATA}/emase.pooled.fullTranscripts.info \
    -M 4 \
    --report-alignment-counts

# Move all output files with sample prefix
for file in gbrs.quantified.diploid.genes.alignment_counts \
            gbrs.quantified.diploid.genes.expected_read_counts \
            gbrs.quantified.diploid.genes.tpm \
            gbrs.quantified.diploid.isoforms.alignment_counts \
            gbrs.quantified.diploid.isoforms.expected_read_counts \
            gbrs.quantified.diploid.isoforms.tpm; do
    if [ -f "$file" ]; then
        mv "$file" "${OUTPUT_DIR}/${SAMPLE}_${file}"
    fi
done

# Clean up temp directory
cd ${WORK_DIR}
rm -rf ${TEMP_DIR}
