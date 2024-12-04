#!/bin/bash

#SBATCH --job-name=quantify
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 08:00:00
#SBATCH --mem=32G
#SBATCH --array=1-246
#SBATCH --output=slurmout/quantify/quant_%A_%a.out
#SBATCH --error=slurmout/quantify/quant_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dpbudke@ucdavis.edu

# Directory settings
GBRS_DATA="GBRS_DATA"  
WORK_DIR="/90daydata/do2_projects/CCPups_redo/EMASE" 
INPUT_DIR="compressed_emase_output"
OUTPUT_DIR="${WORK_DIR}/quantify_output"

# Get the sample name from samples.txt based on array task ID
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

# Input file
COMPRESSED_EMASE_FILE="${WORK_DIR}/${INPUT_DIR}/${SAMPLE}_compressed.h5"

# Create output directory if it doesn't exist
mkdir -p ${WORK_DIR}/quantify_output

# Create sample-specific subdirectory
SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/${SAMPLE}"
mkdir -p ${SAMPLE_OUTPUT_DIR}

# Run GBRS quantify
gbrs quantify \
    -i ${COMPRESSED_EMASE_FILE} \
    -g ${GBRS_DATA}/emase.gene2transcripts.tsv \
    -L ${GBRS_DATA}/emase.pooled.fullTranscripts.info \
    -M 4 \
    --report-alignment-counts \
    -o ${SAMPLE_OUTPUT_DIR}/${SAMPLE}

# Print completion message
echo "Completed processing ${SAMPLE}"
