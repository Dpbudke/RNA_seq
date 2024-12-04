#!/bin/bash

#SBATCH --job-name=compress_emase
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 08:00:00
#SBATCH --mem=32G
#SBATCH --array=1-246
#SBATCH --output=slurmout/compress/compress_%A_%a.out
#SBATCH --error=slurmout/compress/compress_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dpbudke@ucdavis.edu

module load python_3/3.11.1

# Get the sample ID from the samples.txt file based on array task ID
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

# Define input and output paths
INPUT_DIR="common_alignments_output"
OUTPUT_DIR="compressed_emase_output"
mkdir -p ${OUTPUT_DIR}

# Set file paths
EMASE_FILE="${INPUT_DIR}/${SAMPLE}_common.h5"
COMPRESSED_EMASE_FILE="${OUTPUT_DIR}/${SAMPLE}_compressed.h5"

# Print sample being processed
echo "Processing sample: ${SAMPLE}"
echo "Input file: ${EMASE_FILE}"
echo "Output file: ${COMPRESSED_EMASE_FILE}"

# Run GBRS compress
gbrs compress \
    -i ${EMASE_FILE} \
    -o ${COMPRESSED_EMASE_FILE}

# Check exit status
if [ $? -eq 0 ]; then
    echo "Successfully compressed ${SAMPLE}"
else
    echo "Error compressing ${SAMPLE}"
    exit 1
fi
