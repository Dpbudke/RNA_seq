#!/bin/bash

#===============================================================================
# SLURM job configuration
#===============================================================================
#SBATCH --job-name=bam2emase
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 08:00:00
#SBATCH --mem=32G
#SBATCH --array=1-246
#SBATCH --output=slurmout/bam2emase/bam2emase_%A_%a.out
#SBATCH --error=slurmout/bam2emase/bam2emase_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dpbudke@ucdavis.edu

#===============================================================================
# Environment setup
#===============================================================================
# Load necessary modules
module load samtools

# Define directories and settings
BAM_DIR="/90daydata/do2_projects/CCPups_redo/EMASE/sorted_bams"
GBRS_DATA="/90daydata/do2_projects/CCPups_redo/EMASE/GBRS_DATA"
HAPLOTYPE_CODES="A,B,C,D,E,F,G,H"
EMASE_OUTPUT_DIR="emase_output"

# Create output directory if it doesn't exist
mkdir -p "$EMASE_OUTPUT_DIR"

#===============================================================================
# Get sample information
#===============================================================================
# Get sample ID from samples.txt based on array task ID
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

#===============================================================================
# Main processing
#===============================================================================
# Process both R1 and R2 BAM files for the sample
for READ in R1 R2; do
    # Updated BAM filename format to match new structure (.sorted.bam)
    BAM_FILE="${BAM_DIR}/${SAMPLE}_${READ}.sorted.bam"
    EMASE_FILE="${EMASE_OUTPUT_DIR}/${SAMPLE}_${READ}.emase"

    if [ -f "$BAM_FILE" ]; then
        # Run gbrs bam2emase
        echo "Processing sample: $SAMPLE ($READ)"
        gbrs bam2emase \
            -i "$BAM_FILE" \
            -m "$GBRS_DATA"/emase.fullTranscripts.info \
            -h "$HAPLOTYPE_CODES" \
            -o "$EMASE_FILE"
    else
        echo "Error: BAM file not found for sample $SAMPLE ($READ): $BAM_FILE"
        exit 1
    fi
done

echo "Job completed for sample: $SAMPLE"
