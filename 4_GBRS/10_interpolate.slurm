#!/bin/bash

#SBATCH --job-name=interpolate
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 08:00:00
#SBATCH --mem=32G
#SBATCH --array=1-246
#SBATCH --output=slurmout/interpolate/_%A_%a.out
#SBATCH --error=slurmout/interpolate/_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dpbudke@ucdavis.edu

# Set working directory
cd /90daydata/do2_projects/CCPups_redo/EMASE/

# Get the sample ID from the samples.txt file based on SLURM array task ID
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

# Create a temporary directory in work area
TEMP_DIR="/90daydata/do2_projects/CCPups_redo/EMASE/temp_interpolate_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p ${TEMP_DIR}
cd ${TEMP_DIR}

# Define paths (using absolute paths since we're in temp directory)
WORK_DIR="/90daydata/do2_projects/CCPups_redo/EMASE"
export GBRS_DATA="${WORK_DIR}/GBRS_DATA"

# Create output directory if it doesn't exist
mkdir -p ${WORK_DIR}/genome_output

# Interpolate genotypes
gbrs interpolate \
    -i ${WORK_DIR}/reconstruct_output/${SAMPLE}/${SAMPLE}.genoprobs.npz \
    -g ${GBRS_DATA}/ref.genome_grid.GRCm39.tsv \
    -p ${GBRS_DATA}/ref.gene_pos.ordered_ensBuild_105.npz \
    -o ${WORK_DIR}/genome_output/${SAMPLE}.interpolated.genoprobs.npz

# Clean up temp directory
cd ${WORK_DIR}
rm -rf ${TEMP_DIR}
