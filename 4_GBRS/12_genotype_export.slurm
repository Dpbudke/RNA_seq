#!/bin/bash

#SBATCH --job-name=gbrs_export
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 08:00:00
#SBATCH --mem=32G
#SBATCH --array=1-246
#SBATCH --output=slurmout/export/_%A_%a.out
#SBATCH --error=slurmout/export/_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dpbudke@ucdavis.edu

# Set the working directory
cd /90daydata/do2_projects/CCPups_redo/EMASE/

# Set GBRS_DATA environment variable
export GBRS_DATA="/90daydata/do2_projects/CCPups_redo/EMASE/GBRS_DATA"

# Create output directory if it doesn't exist
mkdir -p genotype_probs

# Get list of sample files and extract the current sample based on array task ID
sample_files=(genome_output/*.interpolated.genoprobs.npz)
current_file="${sample_files[$SLURM_ARRAY_TASK_ID-1]}"

# Extract sample name from the file path
sample_name=$(basename "$current_file" .interpolated.genoprobs.npz)

# Run gbrs export command
gbrs export \
    -i "${current_file}" \
    -s A,B,C,D,E,F,G,H \
    -g ${GBRS_DATA}/ref.genome_grid.GRCm39.tsv \
    -o genotype_probs/${sample_name}.gbrs.interpolated.genoprobs.tsv
