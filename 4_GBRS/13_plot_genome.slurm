#!/bin/bash
#SBATCH --job-name=gbrs_plot
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 08:00:00
#SBATCH --mem=32G
#SBATCH --array=1-246
#SBATCH --output=slurmout/plot/_%A_%a.out
#SBATCH --error=slurmout/plot/_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dpbudke@ucdavis.edu

# Set the working directory
cd /90daydata/do2_projects/CCPups_redo/EMASE/

# Set GBRS_DATA environment variable
export GBRS_DATA="/90daydata/do2_projects/CCPups_redo/EMASE/GBRS_DATA"

# Create output directory if it doesn't exist
mkdir -p genome_plots

# Option 1: Process files from genome_output directory
for INPUT_FILE in genome_output/*.interpolated.genoprobs.npz; do
    # Extract sample name from filename
    FILENAME=$(basename "$INPUT_FILE")
    SAMPLE_ID="${FILENAME%%.*}"
    
    echo "Processing sample: $SAMPLE_ID"
    
    # Run gbrs plot command
    gbrs plot \
        -i "$INPUT_FILE" \
        -o "genome_plots/${SAMPLE_ID}.plotted.genome.pdf" \
        -n "$SAMPLE_ID"
done
