#!/bin/bash
#SBATCH --job-name=bam_sort
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 08:00:00
#SBATCH --mem=32000
#SBATCH --array=1-246 
#SBATCH --output=slurmout/bam_sort/bamsort_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/bam_sort/bamsort_%A_%a.err # File to which STDERR will be written
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dpbudke@ucdavis.edu

# Load required modules (modify as needed for your system)
module load samtools

# Get the sample ID from the samples.txt file based on SLURM array index
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

# Set directories
ALIGN_DIR="alignment_output"
OUT_DIR="sorted_bams"

# Create output directory if it doesn't exist
mkdir -p ${OUT_DIR}

# Sort both R1 and R2 BAM files
for READ in R1 R2; do
    INPUT_BAM="${ALIGN_DIR}/${SAMPLE}_${READ}.bam"
    OUTPUT_BAM="${OUT_DIR}/${SAMPLE}_${READ}.sorted.bam"

    echo "Processing ${INPUT_BAM}"
    
    samtools sort \
        -o ${OUTPUT_BAM} \
        ${INPUT_BAM}
    
    # Index the sorted BAM file
    samtools index ${OUTPUT_BAM}
done
