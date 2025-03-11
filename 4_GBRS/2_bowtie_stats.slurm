#!/bin/bash
#SBATCH --job-name=bam_stats
#SBATCH --output=slurmout/bam_stats/_%A_%a.out
#SBATCH --error=slurmout/bam_stats/_%A_%a.err
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 08:00:00
#SBATCH --mem=32G
#SBATCH --array=1-246

# Load necessary modules
module load samtools || echo "WARNING: Failed to load samtools module"

# Directory containing BAM files
BAM_DIR="/90daydata/do2_projects/CCPups_redo/EMASE/alignment_output"
FASTQ_DIR="/90daydata/do2_projects/CCPups_redo/01-HTS_Preproc"
OUTPUT_DIR="./bam_stats_results"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Get the sample ID from the array task ID
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)
echo "Processing sample: ${SAMPLE}"

# Output file for this sample
SAMPLE_OUTPUT="${OUTPUT_DIR}/${SAMPLE}_stats.txt"

# Initialize output file
echo -e "Sample\tTotal_R1_Reads\tTotal_R2_Reads\tMapped_R1_Reads\tMapped_R2_Reads\tTotal_Mapped_Reads\tUnmapped_Reads" > $SAMPLE_OUTPUT

# Process R1
r1_bam="${BAM_DIR}/${SAMPLE}_R1.bam"
echo "Looking for BAM file: $r1_bam"

if [ -f "$r1_bam" ]; then
    echo "Found R1 BAM file"
    
    # Get total reads from R1 fastq
    r1_fastq_path=$(find ${FASTQ_DIR}/${SAMPLE} -name "*R1*.fastq.gz" | head -1)
    echo "R1 fastq path: $r1_fastq_path"
    
    if [ -z "$r1_fastq_path" ] || [ ! -f "$r1_fastq_path" ]; then
        echo "Warning: R1 fastq file not found!"
        total_r1_reads=0
    else
        echo "Counting R1 reads..."
        total_r1_reads=$(zcat $r1_fastq_path | awk 'NR%4==1' | wc -l)
        echo "Total R1 reads: $total_r1_reads"
    fi
    
    # Count uniquely mapped reads by counting unique read IDs
    echo "Counting uniquely mapped R1 reads..."
    r1_mapped=$(samtools view -F 4 $r1_bam | cut -f1 | sort -u | wc -l)
    echo "Uniquely mapped R1 reads: $r1_mapped"
else
    echo "R1 BAM file not found"
    r1_mapped=0
    total_r1_reads=0
fi

# Process R2
r2_bam="${BAM_DIR}/${SAMPLE}_R2.bam"
echo "Looking for BAM file: $r2_bam"

if [ -f "$r2_bam" ]; then
    echo "Found R2 BAM file"
    
    # Get total reads from R2 fastq
    r2_fastq_path=$(find ${FASTQ_DIR}/${SAMPLE} -name "*R2*.fastq.gz" | head -1)
    echo "R2 fastq path: $r2_fastq_path"
    
    if [ -z "$r2_fastq_path" ] || [ ! -f "$r2_fastq_path" ]; then
        echo "Warning: R2 fastq file not found!"
        total_r2_reads=0
    else
        echo "Counting R2 reads..."
        total_r2_reads=$(zcat $r2_fastq_path | awk 'NR%4==1' | wc -l)
        echo "Total R2 reads: $total_r2_reads"
    fi
    
    # Count uniquely mapped reads by counting unique read IDs
    echo "Counting uniquely mapped R2 reads..."
    r2_mapped=$(samtools view -F 4 $r2_bam | cut -f1 | sort -u | wc -l)
    echo "Uniquely mapped R2 reads: $r2_mapped"
else
    echo "R2 BAM file not found"
    r2_mapped=0
    total_r2_reads=0
fi

# Calculate total mapped and unmapped reads
total_mapped=$((r1_mapped + r2_mapped))
total_reads=$((total_r1_reads + total_r2_reads))
unmapped=$((total_reads - total_mapped))

echo "Total mapped reads: $total_mapped"
echo "Total unmapped reads: $unmapped"

# Write results to output file
echo -e "${SAMPLE}\t${total_r1_reads}\t${total_r2_reads}\t${r1_mapped}\t${r2_mapped}\t${total_mapped}\t${unmapped}" >> $SAMPLE_OUTPUT

echo "Statistics for sample ${SAMPLE} written to ${SAMPLE_OUTPUT}"
