#!/bin/bash
#SBATCH --job-name=emase_alignment # Job name
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 08:00:00
#SBATCH --mem=32000
#SBATCH --array=1-246
#SBATCH --output=slurmout/emase_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/emase_%A_%a.err # File to which STDERR will be written
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dpbudke@ucdavis.edu

# Load required modules
module load bowtie
module load samtools

# Set main directories
READS_DIR="/90daydata/do2_projects/CCPups_redo/01-HTS_Preproc"
OUTPUT_DIR="/90daydata/do2_projects/CCPups_redo/EMASE/alignment_output"
BOWTIE_INDEX="/90daydata/do2_projects/CCPups_redo/EMASE/GBRS_DATA/bowtie/bowtie.transcripts"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}
mkdir -p ${OUTPUT_DIR}/logs

# Create a log file
LOGFILE="${OUTPUT_DIR}/logs/alignment_$(date +%Y%m%d_%H%M%S).log"

# Function to log messages
log_message() {
    echo "$(date +"%Y-%m-%d %H:%M:%S") - $1" | tee -a ${LOGFILE}
}

# Read sample name from samples.txt
sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" samples.txt`
log_message "Starting processing of sample: ${sample}"

# Construct sample directory path
sample_dir="${READS_DIR}/${sample}"

# Check if the sample directory exists
if [ ! -d "${sample_dir}" ]; then
    log_message "ERROR: Sample directory does not exist: ${sample_dir}"
    exit 1
fi

# Process R1 files
R1_files=(${sample_dir}/*R1*.fastq.gz)
if [ ${#R1_files[@]} -gt 0 ]; then
    read_file="${R1_files[0]}"
    log_message "Processing R1 file: ${read_file}"
    bam_file="${OUTPUT_DIR}/${sample}_R1.bam"
    
    # Run bowtie alignment
    log_message "Starting bowtie alignment for R1"
    bowtie -q -a --best --strata --sam -v 3 ${BOWTIE_INDEX} ${read_file} 2>>${LOGFILE} | samtools view -bS - > ${bam_file}
    
    if [ -s ${bam_file} ]; then
        log_message "Successfully created ${bam_file}"
    else
        log_message "ERROR: Failed to create or empty file ${bam_file}"
    fi
else
    log_message "WARNING: No R1 files found in ${sample_dir}"
fi

# Process R2 files
R2_files=(${sample_dir}/*R2*.fastq.gz)
if [ ${#R2_files[@]} -gt 0 ]; then
    read_file="${R2_files[0]}"
    log_message "Processing R2 file: ${read_file}"
    bam_file="${OUTPUT_DIR}/${sample}_R2.bam"
    
    # Run bowtie alignment
    log_message "Starting bowtie alignment for R2"
    bowtie -q -a --best --strata --sam -v 3 ${BOWTIE_INDEX} ${read_file} 2>>${LOGFILE} | samtools view -bS - > ${bam_file}
    
    if [ -s ${bam_file} ]; then
        log_message "Successfully created ${bam_file}"
    else
        log_message "ERROR: Failed to create or empty file ${bam_file}"
    fi
else
    log_message "WARNING: No R2 files found in ${sample_dir}"
fi

log_message "Finished processing sample: ${sample}"

