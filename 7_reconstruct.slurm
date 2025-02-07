#!/bin/bash

#SBATCH --job-name=GBRS_reconstruct
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 08:00:00
#SBATCH --mem=32G
#SBATCH --array=1-246
#SBATCH --output=slurmout/reconstruct/recon_%A_%a.out
#SBATCH --error=slurmout/reconstruct/recon_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dpbudke@ucdavis.edu

# Exit on error
set -e 

# Function to print timestamp with message
log_time() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Directory setup
WORKING_DIR="/90daydata/do2_projects/CCPups_redo/EMASE"
export GBRS_DATA="${WORKING_DIR}/GBRS_DATA"
QUANTIFY_OUTPUT="${WORKING_DIR}/quantify_output"

# Create logs directory
mkdir -p logs

# Determine if the current sample is male or female and get sample ID
TOTAL_SAMPLES=$(wc -l < samples.txt)
if [ ${SLURM_ARRAY_TASK_ID} -gt ${TOTAL_SAMPLES} ]; then
    log_time "Error: Task ID exceeds number of samples"
    exit 1
fi

SAMPLE_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)
log_time "Processing sample: ${SAMPLE_ID}"

# Determine sex from metadata
SEX=$(awk -F',' -v sample="$SAMPLE_ID" '$3 ~ sample {print $6}' Data4R_AllPups.csv | tr -d ' ' | head -n 2 | tail -n 1 | cut -c1 | tr '[:lower:]' '[:upper:]')

if [ "$SEX" = "F" ]; then
    TRANPROB_FILE="${GBRS_DATA}/tranprob.CC.G20.F.npz"
    log_time "Sample is female, using female transition probabilities"
elif [ "$SEX" = "M" ]; then
    TRANPROB_FILE="${GBRS_DATA}/tranprob.CC.G20.M.npz"
    log_time "Sample is male, using male transition probabilities"
else
    log_time "Error: Unable to determine sex for sample ${SAMPLE_ID}"
    exit 1
fi

# Setup paths
OUTPUT_DIR="${WORKING_DIR}/reconstruct_output/${SAMPLE_ID}"
TPM_FILE="${QUANTIFY_OUTPUT}/${SAMPLE_ID}/${SAMPLE_ID}.multiway.genes.tpm"
GENE_POS_FILE="${GBRS_DATA}/ref.gene_pos.ordered_ensBuild_105.npz"
OUTBASE="${OUTPUT_DIR}/${SAMPLE_ID}"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Create symlinks if they don't exist
if [ ! -f "${GBRS_DATA}/avecs.npz" ]; then
    log_time "Creating symlink for avecs.npz..."
    ln -sf "${GBRS_DATA}/gbrs_emissions_all_tissues.avecs.npz" "${GBRS_DATA}/avecs.npz"
fi

if [ ! -f "${GBRS_DATA}/ref.fa.fai" ]; then
    log_time "Creating symlink for ref.fa.fai..."
    ln -sf "${GBRS_DATA}/Mus_musculus.GRCm39.dna.primary_assembly.fa.fai" "${GBRS_DATA}/ref.fa.fai"
fi

# Verify input files
log_time "Verifying input files..."
for file in "$TPM_FILE" "$TRANPROB_FILE" "$GENE_POS_FILE"; do
    if [ ! -r "$file" ]; then
        log_time "ERROR: Cannot read file: $file"
        ls -l "$file"
        exit 1
    fi
    size=$(ls -lh "$file" | awk '{print $5}')
    log_time "Verified: $file (Size: $size)"
done

# Run GBRS reconstruction
log_time "Starting GBRS reconstruction..."
log_time "Command to be executed:"
echo "gbrs reconstruct \\"
echo "    -e \"${TPM_FILE}\" \\"
echo "    -t \"${TRANPROB_FILE}\" \\"
echo "    -g \"${GENE_POS_FILE}\" \\"
echo "    -o \"${OUTBASE}\""

# Run the command
gbrs reconstruct \
    -e "${TPM_FILE}" \
    -t "${TRANPROB_FILE}" \
    -g "${GENE_POS_FILE}" \
    -o "${OUTBASE}" 2>&1

EXIT_STATUS=$?
log_time "GBRS reconstruction process finished with exit status: ${EXIT_STATUS}"

# Check for output files
for ext in ".genotypes.tsv" ".genoprobs.npz" ".genotypes.npz"; do
    if [ -f "${OUTBASE}${ext}" ]; then
        size=$(ls -lh "${OUTBASE}${ext}" | awk '{print $5}')
        log_time "Found output file: ${OUTBASE}${ext} (Size: $size)"
    else
        log_time "Missing expected output file: ${OUTBASE}${ext}"
        EXIT_STATUS=1
    fi
done

exit $EXIT_STATUS
