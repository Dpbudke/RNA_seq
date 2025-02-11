#!/bin/bash

# Define directories
LANE1_DIR="/project/do2_projects/AFRI/CCPups_redo/Raw_data/Run_1/Data/vdp8jyst5/Unaligned/Project_BBDB_Nova1118_DTSA968"
LANE2_DIR="/project/do2_projects/AFRI/CCPups_redo/Raw_data/Run_2/Data/6k6wciftf9/Un_DTSA968/Project_BBDB_Nova1118P_Budke"
OUTPUT_DIR="./concatenated_reads"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Function to concatenate files
concatenate_files() {
    local sample=$1
    local read=$2
    
    lane1_file=$(find "${LANE1_DIR}" -maxdepth 1 -name "${sample}*_L008_${read}_001.fastq.gz" -print -quit)
    lane2_file=$(find "${LANE2_DIR}" -maxdepth 1 -name "${sample}*_L004_${read}_001.fastq.gz" -print -quit)
    output_file="${OUTPUT_DIR}/${sample}${read}_concatenated.fastq.gz"
    
    if [[ -f "$lane1_file" && -f "$lane2_file" ]]; then
        echo "Concatenating ${read} files for sample ${sample}..."
        cat "$lane1_file" "$lane2_file" > "$output_file"
        echo "Concatenation complete: $output_file"
    else
        echo "Error: Missing input file(s) for sample ${sample}, ${read}"
        [[ ! -f "$lane1_file" ]] && echo "Lane 1 file not found: $lane1_file"
        [[ ! -f "$lane2_file" ]] && echo "Lane 2 file not found: $lane2_file"
    fi
}

# Main loop
echo "Searching for files in: ${LANE1_DIR}"
mapfile -t sample_files < <(find "${LANE1_DIR}" -maxdepth 1 -name "*_L008_R1_001.fastq.gz")

if [ ${#sample_files[@]} -eq 0 ]; then
    echo "No files found matching the pattern in ${LANE1_DIR}"
    exit 1
fi

for sample_file in "${sample_files[@]}"; do
    # Extract sample name, taking only the part before the first underscore
    sample=$(basename "$sample_file" | cut -d'_' -f1)
    
    concatenate_files "${sample}_" "R1"
    concatenate_files "${sample}_" "R2"
done

echo "All concatenations complete."
