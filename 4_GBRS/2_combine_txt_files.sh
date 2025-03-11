#!/bin/bash
# Simple script to combine all TXT files in bam_stats_results directory

# Input directory and output file
INPUT_DIR="./bam_stats_results"
OUTPUT_FILE="summary_bowtie_alignments.txt"

echo "Combining all files from ${INPUT_DIR} into ${OUTPUT_FILE}..."

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
  echo "Error: Directory $INPUT_DIR not found!"
  exit 1
fi

# Get the first file to use its header
FIRST_FILE=$(ls ${INPUT_DIR}/*_stats.txt | head -1)
if [ -z "$FIRST_FILE" ]; then
  echo "Error: No stats files found in $INPUT_DIR"
  exit 1
fi

# Write header from first file to output
head -n 1 "$FIRST_FILE" > "$OUTPUT_FILE"
echo "Added header from $(basename "$FIRST_FILE")"

# Append data from each file (skipping the header)
for file in ${INPUT_DIR}/*_stats.txt; do
  filename=$(basename "$file")
  echo "Processing $filename..."
  tail -n +2 "$file" >> "$OUTPUT_FILE"
done

echo "Done! All results combined into $OUTPUT_FILE"
echo "Found $(grep -c -v "Sample" $OUTPUT_FILE) samples in the output file"
