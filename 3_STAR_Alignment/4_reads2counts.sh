#!/bin/bash

# Create directories
mkdir -p 03-Counts/tmp

# Create and save list of samples (directories in 02-STAR_alignment)
ls -1 02-STAR_alignment/ > 03-Counts/tmp/sample_list.txt

# Create count files using the sample list
while read sample; do
    echo "Processing ${sample}"
    cat 02-STAR_alignment/${sample}/${sample}_ReadsPerGene.out.tab | tail -n +5 | cut -f4 > 03-Counts/tmp/${sample}.count
done < 03-Counts/tmp/sample_list.txt

# Get gene IDs (using first sample from the list)
FIRST_SAMPLE=$(head -n 1 03-Counts/tmp/sample_list.txt)
tail -n +5 02-STAR_alignment/${FIRST_SAMPLE}/${FIRST_SAMPLE}_ReadsPerGene.out.tab | cut -f1 > 03-Counts/tmp/geneids.txt

# Create the paste command using the same sample order
PASTE_CMD="paste 03-Counts/tmp/geneids.txt"
while read sample; do
    PASTE_CMD="${PASTE_CMD} 03-Counts/tmp/${sample}.count"
done < 03-Counts/tmp/sample_list.txt

# Create header using the same sample order
echo -n "gene_id" > 03-Counts/tmp/header.txt
while read sample; do
    echo -n -e "\t${sample}" >> 03-Counts/tmp/header.txt
done < 03-Counts/tmp/sample_list.txt
echo "" >> 03-Counts/tmp/header.txt

# Execute paste command and combine with header
${PASTE_CMD} > 03-Counts/tmp/tmp.out
cat 03-Counts/tmp/header.txt 03-Counts/tmp/tmp.out > 03-Counts/GRCm39_AFRI_counts.txt

# Clean up
rm -rf 03-Counts/tmp

echo "Done! Created 03-Counts/GRCm39_AFRI_counts.txt"