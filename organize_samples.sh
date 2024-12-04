#!/bin/bash

# Loop through each unique sample name
for sample in $(ls *_R1_concatenated.fastq.gz | sed 's/_R1_concatenated.fastq.gz//' | sort | uniq)
do
    # Create a directory for each sample
    mkdir "$sample"
    
    # Move the corresponding R1 and R2 files into the sample directory
    mv "${sample}_R1_concatenated.fastq.gz" "$sample/"
    mv "${sample}_R2_concatenated.fastq.gz" "$sample/"
done
