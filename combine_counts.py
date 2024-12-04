#!/usr/bin/env python3
import pandas as pd
import csv
import glob
import os

def transform_counts_table(input_file, sample_name):
    # Define strain mapping for alignment columns
    strain_map = {
        'aln_A': 'AJ',
        'aln_B': 'C57',
        'aln_C': '129',
        'aln_D': 'NOD',
        'aln_E': 'NZO',
        'aln_F': 'CAST',
        'aln_G': 'PWK',
        'aln_H': 'WSB'
    }
    
    transformed_data = {}  # Dictionary to store transformed data
    
    with open(input_file, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)  # Skip header
        
        for row in reader:
            transcript_id = row[0]
            # Create strain-specific transcript IDs
            for strain_col, strain in strain_map.items():
                strain_transcript = f"{transcript_id}_{strain}"
                count_idx = header.index(strain_col)
                count = float(row[count_idx]) if row[count_idx] else 0.0
                
                if strain_transcript not in transformed_data:
                    transformed_data[strain_transcript] = {}
                
                transformed_data[strain_transcript][sample_name] = count
    
    return transformed_data

def combine_all_samples():
    # Get all processed.csv files in working directory
    file_list = sorted(glob.glob("*_processed.csv"))  # Sort to ensure consistent order
    
    # First get all samples and transcript-strains
    all_samples = []
    all_transcript_strains = set()
    
    print("First pass: collecting all samples and transcripts...")
    for file_path in file_list:
        sample_name = os.path.basename(file_path).replace('_processed.csv', '')
        all_samples.append(sample_name)
        
        data = transform_counts_table(file_path, sample_name)
        all_transcript_strains.update(data.keys())
    
    # Initialize the final matrix with zeros
    matrix = pd.DataFrame(0.0, 
                         index=sorted(list(all_transcript_strains)),
                         columns=all_samples)
    
    # Fill in the counts
    print("Second pass: filling in counts...")
    for i, file_path in enumerate(file_list):
        sample_name = os.path.basename(file_path).replace('_processed.csv', '')
        print(f"Processing {sample_name}...")
        
        data = transform_counts_table(file_path, sample_name)
        
        # Fill in counts for this sample
        for transcript_strain, counts in data.items():
            matrix.loc[transcript_strain, sample_name] = counts[sample_name]
    
    print("Saving matrix...")
    matrix.to_csv("complete_counts_matrix.csv", float_format='%.1f')
    
    print(f"Final matrix shape: {matrix.shape}")
    return matrix

if __name__ == "__main__":
    final_counts = combine_all_samples()