import os
import pandas as pd
import glob
from collections import defaultdict

def process_files(base_dir):
    # Create output directory if it doesn't exist
    output_dir = os.path.join(base_dir, 'Final_counts')
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize dictionaries to store data
    counts_dict = defaultdict(dict)
    haplotypes_dict = defaultdict(dict)
    
    # Get list of all sample directories
    pattern = os.path.join(base_dir, 'quantify_ASE', '*', '*_gbrs.quantified.diploid.genes.expected_read_counts')
    files = glob.glob(pattern)
    
    # Process each file
    for file_path in files:
        # Extract sample name from directory
        sample_name = os.path.basename(os.path.dirname(file_path))
        print(f"Processing sample: {sample_name}")
        
        # Read the file with explicit tab delimiter
        try:
            df = pd.read_csv(file_path, sep='\t', index_col=False)
            
            # Add data to dictionaries
            for _, row in df.iterrows():
                counts_dict[row['locus']][sample_name] = row['total']
                haplotypes_dict[row['locus']][sample_name] = row['notes']
                
        except Exception as e:
            print(f"Error processing {sample_name}: {str(e)}")
    
    # Convert dictionaries to DataFrames
    counts_df = pd.DataFrame.from_dict(counts_dict, orient='index')
    haplotypes_df = pd.DataFrame.from_dict(haplotypes_dict, orient='index')
    
    # Sort columns (sample names) for consistency
    counts_df.sort_index(axis=1, inplace=True)
    haplotypes_df.sort_index(axis=1, inplace=True)
    
    # Save the results
    counts_output = os.path.join(output_dir, 'combined_counts.csv')
    haplotypes_output = os.path.join(output_dir, 'combined_haplotypes.csv')
    
    counts_df.to_csv(counts_output)
    haplotypes_df.to_csv(haplotypes_output)
    
    # Print summary statistics
    print(f"\nProcessing complete!")
    print(f"Counts table saved to: {counts_output}")
    print(f"Haplotypes table saved to: {haplotypes_output}")
    print(f"\nCounts table shape: {counts_df.shape}")
    print(f"Haplotypes table shape: {haplotypes_df.shape}")
    
    # Print first few rows of each table for verification
    print("\nFirst few rows of counts table:")
    print(counts_df.head())
    print("\nFirst few rows of haplotypes table:")
    print(haplotypes_df.head())
    
    return len(files)

if __name__ == "__main__":
    base_dir = "/90daydata/do2_projects/CCPups_redo/EMASE"
    num_samples = process_files(base_dir)
    print(f"\nTotal number of samples processed: {num_samples}")