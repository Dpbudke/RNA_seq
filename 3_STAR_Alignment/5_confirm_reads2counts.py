#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np

def read_readspergene(file_path):
    """
    Read a STAR ReadsPerGene.out.tab file and return a Series with gene counts
    Skip exactly 4 lines (matching tail -n +5) and use the 1st and 4th columns
    """
    # Read the file, skipping exactly 4 rows (matching tail -n +5)
    df = pd.read_csv(file_path, sep='\t', skiprows=4, header=None, names=['gene_id', 'col2', 'col3', 'count'])
    # Create a Series with gene IDs as index and counts as values
    return pd.Series(data=df['count'].values, index=df['gene_id'].values, name=os.path.basename(os.path.dirname(file_path)))

def compare_counts(star_dir, counts_file):
    """
    Compare individual STAR count files with the consolidated counts table
    Follows the exact same process used in the bash script that generated the counts
    """
    # Read the consolidated counts table
    counts_df = pd.read_csv(counts_file, sep='\t', index_col=0)
    
    # Dictionary to store results
    results = {}
    
    # Get all sample directories
    sample_dirs = [d for d in os.listdir(star_dir) 
                  if os.path.isdir(os.path.join(star_dir, d))]
    sample_dirs.sort()  # Sort for consistent output
    
    print(f"\nChecking {len(sample_dirs)} samples...")
    
    # Check if all samples in the counts table are present in the directories
    missing_samples = set(counts_df.columns) - set(sample_dirs)
    extra_samples = set(sample_dirs) - set(counts_df.columns)
    
    if missing_samples:
        print(f"\nWARNING: The following samples are in the counts table but not in {star_dir}:")
        print(", ".join(sorted(missing_samples)))
    
    if extra_samples:
        print(f"\nWARNING: The following samples are in {star_dir} but not in the counts table:")
        print(", ".join(sorted(extra_samples)))
    
    # Only process samples that exist in both places
    samples_to_check = set(counts_df.columns) & set(sample_dirs)
    
    for sample in sorted(samples_to_check):
        print(f"\nChecking sample: {sample}")
        # Construct path to ReadsPerGene file
        readspergene_file = os.path.join(star_dir, sample, f"{sample}_ReadsPerGene.out.tab")
        
        if not os.path.exists(readspergene_file):
            results[sample] = {
                'status': 'ERROR',
                'message': f'File not found: {readspergene_file}'
            }
            continue
            
        try:
            # Read individual STAR counts - only the 4th column after skipping 4 lines
            star_counts = read_readspergene(readspergene_file)
            
            # Get corresponding column from counts table
            consolidated_counts = counts_df[sample]
            
            # Align the indices (genes) of both series
            star_counts, consolidated_counts = star_counts.align(consolidated_counts, join='outer')
            
            # Compare the counts
            matching = np.array_equal(
                star_counts.fillna(-999).values,  # Use -999 as sentinel for missing values
                consolidated_counts.fillna(-999).values
            )
            
            # Calculate statistics
            diff = star_counts - consolidated_counts
            max_diff = diff.abs().max() if not diff.isna().all() else 0
            num_differences = (diff != 0).sum()
            
            # Store detailed results
            results[sample] = {
                'status': 'PASS' if matching else 'FAIL',
                'matching': matching,
                'num_differences': num_differences,
                'max_difference': max_diff,
                'missing_in_star': consolidated_counts.index[consolidated_counts.notna() & star_counts.isna()].tolist(),
                'missing_in_counts': star_counts.index[star_counts.notna() & consolidated_counts.isna()].tolist(),
                'first_few_diffs': [] if matching else [
                    (idx, star_counts[idx], consolidated_counts[idx])
                    for idx in diff[diff != 0].head().index
                ]
            }
            
        except Exception as e:
            results[sample] = {
                'status': 'ERROR',
                'message': str(e)
            }
    
    return results

def main():
    # Define paths
    star_dir = "02-STAR_alignment"
    counts_file = "03-Counts/GRCm39_AFRI_counts.txt"
    
    if not os.path.exists(counts_file):
        print(f"Error: Counts file not found: {counts_file}")
        return
        
    if not os.path.exists(star_dir):
        print(f"Error: STAR alignment directory not found: {star_dir}")
        return
    
    # Run comparison
    results = compare_counts(star_dir, counts_file)
    
    # Print results summary
    print("\nQC Results Summary:")
    print("-" * 50)
    
    # Count total passes/fails
    total_samples = len(results)
    passes = sum(1 for r in results.values() if r.get('status') == 'PASS')
    fails = sum(1 for r in results.values() if r.get('status') == 'FAIL')
    errors = sum(1 for r in results.values() if r.get('status') == 'ERROR')
    
    print(f"\nTotal samples checked: {total_samples}")
    print(f"Passed: {passes}")
    print(f"Failed: {fails}")
    print(f"Errors: {errors}")
    
    # Print detailed results for non-passing samples
    if fails > 0 or errors > 0:
        print("\nDetailed results for non-passing samples:")
        print("-" * 50)
        for sample, result in results.items():
            if result['status'] != 'PASS':
                print(f"\nSample: {sample}")
                if result['status'] == 'ERROR':
                    print(f"Error: {result['message']}")
                else:
                    print(f"Status: {result['status']}")
                    print(f"Number of differences: {result['num_differences']}")
                    print(f"Maximum difference: {result['max_difference']}")
                    if result['missing_in_star']:
                        print(f"Genes missing in STAR counts: {len(result['missing_in_star'])}")
                    if result['missing_in_counts']:
                        print(f"Genes missing in consolidated counts: {len(result['missing_in_counts'])}")
                    if result['first_few_diffs']:
                        print("\nFirst few differences (Gene, STAR count, Consolidated count):")
                        for gene, star_count, consolidated_count in result['first_few_diffs']:
                            print(f"{gene}: {star_count} vs {consolidated_count}")

if __name__ == "__main__":
    main()