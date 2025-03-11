#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
import glob
from collections import defaultdict

def read_quantification_file(file_path):
    """
    Read an EMASE quantification file and return a Series with gene counts
    """
    try:
        # Read the file with explicit tab delimiter
        df = pd.read_csv(file_path, sep='\t', index_col=False)
        # Create a Series with gene IDs as index and total counts as values
        return pd.Series(data=df['total'].values, index=df['locus'].values, 
                         name=os.path.basename(os.path.dirname(file_path)))
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return pd.Series(name=os.path.basename(os.path.dirname(file_path)))

def compare_counts(emase_dir, counts_file):
    """
    Compare individual EMASE quantification files with the consolidated counts table
    """
    # Read the consolidated counts table
    counts_df = pd.read_csv(counts_file, index_col=0)
    
    # Dictionary to store results
    results = {}
    
    # Get all sample files
    pattern = os.path.join(emase_dir, 'quantify_ASE', '*', '*_gbrs.quantified.diploid.genes.expected_read_counts')
    sample_files = glob.glob(pattern)
    
    # Extract sample names from directory paths
    sample_dirs = [os.path.basename(os.path.dirname(file)) for file in sample_files]
    sample_to_file = {os.path.basename(os.path.dirname(file)): file for file in sample_files}
    
    print(f"\nChecking {len(sample_dirs)} samples...")
    
    # Check if all samples in the counts table are present in the directories
    missing_samples = set(counts_df.columns) - set(sample_dirs)
    extra_samples = set(sample_dirs) - set(counts_df.columns)
    
    if missing_samples:
        print(f"\nWARNING: The following samples are in the counts table but not in {emase_dir}:")
        print(", ".join(sorted(missing_samples)))
    
    if extra_samples:
        print(f"\nWARNING: The following samples are in {emase_dir} but not in the counts table:")
        print(", ".join(sorted(extra_samples)))
    
    # Only process samples that exist in both places
    samples_to_check = set(counts_df.columns) & set(sample_dirs)
    
    for sample in sorted(samples_to_check):
        print(f"\nChecking sample: {sample}")
        # Get path to quantification file
        quant_file = sample_to_file.get(sample)
        
        if not quant_file or not os.path.exists(quant_file):
            results[sample] = {
                'status': 'ERROR',
                'message': f'File not found: {quant_file}'
            }
            continue
            
        try:
            # Read individual EMASE counts
            emase_counts = read_quantification_file(quant_file)
            
            # Get corresponding column from counts table
            consolidated_counts = counts_df[sample]
            
            # Align the indices (genes) of both series
            emase_counts, consolidated_counts = emase_counts.align(consolidated_counts, join='outer')
            
            # Compare the counts - need to handle floating point comparisons
            # Convert to numpy arrays for comparison
            emase_array = np.nan_to_num(emase_counts.values, nan=-999)
            consolidated_array = np.nan_to_num(consolidated_counts.values, nan=-999)
            
            # Check if the arrays are almost equal (handle floating point issues)
            is_close = np.allclose(
                emase_array,
                consolidated_array,
                rtol=1e-05,  # Relative tolerance
                atol=1e-08,  # Absolute tolerance
                equal_nan=True
            )
            
            # Calculate statistics
            diff = emase_counts - consolidated_counts
            max_diff = diff.abs().max() if not diff.isna().all() else 0
            num_differences = (~np.isclose(
                np.nan_to_num(emase_counts.values),
                np.nan_to_num(consolidated_counts.values),
                rtol=1e-05, atol=1e-08, equal_nan=True
            )).sum()
            
            # Store detailed results
            results[sample] = {
                'status': 'PASS' if is_close else 'FAIL',
                'matching': is_close,
                'num_differences': num_differences,
                'max_difference': max_diff,
                'missing_in_emase': consolidated_counts.index[consolidated_counts.notna() & emase_counts.isna()].tolist(),
                'missing_in_counts': emase_counts.index[emase_counts.notna() & consolidated_counts.isna()].tolist(),
                'first_few_diffs': [] if is_close else [
                    (idx, emase_counts[idx], consolidated_counts[idx])
                    for idx in diff[abs(diff) > 1e-05].head().index
                ]
            }
            
        except Exception as e:
            results[sample] = {
                'status': 'ERROR',
                'message': str(e)
            }
    
    return results

def compare_haplotypes(emase_dir, haplotypes_file):
    """
    Compare individual EMASE haplotype information with the consolidated haplotypes table
    """
    # Read the consolidated haplotypes table
    haplotypes_df = pd.read_csv(haplotypes_file, index_col=0)
    
    # Dictionary to store results
    results = {}
    
    # Get all sample files
    pattern = os.path.join(emase_dir, 'quantify_ASE', '*', '*_gbrs.quantified.diploid.genes.expected_read_counts')
    sample_files = glob.glob(pattern)
    
    # Extract sample names from directory paths
    sample_dirs = [os.path.basename(os.path.dirname(file)) for file in sample_files]
    sample_to_file = {os.path.basename(os.path.dirname(file)): file for file in sample_files}
    
    print(f"\nChecking haplotype information for {len(sample_dirs)} samples...")
    
    # Check if all samples in the haplotypes table are present in the directories
    missing_samples = set(haplotypes_df.columns) - set(sample_dirs)
    extra_samples = set(sample_dirs) - set(haplotypes_df.columns)
    
    if missing_samples:
        print(f"\nWARNING: The following samples are in the haplotypes table but not in {emase_dir}:")
        print(", ".join(sorted(missing_samples)))
    
    if extra_samples:
        print(f"\nWARNING: The following samples are in {emase_dir} but not in the haplotypes table:")
        print(", ".join(sorted(extra_samples)))
    
    # Only process samples that exist in both places
    samples_to_check = set(haplotypes_df.columns) & set(sample_dirs)
    
    for sample in sorted(samples_to_check):
        print(f"\nChecking haplotypes for sample: {sample}")
        # Get path to quantification file
        quant_file = sample_to_file.get(sample)
        
        if not quant_file or not os.path.exists(quant_file):
            results[sample] = {
                'status': 'ERROR',
                'message': f'File not found: {quant_file}'
            }
            continue
            
        try:
            # Read the original file with explicit tab delimiter
            emase_df = pd.read_csv(quant_file, sep='\t', index_col=False)
            
            # Create a Series with gene IDs as index and notes (haplotypes) as values
            emase_haplotypes = pd.Series(data=emase_df['notes'].values, index=emase_df['locus'].values, 
                                       name=sample)
            
            # Get corresponding column from haplotypes table
            consolidated_haplotypes = haplotypes_df[sample]
            
            # Align the indices (genes) of both series
            emase_haplotypes, consolidated_haplotypes = emase_haplotypes.align(consolidated_haplotypes, join='outer')
            
            # Count mismatches - FIX: properly handle NaN values
            # First, identify positions where both series have NaN
            both_nan = emase_haplotypes.isna() & consolidated_haplotypes.isna()
            # Then do the comparison, but consider NaN==NaN to be True (not different)
            diff = (emase_haplotypes != consolidated_haplotypes) & ~both_nan
            num_differences = diff.sum()
            
            # Store detailed results
            results[sample] = {
                'status': 'PASS' if num_differences == 0 else 'FAIL',
                'matching': num_differences == 0,
                'num_differences': num_differences,
                'missing_in_emase': consolidated_haplotypes.index[consolidated_haplotypes.notna() & emase_haplotypes.isna()].tolist(),
                'missing_in_counts': emase_haplotypes.index[emase_haplotypes.notna() & consolidated_haplotypes.isna()].tolist(),
                'first_few_diffs': [] if num_differences == 0 else [
                    (idx, emase_haplotypes[idx], consolidated_haplotypes[idx])
                    for idx in diff[diff].head().index
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
    base_dir = "/90daydata/do2_projects/CCPups_redo/EMASE"
    emase_dir = base_dir
    counts_file = os.path.join(base_dir, "Final_counts", "combined_counts.csv")
    haplotypes_file = os.path.join(base_dir, "Final_counts", "combined_haplotypes.csv")
    
    if not os.path.exists(counts_file):
        print(f"Error: Counts file not found: {counts_file}")
        return
        
    if not os.path.exists(haplotypes_file):
        print(f"Error: Haplotypes file not found: {haplotypes_file}")
        return
        
    if not os.path.exists(emase_dir):
        print(f"Error: EMASE directory not found: {emase_dir}")
        return
    
    # Run comparison for counts
    print("\n=== Checking Read Counts ===")
    counts_results = compare_counts(emase_dir, counts_file)
    
    # Run comparison for haplotypes
    print("\n=== Checking Haplotype Information ===")
    haplotypes_results = compare_haplotypes(emase_dir, haplotypes_file)
    
    # Print results summary for counts
    print("\nRead Counts QC Results Summary:")
    print("-" * 50)
    
    # Count total passes/fails for counts
    total_samples = len(counts_results)
    passes = sum(1 for r in counts_results.values() if r.get('status') == 'PASS')
    fails = sum(1 for r in counts_results.values() if r.get('status') == 'FAIL')
    errors = sum(1 for r in counts_results.values() if r.get('status') == 'ERROR')
    
    print(f"\nTotal samples checked: {total_samples}")
    print(f"Passed: {passes}")
    print(f"Failed: {fails}")
    print(f"Errors: {errors}")
    
    # Print detailed results for non-passing samples (counts)
    if fails > 0 or errors > 0:
        print("\nDetailed results for non-passing count samples:")
        print("-" * 50)
        for sample, result in counts_results.items():
            if result['status'] != 'PASS':
                print(f"\nSample: {sample}")
                if result['status'] == 'ERROR':
                    print(f"Error: {result['message']}")
                else:
                    print(f"Status: {result['status']}")
                    print(f"Number of differences: {result['num_differences']}")
                    print(f"Maximum difference: {result['max_difference']}")
                    if result['missing_in_emase']:
                        print(f"Genes missing in EMASE counts: {len(result['missing_in_emase'])}")
                    if result['missing_in_counts']:
                        print(f"Genes missing in consolidated counts: {len(result['missing_in_counts'])}")
                    if result['first_few_diffs']:
                        print("\nFirst few differences (Gene, EMASE count, Consolidated count):")
                        for gene, emase_count, consolidated_count in result['first_few_diffs']:
                            print(f"{gene}: {emase_count} vs {consolidated_count}")
    
    # Print results summary for haplotypes
    print("\nHaplotype QC Results Summary:")
    print("-" * 50)
    
    # Count total passes/fails for haplotypes
    total_samples = len(haplotypes_results)
    passes = sum(1 for r in haplotypes_results.values() if r.get('status') == 'PASS')
    fails = sum(1 for r in haplotypes_results.values() if r.get('status') == 'FAIL')
    errors = sum(1 for r in haplotypes_results.values() if r.get('status') == 'ERROR')
    
    print(f"\nTotal samples checked: {total_samples}")
    print(f"Passed: {passes}")
    print(f"Failed: {fails}")
    print(f"Errors: {errors}")
    
    # Print detailed results for non-passing samples (haplotypes)
    if fails > 0 or errors > 0:
        print("\nDetailed results for non-passing haplotype samples:")
        print("-" * 50)
        for sample, result in haplotypes_results.items():
            if result['status'] != 'PASS':
                print(f"\nSample: {sample}")
                if result['status'] == 'ERROR':
                    print(f"Error: {result['message']}")
                else:
                    print(f"Status: {result['status']}")
                    print(f"Number of differences: {result['num_differences']}")
                    if result['missing_in_emase']:
                        print(f"Genes missing in EMASE haplotypes: {len(result['missing_in_emase'])}")
                    if result['missing_in_counts']:
                        print(f"Genes missing in consolidated haplotypes: {len(result['missing_in_counts'])}")
                    if result['first_few_diffs']:
                        print("\nFirst few differences (Gene, EMASE haplotype, Consolidated haplotype):")
                        for gene, emase_haplotype, consolidated_haplotype in result['first_few_diffs']:
                            print(f"{gene}: '{emase_haplotype}' vs '{consolidated_haplotype}'")

if __name__ == "__main__":
    main()