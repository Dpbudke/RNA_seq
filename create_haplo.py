#!/usr/bin/env python3

import pandas as pd
import os
from pathlib import Path
import sys

def get_strain_map():
    """Create a mapping of exact strain names to their alignment columns."""
    strain_map = {
        'A/J': 'aln_A',
        'C57BL/6J': 'aln_B',
        '129S1/SvImJ': 'aln_C',
        'NOD/LtJ': 'aln_D',
        'NZO/HILtJ': 'aln_E',
        'CAST/EiJ': 'aln_F',
        'PWK/PhJ': 'aln_G',
        'WSB/EiJ': 'aln_H'
    }
    return strain_map

def unique_non_null(lst):
    """Return unique non-null values from a list while preserving order."""
    seen = set()
    return [x for x in lst if pd.notna(x) and not (x in seen or seen.add(x))]

def load_metadata(metadata_path):
    """Load and process metadata file from working directory."""
    try:
        # Read the CSV file
        df = pd.read_csv(metadata_path)
        
        # Strip whitespace from column names
        df.columns = df.columns.str.strip()
        
        print("\nAvailable columns in metadata file:")
        for col in df.columns:
            print(f"  - {col}")
        
        # Look for Animal ID column (after stripping whitespace)
        if 'Animal ID' not in df.columns:
            print("\nError: Could not find 'Animal ID' column")
            sys.exit(1)
            
        # Look for Cross column (after stripping whitespace)
        if 'Cross' not in df.columns:
            print("\nError: Could not find 'Cross' column")
            sys.exit(1)
            
        # Create a standardized DataFrame with renamed columns
        standardized_df = df.rename(columns={
            'Animal ID': 'Animal_ID',
            'Cross': 'Cross'
        })
        
        print(f"\nFound required columns:")
        print(f"  - Using 'Animal ID' for sample identification")
        print(f"  - Using 'Cross' for strain cross information")
        
        # Verify we have data
        if len(standardized_df) == 0:
            print("Error: Metadata file contains no data")
            sys.exit(1)
            
        return standardized_df
        
    except FileNotFoundError:
        print(f"Error: Metadata file '{metadata_path}' not found in current directory")
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"Error: Metadata file '{metadata_path}' is empty")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading metadata file: {str(e)}")
        sys.exit(1)

def load_haplotype_file(cc_strain):
    """Load haplotype determination file for a given CC strain from working directory."""
    filename = f"{cc_strain}_final.csv"
    try:
        return pd.read_csv(filename)
    except FileNotFoundError:
        print(f"Error: Haplotype file '{filename}' not found in current directory")
        return None

def get_haplotype_files_for_cross(cross):
    """Determine which haplotype files to use based on the cross."""
    strains = cross.split('x')
    return strains

def is_ccrix(cross):
    """Determine if a cross is CC-RIX by checking if the strains are different."""
    strains = cross.split('x')
    return len(strains) == 2 and strains[0] != strains[1]

def process_sample(sample_id, cross, counts_base_path):
    """Process a single sample's data."""
    cc_strains = get_haplotype_files_for_cross(cross)
    print(f"\nProcessing {cross} sample:")
    
    # Determine if this is a CC-RIX sample
    ccrix = is_ccrix(cross)
    if ccrix:
        print("CC-RIX sample detected")
    else:
        print("Standard cross sample detected")
    
    print(f"Looking for haplotype files: {[f'{strain}_final.csv' for strain in cc_strains]}")
    
    # Load all necessary haplotype files from working directory
    haplotype_dfs = []
    # For non-CC-RIX, we only need one copy of the haplotype file
    seen_strains = set()
    for strain in cc_strains:
        if not ccrix and strain in seen_strains:
            continue
        seen_strains.add(strain)
        
        hap_df = load_haplotype_file(strain)
        if hap_df is not None:
            hap_df['source_strain'] = strain
            haplotype_dfs.append(hap_df)
    
    if not haplotype_dfs:
        print(f"Error: No valid haplotype files found for sample {sample_id}")
        return None
    
    # Combine haplotype information if necessary
    if len(haplotype_dfs) > 1:
        combined_hap = pd.concat(haplotype_dfs, axis=0)
    else:
        combined_hap = haplotype_dfs[0]
    
    # Load counts file
    counts_file = Path(counts_base_path) / str(sample_id).strip() / f"{str(sample_id).strip()}.multiway.isoforms.alignment_counts"
    try:
        counts_df = pd.read_csv(counts_file, sep='\t')
    except FileNotFoundError:
        print(f"Error: Counts file not found at {counts_file}")
        return None
    except PermissionError:
        print(f"Error: Permission denied accessing counts file at {counts_file}")
        return None
    
    counts_df = counts_df.rename(columns={'locus': 'transcript_id'})
    
    # Merge haplotype and counts data
    merged_df = pd.merge(combined_hap, counts_df, on='transcript_id', how='inner')
    
    strain_map = get_strain_map()
    
    # Create the output DataFrame with basic columns
    final_df = pd.DataFrame()
    final_df['transcript_id'] = merged_df['transcript_id'].unique()
    
    # Initialize strain columns
    for strain_col in ['Primary_Strain', 'Secondary_Strain']:
        final_df[strain_col] = ''
    
    # Initialize all alignment columns with NA
    all_aln_cols = [col for col in merged_df.columns if col.startswith('aln_')]
    for col in all_aln_cols:
        final_df[col] = pd.NA
    
    # Process each transcript
    for transcript_id in final_df['transcript_id']:
        transcript_rows = merged_df[merged_df['transcript_id'] == transcript_id]
        
        # Collect all strains for this transcript
        primary_strains = []
        secondary_strains = []
        alignment_values = {}
        
        for _, row in transcript_rows.iterrows():
            if pd.notna(row['Primary_Strain']):
                primary_strains.append(row['Primary_Strain'])
                if row['Primary_Strain'] in strain_map:
                    aln_col = strain_map[row['Primary_Strain']]
                    if aln_col in row:
                        alignment_values[aln_col] = row[aln_col]
            
            if pd.notna(row['Secondary_Strain']):
                secondary_strains.append(row['Secondary_Strain'])
                if row['Secondary_Strain'] in strain_map:
                    aln_col = strain_map[row['Secondary_Strain']]
                    if aln_col in row:
                        alignment_values[aln_col] = row[aln_col]
        
        # Update final DataFrame
        idx = final_df[final_df['transcript_id'] == transcript_id].index[0]
        final_df.at[idx, 'Primary_Strain'] = '; '.join(filter(None, unique_non_null(primary_strains)))
        final_df.at[idx, 'Secondary_Strain'] = '; '.join(filter(None, unique_non_null(secondary_strains)))
        
        # Set only the relevant alignment values
        for aln_col, value in alignment_values.items():
            final_df.at[idx, aln_col] = value
    
    # Print debug information
    print(f"\nDebug information for sample {sample_id}:")
    print(f"Number of unique transcripts: {len(final_df)}")
    print(f"Number of transcripts with alignments: {len(final_df.dropna(subset=[col for col in final_df.columns if col.startswith('aln_')], how='all'))}")
    print("\nExample rows:")
    print(final_df.head())
    
    return final_df

def main():
    # Configuration
    metadata_path = "filtered_Data4R_AllPups.csv"
    counts_base_path = "/90daydata/do2_projects/CCPups_redo/EMASE/quantify_output/"
    output_dir = "processed_output"
    
    # Print current working directory for verification
    print(f"Working directory: {os.getcwd()}")
    
    # Verify existence of required files in working directory
    required_files = ["filtered_Data4R_AllPups.csv", "CC001_final.csv", "CC019_final.csv", "CC040_final.csv"]
    missing_files = [f for f in required_files if not os.path.exists(f)]
    if missing_files:
        print("Error: The following required files are missing from the working directory:")
        for file in missing_files:
            print(f"  - {file}")
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Load metadata
    print("Loading metadata file...")
    metadata_df = load_metadata(metadata_path)
    
    # Process each sample
    total_samples = len(metadata_df)
    for idx, row in metadata_df.iterrows():
        sample_id = row['Animal_ID'].strip()
        cross = row['Cross'].strip()
        
        print(f"Processing sample {sample_id} ({cross})... [{idx + 1}/{total_samples}]")
        
        result_df = process_sample(sample_id, cross, counts_base_path)
        
        if result_df is not None:
            # Save processed data
            output_path = Path(output_dir) / f"{sample_id}_processed.csv"
            result_df.to_csv(output_path, index=False)
            print(f"Successfully saved processed data to {output_path}")
        else:
            print(f"Failed to process sample {sample_id}")

if __name__ == "__main__":
    main()