#!/bin/bash

#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 08:00:00
#SBATCH --mem=32000
#SBATCH --array=1-246 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dpbudke@ucdavis.edu
#SBATCH --job-name=emase_common_align
#SBATCH --output=slurmout/common_align_merge/common_align_%A_%a.out
#SBATCH --error=slurmout/common_align_merge/common_align_%A_%a.err

module load h5py/3.8.0  

python - <<'EOF'
import h5py
import numpy as np
import os
import sys
import logging
from datetime import datetime
import gc

# Set working directory
os.chdir('/90daydata/do2_projects/CCPups_redo/EMASE')

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def load_emase_file(filename):
    """Load an EMASE HDF5 file and return alignment data organized by strain."""
    logging.info(f"Loading EMASE file: {filename}")
    try:
        with h5py.File(filename, 'r') as f:
            # Check for both lname and rname, prioritize lname
            if '/lname' in f:
                read_names = f['/lname'][:]
            elif '/rname' in f:
                read_names = f['/rname'][:]
            else:
                raise KeyError("Neither /lname nor /rname found in file")
                
            if isinstance(read_names[0], bytes):
                read_names = [name.decode('utf-8') for name in read_names]
            logging.info(f"Loaded {len(read_names)} read names from {filename}")

            alignments_by_strain = {}
            for i in range(8):  # h0 to h7 for each progenitor strain
                strain_group = f'/h{i}'
                if strain_group in f:
                    indices = f[f'{strain_group}/indices'][:]
                    indptr = f[f'{strain_group}/indptr'][:]
                    alignments_by_strain[f'strain_{i}'] = (indices, indptr)
                    logging.info(f"Loaded alignments for {strain_group}")
                    gc.collect()
            
            # Copy all attributes from the input file
            file_attrs = dict(f.attrs)
            return read_names, alignments_by_strain, file_attrs
    except Exception as e:
        logging.error(f"Error loading EMASE file {filename}: {str(e)}")
        raise

def save_emase_file(filename, common_reads, alignments_by_strain, original_attrs):
    """Save common alignments to a new EMASE HDF5 file."""
    logging.info(f"Saving to: {filename}")
    try:
        with h5py.File(filename, 'w') as f:
            # Save read names under both lname and rname for compatibility
            read_names_encoded = np.array([name.encode('utf-8') for name in common_reads])
            f.create_dataset('lname', data=read_names_encoded)  # Primary name for GBRS
            f.create_dataset('rname', data=read_names_encoded)  # Backup for compatibility
            del read_names_encoded
            gc.collect()
            
            # Save alignments for each strain
            for strain, (indices, indptr) in alignments_by_strain.items():
                strain_num = int(strain.split('_')[1])
                strain_group = f.create_group(f'h{strain_num}')
                strain_group.create_dataset('indices', data=indices)
                strain_group.create_dataset('indptr', data=indptr)
                gc.collect()
            
            # Copy original attributes
            for key, value in original_attrs.items():
                f.attrs[key] = value
                
            # Update or add new attributes
            f.attrs['creation_date'] = datetime.now().isoformat()
            f.attrs['number_of_reads'] = len(common_reads)
            
        logging.info(f"Successfully saved to {filename}")
    except Exception as e:
        logging.error(f"Error saving file {filename}: {str(e)}")
        raise

def get_read_alignments(indices, indptr, read_idx):
    """Extract alignments for a specific read from CSR format data."""
    start = indptr[read_idx]
    end = indptr[read_idx + 1]
    return indices[start:end]

def process_strain_alignments(strain, r1_alignments, r2_alignments, read_indices):
    """Process alignments for a single strain to reduce memory usage."""
    logging.info(f"Processing alignments for strain: {strain}")
    try:
        r1_indices, r1_indptr = r1_alignments[strain]
        r2_indices, r2_indptr = r2_alignments[strain]
        
        # Initialize output arrays
        all_indices = []
        indptr = [0]
        
        # Process each read
        for read_idx in read_indices:
            # Get alignments for both R1 and R2
            r1_aligns = get_read_alignments(r1_indices, r1_indptr, read_idx)
            r2_aligns = get_read_alignments(r2_indices, r2_indptr, read_idx)
            
            # Combine alignments
            combined = np.unique(np.concatenate([r1_aligns, r2_aligns]))
            all_indices.extend(combined)
            indptr.append(len(all_indices))
            
            # Garbage collection every 10000 reads
            if len(indptr) % 10000 == 0:
                gc.collect()
                logging.info(f"Processed {len(indptr)-1} reads for {strain}")
        
        return np.array(all_indices), np.array(indptr)
        
    except Exception as e:
        logging.error(f"Error processing strain {strain}: {str(e)}")
        raise

def find_common_alignments(r1_file, r2_file, output_file):
    logging.info(f"Starting alignment processing for files: {r1_file}, {r2_file}")
    r1_reads, r1_alignments, r1_attrs = load_emase_file(r1_file)
    r2_reads, r2_alignments, _ = load_emase_file(r2_file)
    logging.info("Successfully loaded R1 and R2 reads and alignments.")

    # Create index mappings for faster lookups
    r1_read_dict = {read: idx for idx, read in enumerate(r1_reads)}
    common_reads = []
    common_indices = []
    
    # Find common reads and their indices
    for read in r2_reads:
        if read in r1_read_dict:
            common_reads.append(read)
            common_indices.append(r1_read_dict[read])
    
    common_indices = np.array(common_indices)
    logging.info(f"Found {len(common_reads)} common reads")

    # Process each strain separately
    common_alignments = {}
    for strain in r1_alignments.keys():
        if strain in r2_alignments:
            indices, indptr = process_strain_alignments(
                strain, 
                r1_alignments, 
                r2_alignments, 
                common_indices
            )
            common_alignments[strain] = (indices, indptr)
            gc.collect()

    logging.info(f"Saving to output file: {output_file}")
    save_emase_file(output_file, common_reads, common_alignments, r1_attrs)
    
    return {
        'total_r1_reads': len(r1_reads),
        'total_r2_reads': len(r2_reads),
        'common_reads': len(common_reads)
    }

def main():
    task_id = int(os.environ.get('SLURM_ARRAY_TASK_ID', 0))
    logging.info(f"Starting processing for SLURM array task ID: {task_id}")

    try:
        with open('samples.txt', 'r') as f:
            samples = f.readlines()
        sample = samples[task_id - 1].strip()
        logging.info(f"Processing sample: {sample}")
    except Exception as e:
        logging.error(f"Error reading sample name: {str(e)}")
        sys.exit(1)

    input_dir = "emase_output"
    output_dir = "common_alignments_output"
    os.makedirs(output_dir, exist_ok=True)
    
    r1_file = os.path.join(input_dir, f"{sample}_R1.emase")
    r2_file = os.path.join(input_dir, f"{sample}_R2.emase")
    output_file = os.path.join(output_dir, f"{sample}_common.h5")

    if not all(os.path.exists(f) for f in [r1_file, r2_file]):
        logging.error(f"Input files missing for sample {sample}")
        sys.exit(1)
    
    try:
        stats = find_common_alignments(r1_file, r2_file, output_file)
        
        stats_file = os.path.join(output_dir, f"{sample}_statistics.txt")
        with open(stats_file, 'w') as f:
            f.write(f"Sample: {sample}\n")
            f.write(f"Total R1 reads: {stats['total_r1_reads']}\n")
            f.write(f"Total R2 reads: {stats['total_r2_reads']}\n")
            f.write(f"Common reads: {stats['common_reads']}\n")
        
        logging.info(f"Completed processing for sample: {sample}")
    except Exception as e:
        logging.error(f"Error processing sample {sample}: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
EOF