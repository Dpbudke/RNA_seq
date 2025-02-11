import os
import sys
import numpy as np
import h5py
from itertools import combinations_with_replacement
from optparse import OptionParser
import time

def expand_transition_matrix(tprob_c, haplotypes, geno_id, chr_num):
    """
    Expand transition matrix to 36x36 diploid matrix, verified version
    """
    start_time = time.time()
    num_genos = len(geno_id)
    num_positions = tprob_c.shape[0]
    
    print(f"\tStarting expansion for Chr {chr_num}")
    print(f"\tInput matrix shape: {tprob_c.shape}")
    print(f"\tOutput matrix shape will be: ({num_positions}, {num_genos}, {num_genos})")
    
    # Initialize with log of small probability
    min_prob = np.log(np.nextafter(0,1))
    tprob_c_big = np.full((num_positions, num_genos, num_genos), min_prob)
    
    # Check if this is X chromosome (5x5 matrix)
    is_x_chr = chr_num == 'X' and tprob_c.shape[1] == 5
    x_founder_map = {'A': 0, 'B': 1, 'C': 2, 'E': 3, 'F': 4} if is_x_chr else {}
    
    # Print initial homozygous transitions for verification
    print("\nChecking first position homozygous transitions:")
    if is_x_chr:
        for h1 in x_founder_map:
            idx = x_founder_map[h1]
            print(f"{h1}{h1}->{h1}{h1}: {tprob_c[0][idx, idx]:.10f}")
    else:
        for h1_idx, h1 in enumerate(haplotypes):
            print(f"{h1}{h1}->{h1}{h1}: {tprob_c[0][h1_idx, h1_idx]:.10f}")
    
    if is_x_chr:
        # Process X chromosome
        for pos in range(num_positions):
            if pos % (num_positions // 10) == 0:
                elapsed = time.time() - start_time
                print(f"\tProcessing position {pos}/{num_positions} "
                      f"({(pos/num_positions)*100:.1f}%) [Elapsed: {elapsed:.1f}s]")
            
            # Direct mapping of probabilities
            for h1 in x_founder_map:
                for h2 in x_founder_map:
                    idx1_5x5 = x_founder_map[h1]
                    idx2_5x5 = x_founder_map[h2]
                    prob = tprob_c[pos][idx1_5x5, idx2_5x5]
                    
                    geno1 = f"{h1}{h1}"
                    geno2 = f"{h2}{h2}"
                    
                    if geno1 in geno_id and geno2 in geno_id:
                        idx1_36x36 = geno_id[geno1]
                        idx2_36x36 = geno_id[geno2]
                        tprob_c_big[pos][idx1_36x36, idx2_36x36] = prob
    else:
        # Process autosomal chromosomes
        for pos in range(num_positions):
            if pos % (num_positions // 10) == 0:
                elapsed = time.time() - start_time
                print(f"\tProcessing position {pos}/{num_positions} "
                      f"({(pos/num_positions)*100:.1f}%) [Elapsed: {elapsed:.1f}s]")
            
            # Direct mapping of probabilities
            for h1_idx, h1 in enumerate(haplotypes):
                for h2_idx, h2 in enumerate(haplotypes):
                    geno1 = f"{h1}{h1}"
                    geno2 = f"{h2}{h2}"
                    
                    if geno1 in geno_id and geno2 in geno_id:
                        idx1 = geno_id[geno1]
                        idx2 = geno_id[geno2]
                        tprob_c_big[pos][idx1, idx2] = tprob_c[pos][h1_idx, h2_idx]
    
    # Print final matrix statistics
    print(f"\nMatrix statistics for Chr {chr_num}:")
    print(f"Value ranges:")
    print(f"Input  - min: {np.min(tprob_c):.10f}, max: {np.max(tprob_c):.10f}")
    print(f"Output - min: {np.min(tprob_c_big):.10f}, max: {np.max(tprob_c_big):.10f}")
    
    # Quick diagonal check on final matrix
    print("\nFinal diagonal check (first position):")
    if is_x_chr:
        for h1 in x_founder_map:
            geno = f"{h1}{h1}"
            if geno in geno_id:
                idx = geno_id[geno]
                print(f"{geno}->{geno}: {tprob_c_big[0][idx, idx]:.10f}")
    else:
        for h1 in haplotypes:
            geno = f"{h1}{h1}"
            if geno in geno_id:
                idx = geno_id[geno]
                print(f"{geno}->{geno}: {tprob_c_big[0][idx, idx]:.10f}")
    
    return tprob_c_big

def female_file_parse(options, haplotypes, num_haps, num_genos, geno_id):
    total_start = time.time()
    chrs = list(map(str, np.arange(19) + 1)) + ['X']
    print(f"\nOpening h5 file for female parsing...")
    
    try:
        h5fh = h5py.File(options.transition_h5, 'r')
    except Exception as e:
        print(f"Error opening h5 file: {e}")
        return
    
    gen = str(options.generation)
    tprob_out = {}
    
    print(f"Beginning processing for {len(chrs)} chromosomes...")
    for idx, c in enumerate(chrs, 1):
        chr_start = time.time()
        print(f"\nProcessing chromosome {c} ({idx}/{len(chrs)})")
        
        keystr = f"{c}:{gen}:F"
        try:
            tprob_c = h5fh[keystr]
            if c == 'X':
                print(f"\tX chromosome matrix shape: {tprob_c.shape}")
            tprob_out[c] = expand_transition_matrix(tprob_c, haplotypes, geno_id, c)
            
            chr_elapsed = time.time() - chr_start
            print(f"Completed chromosome {c} in {chr_elapsed:.1f} seconds")
        except Exception as e:
            print(f"Error processing chromosome {c}: {e}")
            continue

    output_file = f'tranprob.CC.G{gen}.F.npz'
    print(f"\nSaving compressed npz file: {output_file}")
    try:
        np.savez_compressed(output_file, **tprob_out)
        print("File saved successfully")
    except Exception as e:
        print(f"Error saving npz file: {e}")
    
    h5fh.close()
    total_elapsed = time.time() - total_start
    print(f"Completed all female processing in {total_elapsed:.1f} seconds")

def male_file_parse(options, haplotypes, num_haps, num_genos, geno_id):
    total_start = time.time()
    chrs = list(map(str, np.arange(19) + 1)) + ['X']
    print(f"\nOpening h5 file for male parsing...")
    
    try:
        h5fh = h5py.File(options.transition_h5, 'r')
    except Exception as e:
        print(f"Error opening h5 file: {e}")
        return
    
    gen = str(options.generation)
    tprob_out = {}
    
    print(f"Beginning processing for {len(chrs)} chromosomes...")
    for idx, c in enumerate(chrs, 1):
        chr_start = time.time()
        print(f"\nProcessing chromosome {c} ({idx}/{len(chrs)})")
        
        keystr = f"{c}:{gen}:M"
        try:
            tprob_c = h5fh[keystr]
            if c == 'X':
                print(f"\tX chromosome matrix shape: {tprob_c.shape}")
            tprob_out[c] = expand_transition_matrix(tprob_c, haplotypes, geno_id, c)
            
            chr_elapsed = time.time() - chr_start
            print(f"Completed chromosome {c} in {chr_elapsed:.1f} seconds")
        except Exception as e:
            print(f"Error processing chromosome {c}: {e}")
            continue

    output_file = f'tranprob.CC.G{gen}.M.npz'
    print(f"\nSaving compressed npz file: {output_file}")
    try:
        np.savez_compressed(output_file, **tprob_out)
        print("File saved successfully")
    except Exception as e:
        print(f"Error saving npz file: {e}")
    
    h5fh.close()
    total_elapsed = time.time() - total_start
    print(f"Completed all male processing in {total_elapsed:.1f} seconds")

if __name__ == "__main__":
    print("Starting transition probability parsing...")
    start_time = time.time()
    
    parser = OptionParser()
    parser.add_option("-t", "--transition_h5", dest="transition_h5",
                    help="transition prob h5 file from R", metavar="TRANSPROB_h5")
    parser.add_option("-l", "--haplotype_list", dest="haplotype_list", 
                    default="A, B, C, D, E, F, G, H",
                    help="haplotype list", metavar="HAP_LIST")
    parser.add_option("-s", "--sex", dest="sex", default="F",
                    help="sex (F or M)", metavar="SEX")
    parser.add_option("-g", "--generation", dest="generation", default="20",
                    help="number of generations", metavar="GEN_NUM")

    (options, args) = parser.parse_args()

    if not os.path.exists(options.transition_h5):
        sys.exit(f"Error: Input file {options.transition_h5} does not exist")

    print("\nParsing arguments...")
    haplotypes = [h.strip() for h in options.haplotype_list.split(",")]
    
    if not any([x in options.sex for x in ['F', 'M']]):
        sys.exit("'-s' or '--sex' must be either F or M, not: " + options.sex)

    num_haps = len(haplotypes)
    genotypes = [h1+h2 for h1, h2 in combinations_with_replacement(haplotypes, 2)]
    num_genos = len(genotypes)
    geno_id = dict(zip(genotypes, np.arange(num_genos)))

    print(f"Processing with {num_haps} haplotypes: {haplotypes}")
    print(f"Generated {num_genos} genotypes")
    print("Homozygous genotypes:", [g for g in genotypes if g[0] == g[1]])

    if 'F' in options.sex:
        female_file_parse(options, haplotypes, num_haps, num_genos, geno_id)
    if 'M' in options.sex:
        male_file_parse(options, haplotypes, num_haps, num_genos, geno_id)
        
    total_elapsed = time.time() - start_time
    print(f"\nTotal processing completed in {total_elapsed:.1f} seconds")
