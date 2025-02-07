import pandas as pd
import os
from pathlib import Path

# Base directory and file paths
BASE_DIR = "/90daydata/do2_projects/CCPups_redo/EMASE"
GBRS_DATA = os.path.join(BASE_DIR, "GBRS_DATA")
FINAL_DIR = os.path.join(GBRS_DATA, "Haplotype**/GRCm39/transcript_mapped/liftover_output/transcript_mapped/Final/")
RECONSTRUCT_DIR = os.path.join(BASE_DIR, "reconstruct_output")

# Sample to strain mapping
VALIDATIONS = {
   "113C": "CC001", 
   "299C": "CC019",
   "154C": "CC040"
}

def load_gene2transcripts(filepath):
   """Load and parse the gene to transcripts mapping file."""
   gene2trans = {}
   with open(filepath, 'r') as f:
       for line in f:
           parts = line.strip().split('\t')
           gene2trans[parts[0]] = parts[1:]
   return gene2trans

def load_final_csv(strain):
   """Load the final CSV file for a given strain."""
   filepath = os.path.join(FINAL_DIR, f"{strain}_final.csv")
   return pd.read_csv(filepath)

def extract_diplotype(prob_state):
   """Extract diplotype code (e.g., 'GG' from 'GG=0.999')."""
   return prob_state.split('=')[0]

def get_prob_state(final_df, transcript_id):
   """Get Primary_Prob_State diplotype code for a transcript."""
   row = final_df[final_df['transcript_id'] == transcript_id]
   if not row.empty:
       prob_state = row['Primary_Prob_State'].iloc[0]
       return extract_diplotype(prob_state)
   return None

def validate_sample(sample_id, gene2trans, strain_final_df):
   """Validate genotypes for a sample."""
   genotypes_file = os.path.join(RECONSTRUCT_DIR, sample_id, f"{sample_id}.genotypes.tsv")
   genotypes_df = pd.read_csv(genotypes_file, sep='\t')
   
   mismatch_count = 0
   for _, row in genotypes_df.iterrows():
       gene_id = row['#Gene_ID']
       diplotype = row['Diplotype']
       
       if gene_id in gene2trans:
           transcripts = gene2trans[gene_id]
           for transcript in transcripts:
               prob_state = get_prob_state(strain_final_df, transcript)
               if prob_state and prob_state != diplotype:
                   mismatch_count += 1
                   break  # Count one mismatch per gene, even if multiple transcripts mismatch
   
   return mismatch_count, len(genotypes_df)

def main():
   # Load gene to transcript mappings
   gene2trans = load_gene2transcripts(os.path.join(GBRS_DATA, "emase.gene2transcripts.tsv"))
   
   # Process each validation sample
   for sample_id, strain in VALIDATIONS.items():
       print(f"\nValidating {sample_id} against {strain}")
       
       try:
           # Load strain final CSV
           strain_final_df = load_final_csv(strain)
           
           # Validate
           mismatch_count, total_genes = validate_sample(sample_id, gene2trans, strain_final_df)
           
           # Report results
           print(f"Results: {mismatch_count}/{total_genes} genes mismatched")
           print(f"Match rate: {((total_genes - mismatch_count)/total_genes)*100:.2f}%")
       
       except Exception as e:
           print(f"Error processing {sample_id}: {str(e)}")

if __name__ == "__main__":
   main()