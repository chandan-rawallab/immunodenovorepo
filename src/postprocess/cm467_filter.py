import pandas as pd
from pathlib import Path
from pyteomics import fasta
import sys

input_tsv = "/home/amity/Documents/experiments/results/CM467_annotated_alleles.tsv"
output_tsv = "/home/amity/Documents/experiments/results/CM467_filtered_neoepitopes.tsv"
fasta_path = "/home/amity/Documents/experiments/data/reference/uniprot_human_reviewed.fasta"

def load_reference_peptides(fasta_path, lengths=[8, 9, 10, 11, 12, 13, 14, 15]):
    peptides = set()
    print(f"Indexing reference proteome from {fasta_path}...")
    for description, sequence in fasta.read(str(fasta_path)):
        for length in lengths:
            for i in range(len(sequence) - length + 1):
                peptides.add(sequence[i:i+length])
    print(f"Indexed {len(peptides)} unique reference peptides.")
    return peptides

def find_mutation(peptide, reference_peptides):
    aas = "ACDEFGHIKLMNPQRSTVWY"
    for i in range(len(peptide)):
        mut_aa = peptide[i]
        for wt_aa in aas:
            if wt_aa == mut_aa:
                continue
            potential_wt = peptide[:i] + wt_aa + peptide[i+1:]
            if potential_wt in reference_peptides:
                return potential_wt, i + 1, wt_aa, mut_aa
    return None

def main():
    print(f"Loading annotated data from {input_tsv}...")
    df = pd.read_csv(input_tsv, sep='\t')
    
    # Drop rows without a sequence
    df = df.dropna(subset=['Sequence'])
    
    ref_peptides = load_reference_peptides(fasta_path)
    
    results = []
    print("Performing missense mutation detection (Levenshtein distance 1)...")
    for _, row in df.iterrows():
        peptide = row['Sequence']
        
        if peptide in ref_peptides:
            # Self-peptide, not a missense neoepitope
            continue
            
        mutation_info = find_mutation(peptide, ref_peptides)
        if mutation_info:
            wt_seq, pos, wt_aa, mut_aa = mutation_info
            row_dict = row.to_dict()
            row_dict.update({
                'wildtype_peptide': wt_seq,
                'mutation_pos': pos,
                'wt_aa': wt_aa,
                'mut_aa': mut_aa,
                'mutation_type': 'missense'
            })
            results.append(row_dict)
            
    df_final = pd.DataFrame(results)
    if not df_final.empty:
        # Sort by binding rank
        df_final = df_final.sort_values(by='Binding_Rank')
        df_final.to_csv(output_tsv, sep="\t", index=False)
        print(f"Saved {len(df_final)} missense neoepitope candidates to {output_tsv}")
        print("\nTop 5 candidates by binding affinity:")
        print(df_final[['Sequence', 'wildtype_peptide', 'Best_Allele', 'Binding_Rank']].head())
    else:
        print("No missense neoepitope candidates found after filtering.")

if __name__ == "__main__":
    main()
