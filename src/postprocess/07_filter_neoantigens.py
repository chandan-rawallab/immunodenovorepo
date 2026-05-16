#!/usr/bin/env python3
"""Activity 3: FDR, length, and missense mutation filters for neoantigen candidates."""

import argparse
import pandas as pd
import re
from pathlib import Path
from pyteomics import fasta
from collections import defaultdict

# Canonical HLA-I length 8-11
CANONICAL_HLA_I_PATTERN = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]{8,11}$")

def load_reference_peptides(fasta_path, lengths=[8, 9, 10, 11]):
    """Index all possible peptides of specific lengths from the reference proteome."""
    peptides = set()
    print(f"Indexing reference proteome from {fasta_path}...")
    # Using pyteomics for fasta reading
    for description, sequence in fasta.read(str(fasta_path)):
        for length in lengths:
            for i in range(len(sequence) - length + 1):
                peptides.add(sequence[i:i+length])
    print(f"Indexed {len(peptides)} unique reference peptides.")
    return peptides

def find_mutation(peptide, reference_peptides):
    """
    Check if the peptide has exactly one mutation from a reference peptide.
    Returns (wt_peptide, pos, wt_aa, mut_aa) or None.
    """
    aas = "ACDEFGHIKLMNPQRSTVWY"
    for i in range(len(peptide)):
        mut_aa = peptide[i]
        for wt_aa in aas:
            if wt_aa == mut_aa:
                continue
            # Try substituting the current AA with wt_aa to see if it matches a reference peptide
            potential_wt = peptide[:i] + wt_aa + peptide[i+1:]
            if potential_wt in reference_peptides:
                return potential_wt, i + 1, wt_aa, mut_aa
    return None

def main():
    parser = argparse.ArgumentParser(description="Filter de novo candidates for neoantigens.")
    parser.add_argument("--input", type=Path, required=True, help="Path to de_novo_candidates.tsv")
    parser.add_argument("--psms", type=Path, required=True, help="Path to immunopeptidome_psms.tsv (database search hits)")
    parser.add_argument("--fasta", type=Path, required=True, help="Path to reference proteome FASTA")
    parser.add_argument("--output", type=Path, required=True, help="Output filtered TSV")
    parser.add_argument("--score_cutoff", type=float, default=0.7, help="De novo score cutoff")
    
    args = parser.parse_args()
    
    # 1. Load data
    print(f"Loading candidates from {args.input}...")
    df_candidates = pd.read_csv(args.input, sep="\t")
    
    print(f"Loading database search PSMs from {args.psms}...")
    df_psms = pd.read_csv(args.psms, sep="\t")
    # Identify peptides found in database search per sample to subtract them
    db_peptides = defaultdict(set)
    for _, row in df_psms.iterrows():
        db_peptides[str(row['sample_id'])].add(str(row['peptide']))
    
    # 2. Basic filters (Length and Score)
    initial_count = len(df_candidates)
    df = df_candidates[df_candidates['de_novo_score'] >= args.score_cutoff].copy()
    print(f"Filtered by score >= {args.score_cutoff}: {initial_count} -> {len(df)}")
    
    df = df[df['peptide'].str.match(CANONICAL_HLA_I_PATTERN)].copy()
    print(f"Filtered by length (8-11 AA): {len(df)}")
    
    # 3. Database subtraction (Self-peptides from this sample)
    def is_in_db(row):
        return row['peptide'] in db_peptides.get(str(row['sample_id']), set())
    
    mask_in_db = df.apply(is_in_db, axis=1)
    df = df[~mask_in_db].copy()
    print(f"Filtered out database hits: {len(df)}")
    
    # 4. Missense mutation detection
    ref_peptides = load_reference_peptides(args.fasta)
    
    results = []
    print("Performing missense mutation detection (Levenshtein distance 1)...")
    for _, row in df.iterrows():
        peptide = row['peptide']
        
        # Check if it's a perfect match in the whole proteome (Self-peptide from other proteins/samples)
        if peptide in ref_peptides:
            continue
            
        # Check for 1-off mutation
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
        else:
            # Optionally keep as Class C (non-missense, maybe frameshift or other)
            # For now, following the specific "single missense mutation" requirement
            pass
            
    df_final = pd.DataFrame(results)
    if not df_final.empty:
        df_final.to_csv(args.output, sep="\t", index=False)
        print(f"Saved {len(df_final)} neoantigen candidates to {args.output}")
    else:
        print("No neoantigen candidates found after filtering.")

if __name__ == "__main__":
    main()
