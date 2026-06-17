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
    """Index all 8-11mer peptides from human proteins, or all entries if the FASTA has no organism tags."""
    records = list(fasta.read(str(fasta_path)))
    has_organism_tags = any("OS=" in desc or "OX=" in desc for desc, _ in records)
    peptides = set()
    print(f"Indexing reference proteome from {fasta_path}...")
    for desc, sequence in records:
        if has_organism_tags and "OS=Homo sapiens" not in desc and "OX=9606" not in desc:
            continue
        for length in lengths:
            seq_len = len(sequence)
            for i in range(seq_len - length + 1):
                peptides.add(sequence[i:i+length])
    print(f"Indexed {len(peptides)} unique reference peptides.")
    return peptides

def find_mutation(peptide, reference_set):
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
            potential_wt = peptide[:i] + wt_aa + peptide[i+1:]
            if potential_wt in reference_set:
                return potential_wt, i + 1, wt_aa, mut_aa
    return None

def find_source_proteins(peptide, fasta_path):
    """Find Uniprot IDs of proteins containing the exact peptide."""
    records = list(fasta.read(str(fasta_path)))
    has_organism_tags = any("OS=" in desc or "OX=" in desc for desc, _ in records)
    sources = []
    for description, sequence in records:
        if has_organism_tags and "OS=Homo sapiens" not in description and "OX=9606" not in description:
            continue
        if peptide in sequence:
            parts = description.split('|')
            protein_id = parts[1] if len(parts) > 1 else description.split()[0]
            sources.append(protein_id)
    return sources

def main():
    parser = argparse.ArgumentParser(description="Filter de novo candidates for neoantigens.")
    parser.add_argument("--input", type=Path, required=True, help="Path to de_novo_candidates.tsv")
    parser.add_argument("--psms", type=Path, required=True, help="Path to immunopeptidome_psms.tsv (database search hits)")
    parser.add_argument("--fasta", type=Path, required=True, help="Path to reference proteome FASTA")
    parser.add_argument("--manifest", type=Path, default="configs/sample_manifest.tsv", help="Sample manifest for run_id to sample_id mapping")
    parser.add_argument("--output", type=Path, required=True, help="Output filtered TSV")
    parser.add_argument("--score_cutoff", type=float, default=-0.5, help="De novo score cutoff (default: -0.5 for log-probability scores)")
    parser.add_argument("--min_psm_support", type=int, default=2,
                        help="Minimum number of spectra supporting a peptide (proposal: >= 2)")
    parser.add_argument("--allow_flanking_mutations", action="store_true",
                        help="If set, retain missense mutations at flanking positions (pos 1 or last). Default: exclude.")
    
    args = parser.parse_args()
    
    # 1. Load data
    print(f"Loading candidates from {args.input}...")
    df_candidates = pd.read_csv(args.input, sep="\t")
    
    # Map run_id to sample_id via manifest
    manifest = pd.read_csv(args.manifest, sep="\t")
    run_to_sample = dict(zip(manifest['run_id'], manifest['patient_id']))
    df_candidates['sample_id'] = df_candidates['run_id'].map(run_to_sample).fillna("Unknown")
    
    if 'score' not in df_candidates.columns and 'de_novo_score' in df_candidates.columns:
        df_candidates = df_candidates.rename(columns={'de_novo_score': 'score'})
    if 'score' not in df_candidates.columns:
        raise ValueError("Candidate table must contain `score` or `de_novo_score`.")

    print(f"Loading database search PSMs from {args.psms}...")
    df_psms = pd.read_csv(args.psms, sep="\t")
    # Identify peptides found in database search per sample to subtract them
    db_peptides = defaultdict(set)
    for _, row in df_psms.iterrows():
        db_peptides[str(row['sample_id'])].add(str(row['peptide']))
    
    # 2. Basic filters (Length and Score)
    initial_count = len(df_candidates)
    df = df_candidates[df_candidates['score'] >= args.score_cutoff].copy()
    print(f"Filtered by score >= {args.score_cutoff}: {initial_count} -> {len(df)}")
    
    df = df[df['peptide'].str.match(CANONICAL_HLA_I_PATTERN)].copy()
    print(f"Filtered by length (8-11 AA): {len(df)}")
    
    # 3. Database subtraction (Self-peptides from this sample)
    def is_in_db(row):
        return row['peptide'] in db_peptides.get(str(row['sample_id']), set())
    
    mask_in_db = df.apply(is_in_db, axis=1)
    df = df[~mask_in_db].copy()
    print(f"Filtered out database hits: {len(df)}")
    
    # 3b. PSM support count >= min_psm_support (proposal requirement)
    # A peptide must appear in at least N independent spectra to be retained.
    if 'run_id' in df.columns:
        psm_counts = df.groupby(['sample_id', 'peptide'])['run_id'].count().rename('psm_count')
    else:
        psm_counts = df.groupby(['sample_id', 'peptide']).size().rename('psm_count')
    df = df.join(psm_counts, on=['sample_id', 'peptide'])
    before = len(df)
    df = df[df['psm_count'] >= args.min_psm_support].copy()
    print(f"Filtered by PSM support >= {args.min_psm_support}: {before} -> {len(df)}")
    
    # 4. Missense mutation detection
    print("Indexing reference proteome (this may take a moment)...")
    ref_set = load_reference_peptides(args.fasta)

    results = []
    print("Performing missense mutation detection (Levenshtein distance 1)...")
    for _, row in df.iterrows():
        peptide = row['peptide']
        
        # Check if it's a perfect match in the whole proteome (Self-peptide from other proteins/samples)
        if peptide in ref_set:
            continue

        # Check for 1-off mutation
        mutation_info = find_mutation(peptide, ref_set)
        row_dict = row.to_dict()

        if mutation_info:
            wt_seq, pos, wt_aa, mut_aa = mutation_info

            # Proposal filter: exclude mutations at flanking positions (pos 1 or pos len)
            # Flanking AA changes often affect proteasomal cleavage, not TCR contact
            is_flanking = (pos == 1 or pos == len(peptide))
            if is_flanking and not args.allow_flanking_mutations:
                continue

            # Find source proteins (lazy — only for the few matched candidates)
            source_proteins = find_source_proteins(wt_seq, args.fasta)

            row_dict.update({
                'wildtype_peptide': wt_seq,
                'mutation_pos': pos,
                'wt_aa': wt_aa,
                'mut_aa': mut_aa,
                'mutation_type': 'missense',
                'source_protein': ";".join(source_proteins)
            })
            results.append(row_dict)
        else:
            # Keep as Class C (non-missense, maybe frameshift or other)
            row_dict.update({
                'wildtype_peptide': '',
                'mutation_pos': '',
                'wt_aa': '',
                'mut_aa': '',
                'mutation_type': 'other',
                'source_protein': ''
            })
            results.append(row_dict)
            
    df_final = pd.DataFrame(results)
    n_missense = (df_final['mutation_type'] == 'missense').sum() if not df_final.empty else 0
    n_other    = (df_final['mutation_type'] == 'other').sum() if not df_final.empty else 0
    
    print(f"\n=== Filter Summary ===")
    print(f"  Total candidates retained: {len(df_final)}")
    print(f"  Missense (Class A/B eligible): {n_missense}")
    print(f"  Other / Class C:               {n_other}")
    
    df_final.to_csv(args.output, sep="\t", index=False)
    print(f"\nSaved {len(df_final)} neoantigen candidates to {args.output}")


if __name__ == "__main__":
    main()
