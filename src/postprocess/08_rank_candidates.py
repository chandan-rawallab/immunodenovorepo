#!/usr/bin/env python3
"""Activity 3: Ranking neoantigen candidates by binding affinity and expression."""

import argparse
import pandas as pd
from pathlib import Path
import subprocess
import tempfile
import csv
import shutil

def normalize_hla(value: str) -> str:
    """Normalize HLA allele to MHCflurry format (e.g., HLA-A*02:01)."""
    value = value.strip()
    if not value:
        return ""
    
    # If it's already in a good format, return it
    if "*" in value and ":" in value:
        return value
    
    # Handle A0201 or HLA-A0201 formats
    clean = value.replace("HLA-", "").replace("*", "").replace(":", "")
    if len(clean) >= 5 and clean[0].isalpha():
        # A0201 -> HLA-A*02:01
        return f"HLA-{clean[0]}*{clean[1:3]}:{clean[3:5]}"
    
    return value

def run_mhcflurry(peptides, alleles, output_path):
    """Run MHCflurry prediction for a list of peptides and alleles."""
    predictor = shutil.which("mhcflurry-predict")
    if not predictor:
        print("Warning: mhcflurry-predict not found. Skipping binding prediction.")
        return False
    
    with tempfile.NamedTemporaryDirectory() as tmpdir:
        input_csv = Path(tmpdir) / "input.csv"
        with open(input_csv, "w") as f:
            f.write("allele,peptide\n")
            for allele in alleles:
                for peptide in peptides:
                    f.write(f"{allele},{peptide}\n")
        
        try:
            subprocess.run([predictor, str(input_csv), "--out", str(output_path)], check=True, capture_output=True)
            return True
        except subprocess.CalledProcessError as e:
            print(f"Error running MHCflurry: {e.stderr.decode()}")
            return False

def main():
    parser = argparse.ArgumentParser(description="Rank neoantigen candidates.")
    parser.add_argument("--input", type=Path, required=True, help="Path to filtered_neoantigens.tsv")
    parser.add_argument("--manifest", type=Path, required=True, help="Path to sample_manifest.tsv")
    parser.add_argument("--output", type=Path, required=True, help="Output ranked TSV")
    parser.add_argument("--binding_rank_cutoff", type=float, default=2.0)
    parser.add_argument("--tpm_cutoff", type=float, default=1.0)
    
    args = parser.parse_args()
    
    # 1. Load data
    df = pd.read_csv(args.input, sep="\t")
    manifest = pd.read_csv(args.manifest, sep="\t")
    
    # 2. Binding Prediction (MHCflurry)
    # Group by sample to get alleles
    results = []
    for sample_id, group in df.groupby("sample_id"):
        sample_meta = manifest[manifest['patient_id'] == sample_id] # Adjust if patient_id != sample_id
        if sample_meta.empty:
            sample_meta = manifest[manifest['run_id'] == sample_id]
            
        if sample_meta.empty:
            print(f"Warning: No metadata found for sample {sample_id}")
            continue
            
        # Get alleles from hla_alleles column or individual columns if they exist
        raw_alleles = str(sample_meta.iloc[0]['hla_alleles']).split(",")
        alleles = [normalize_hla(a) for a in raw_alleles if a and a.strip().lower() != "tbd"]
        
        if not alleles:
            print(f"Warning: No valid HLA alleles for sample {sample_id}")
            group_copy = group.copy()
            group_copy['best_hla'] = ""
            group_copy['binding_rank'] = None
            results.append(group_copy)
            continue
            
        peptides = group['peptide'].unique().tolist()
        with tempfile.NamedTemporaryDirectory() as tmpdir:
            out_csv = Path(tmpdir) / "mhcflurry_out.csv"
            if run_mhcflurry(peptides, alleles, out_csv):
                mhc_df = pd.read_csv(out_csv)
                # Pick best allele per peptide
                best_mhc = mhc_df.sort_values("presentation_percentile").groupby("peptide").first().reset_index()
                merged = group.merge(best_mhc[['peptide', 'allele', 'presentation_percentile']], on='peptide', how='left')
                merged.rename(columns={'allele': 'best_hla', 'presentation_percentile': 'binding_rank'}, inplace=True)
                results.append(merged)
            else:
                group_copy = group.copy()
                group_copy['best_hla'] = ""
                group_copy['binding_rank'] = None
                results.append(group_copy)
                
    if not results:
        print("No results to rank.")
        return
        
    df_ranked = pd.concat(results)
    
    # 3. Expression (Placeholder logic - linking by gene if available)
    # Assuming expression data might be added later or is in a separate file
    # For now, we'll set it to a default or check if 'expression_tpm' exists
    if 'expression_tpm' not in df_ranked.columns:
        df_ranked['expression_tpm'] = 0.0 # Placeholder
        
    # 4. Evidence Class
    def assign_class(row):
        has_binding = pd.notnull(row['binding_rank']) and row['binding_rank'] <= args.binding_rank_cutoff
        has_expression = row['expression_tpm'] >= args.tpm_cutoff
        
        if row.get('mutation_type') == 'missense' and has_binding and has_expression:
            return 'A'
        elif has_binding and has_expression:
            return 'B'
        else:
            return 'C'
            
    df_ranked['evidence_class'] = df_ranked.apply(assign_class, axis=1)
    
    # Sort by class and then by rank/score
    class_order = {'A': 0, 'B': 1, 'C': 2}
    df_ranked['class_sort'] = df_ranked['evidence_class'].map(class_order)
    df_ranked.sort_values(['class_sort', 'binding_rank'], ascending=[True, True], inplace=True)
    df_ranked.drop(columns=['class_sort'], inplace=True)
    
    df_ranked.to_csv(args.output, sep="\t", index=False)
    print(f"Saved {len(df_ranked)} ranked candidates to {args.output}")

if __name__ == "__main__":
    main()
