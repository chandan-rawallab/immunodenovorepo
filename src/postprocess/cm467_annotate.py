import pandas as pd
import os
import subprocess
import tempfile
from pathlib import Path

# Define paths
extracted_data_path = "/home/amity/Documents/experiments/results/CM467_extracted_data.tsv"
annotated_output_path = "/home/amity/Documents/experiments/results/CM467_annotated_alleles.tsv"
mhcflurry_bin = "/home/amity/Documents/experiments/.venv/bin/mhcflurry-predict"

# Alleles for CM467 from s1_supporting_info.txt
raw_alleles = ["A01:01", "A24:02", "B13:02", "B39:06", "C06:02", "C12:03"]

def normalize_hla(value):
    value = value.strip()
    if "*" in value and ":" in value:
        return value
    if len(value) >= 5:
        # A01:01 -> HLA-A*01:01
        return f"HLA-{value[0]}*{value[1:]}"
    return value

alleles = [normalize_hla(a) for a in raw_alleles]
print(f"Normalized alleles: {alleles}")

def main():
    if not os.path.exists(extracted_data_path):
        print(f"Error: {extracted_data_path} not found.")
        return

    print(f"Reading extracted data from {extracted_data_path}...")
    df = pd.read_csv(extracted_data_path, sep='\t')
    
    # Filter for standard amino acids only (MHCflurry requirement)
    # Peptides must be 8-15 AA and only standard AA
    peptides = df['Sequence'].unique().tolist()
    valid_peptides = [p for p in peptides if len(p) >= 8 and len(p) <= 15 and all(c in "ACDEFGHIKLMNPQRSTVWY" for c in p)]
    
    print(f"Running MHCflurry for {len(valid_peptides)} unique peptides across {len(alleles)} alleles...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        input_csv = Path(tmpdir) / "input.csv"
        output_csv = Path(tmpdir) / "output.csv"
        
        # Prepare input for MHCflurry
        input_data = []
        for allele in alleles:
            for peptide in valid_peptides:
                input_data.append({"allele": allele, "peptide": peptide})
        
        pd.DataFrame(input_data).to_csv(input_csv, index=False)
        
        # Run MHCflurry
        try:
            subprocess.run([mhcflurry_bin, str(input_csv), "--out", str(output_csv)], check=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running MHCflurry: {e.stderr.decode()}")
            return
        
        # Read results
        mhc_df = pd.read_csv(output_csv)
        
        # Pick best allele per peptide (lowest presentation_percentile or affinity_percentile)
        # Using mhcflurry_presentation_percentile if available, else mhcflurry_affinity_percentile
        score_col = 'mhcflurry_presentation_percentile' if 'mhcflurry_presentation_percentile' in mhc_df.columns else 'mhcflurry_affinity_percentile'
        print(f"Using {score_col} for ranking.")
        
        best_mhc = mhc_df.sort_values(score_col).groupby("peptide").first().reset_index()
        
        # Merge back to original data
        result_df = df.merge(best_mhc[['peptide', 'allele', score_col]], left_on='Sequence', right_on='peptide', how='left')
        result_df.rename(columns={'allele': 'Best_Allele', score_col: 'Binding_Rank'}, inplace=True)
        result_df.drop(columns=['peptide'], inplace=True)
        
        # Save results
        result_df.to_csv(annotated_output_path, sep='\t', index=False)
        print(f"Saved annotated data to {annotated_output_path}")
        
        # Summary stats
        print("\nAnnotation Summary:")
        print(f"Peptides successfully annotated: {result_df['Best_Allele'].notna().sum()} / {len(result_df)}")
        print("\nAllele distribution:")
        print(result_df['Best_Allele'].value_counts())

if __name__ == "__main__":
    main()
