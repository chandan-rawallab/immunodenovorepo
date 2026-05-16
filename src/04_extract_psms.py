import os
import pandas as pd
import glob
from tqdm import tqdm
import argparse

def extract_psms(input_dir, output_file, manifest_file):
    """
    Extracts PSMs from MaxQuant msms.txt files and formats them for Objective 3 training.
    
    Args:
        input_dir: Directory containing MaxQuant msms.txt files (e.g. data/psms)
        output_file: Path to save the extracted PSMs (e.g. results/immunopeptidome_psms.tsv)
        manifest_file: Path to sample manifest to map filenames to patient_ids
    """
    if not os.path.exists(manifest_file):
        print(f"ERROR: Manifest file {manifest_file} not found.")
        return

    # Load manifest
    manifest = pd.read_csv(manifest_file, sep='\t')
    run_to_patient = dict(zip(manifest['run_id'], manifest['patient_id']))
    
    # Discover all msms.txt files
    psm_files = glob.glob(os.path.join(input_dir, "msms_*.txt"))
    if not psm_files:
        print(f"No PSM files found in {input_dir}")
        return

    print(f"Processing {len(psm_files)} PSM files...")
    
    all_psms = []
    
    for f in tqdm(psm_files, desc="Extracting PSMs"):
        # Pattern usually msms_<run_id>.raw.txt
        # Extract run_id from filename
        fname = os.path.basename(f)
        run_id = fname.replace("msms_", "").replace(".raw.txt", "").replace(".txt", "")
        patient_id = run_to_patient.get(run_id, "Unknown")
        
        try:
            df = pd.read_csv(f, sep='\t', low_memory=False)
            
            # Required columns
            required = ['Sequence', 'PEP', 'Scan number']
            if not all(col in df.columns for col in required):
                print(f"Skipping {fname}: Missing required columns {required}")
                continue
                
            # Filters
            # 1. Decoy != "+"
            if 'Decoy' in df.columns:
                df = df[df['Decoy'] != '+']
            
            # 2. PEP <= 0.01 (1% FDR proxy)
            df = df[df['PEP'] <= 0.01]
            
            # 3. Length 8-11 (Canonical HLA-I)
            df['Length'] = df['Sequence'].str.len()
            df = df[(df['Length'] >= 8) & (df['Length'] <= 11)]
            
            # 4. Remove empty sequences
            df = df.dropna(subset=['Sequence'])
            
            # Create output format
            subset = pd.DataFrame({
                'sample_id': patient_id,
                'run_id': run_id,
                'spectrum_id': df['Scan number'].astype(int).astype(str),
                'peptide': df['Sequence'],
                'pep_score': df['PEP']
            })
            
            all_psms.append(subset)
            
        except Exception as e:
            print(f"Error processing {fname}: {e}")

    if not all_psms:
        print("No PSMs extracted after filtering.")
        return

    final_df = pd.concat(all_psms, ignore_index=True)
    
    # Final cleanup: ensure spectrum_id is just the number (string)
    # The MGF reader might expect "scan:N" or "SCANS=N", we'll handle that in the Dataset class.
    
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    final_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nSUCCESS: Extracted {len(final_df)} PSMs.")
    print(f"Output saved to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract PSMs from MaxQuant results.")
    parser.add_argument("--input-dir", default="data/psms", help="Directory with MaxQuant msms.txt files")
    parser.add_argument("--output-file", default="results/immunopeptidome_psms.tsv", help="Output path")
    parser.add_argument("--manifest", default="configs/sample_manifest.tsv", help="Sample manifest path")
    
    args = parser.parse_args()
    
    extract_psms(args.input_dir, args.output_file, args.manifest)
