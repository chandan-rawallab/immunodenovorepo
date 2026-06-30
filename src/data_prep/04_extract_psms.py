import os
import pandas as pd
import glob
from tqdm import tqdm
import argparse

def extract_psms(input_dir, output_file, manifest_file, force=False):
    """
    Extracts PSMs from MaxQuant msms.txt files and formats them for Objective 3 training.
    
    Args:
        input_dir: Directory containing MaxQuant msms.txt files (e.g. data/psms)
        output_file: Path to save the extracted PSMs (e.g. results/immunopeptidome_psms.tsv)
        manifest_file: Path to sample manifest to map filenames to patient_ids
        force: Overwrite output file if it exists
    """
    if os.path.exists(output_file) and not force:
        print(f"Output file {output_file} already exists. Skipping extraction. Use --force to overwrite.")
        return

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

    manifest_runs = set(str(x) for x in manifest['run_id'].dropna())
    strict_manifest = os.environ.get("ALLOW_UNKNOWN_RUNS", "").lower() not in {"1", "true", "yes"}

    print(f"Processing {len(psm_files)} PSM files...")
    if strict_manifest:
        print("Strict manifest mode enabled: runs missing from manifest will be skipped.")
    
    all_psms = []
    skipped_unmapped = []
    
    for f in tqdm(psm_files, desc="Extracting PSMs"):
        # Pattern usually msms_<run_id>.raw.txt
        # Extract run_id from filename
        fname = os.path.basename(f)
        run_id = fname.replace("msms_", "").replace(".raw.txt", "").replace(".txt", "")
        if strict_manifest and run_id not in manifest_runs:
            skipped_unmapped.append(run_id)
            continue
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
        if skipped_unmapped:
            print(f"Skipped {len(skipped_unmapped)} PSM files not present in manifest.")
        print("No PSMs extracted after filtering.")
        return

    final_df = pd.concat(all_psms, ignore_index=True)
    
    # Final cleanup: ensure spectrum_id is just the number (string)
    # The MGF reader might expect "scan:N" or "SCANS=N", we'll handle that in the Dataset class.
    
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    final_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nSUCCESS: Extracted {len(final_df)} PSMs.")
    if skipped_unmapped:
        print(f"Skipped {len(skipped_unmapped)} PSM files not present in manifest.")
        print("First skipped run IDs:", ", ".join(skipped_unmapped[:10]))
    print(f"Output saved to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract PSMs from MaxQuant results.")
    parser.add_argument("--input-dir", default="data/psms", help="Directory with MaxQuant msms.txt files")
    parser.add_argument("--output-file", default="results/immunopeptidome_psms.tsv", help="Output path")
    parser.add_argument("--manifest", default="configs/sample_manifest.tsv", help="Sample manifest path")
    parser.add_argument("--force", action="store_true", help="Overwrite output file if it exists")
    
    args = parser.parse_args()
    
    extract_psms(args.input_dir, args.output_file, args.manifest, args.force)
