import os
import glob
import argparse
import pandas as pd
from tqdm import tqdm
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Go up two levels: src/data_prep/ -> src/ -> project_root/
SRC_DIR = os.path.dirname(SCRIPT_DIR)
WORKSPACE_ROOT = os.path.dirname(SRC_DIR)
sys.path.append(SRC_DIR)

try:
    from pyteomics import mgf
except ImportError:
    print("ERROR: pyteomics is required to write MGF files.")
    sys.exit(1)

def extract_unlabeled(mgf_dir, psm_file, output_dir, manifest_file=None, clean_output=False):
    """
    Reads all MGF files, filters out labelled spectra found in psm_file, 
    and writes the unlabeled spectra to new MGF files in output_dir.
    """
    # 1. Load labelled spectra
    labelled_spectra = set()
    if os.path.exists(psm_file):
        print(f"Loading PSM file: {psm_file}")
        psm_df = pd.read_csv(psm_file, sep='\t')
        for _, row in psm_df.iterrows():
            labelled_spectra.add((str(row['run_id']), str(row['spectrum_id'])))
        print(f"Loaded {len(labelled_spectra)} labelled spectra to exclude.")
    else:
        print(f"No PSM file found at {psm_file}. All spectra will be considered unlabelled.")

    os.makedirs(output_dir, exist_ok=True)
    
    # 2. Process MGF files
    mgf_files = glob.glob(os.path.join(mgf_dir, "*.mgf"))
    if not mgf_files:
        print(f"No MGF files found in {mgf_dir}")
        return

    # Filter by manifest if provided
    if manifest_file and os.path.exists(manifest_file):
        print(f"Filtering MGF files based on manifest: {manifest_file}")
        manifest_df = pd.read_csv(manifest_file, sep='\t')
        active_run_ids = set(manifest_df['run_id'].astype(str).tolist())
        mgf_files = [f for f in mgf_files if os.path.splitext(os.path.basename(f))[0] in active_run_ids]
        print(f"Filtered to {len(mgf_files)} MGF files matching manifest run_ids.")

        if clean_output:
            removed = 0
            for existing in glob.glob(os.path.join(output_dir, "*_unlabeled.mgf")):
                existing_run = os.path.basename(existing).replace("_unlabeled.mgf", "")
                if existing_run not in active_run_ids:
                    os.remove(existing)
                    removed += 1
            if removed:
                print(f"Removed {removed} stale unlabeled MGF files outside the active manifest.")

    if not mgf_files:
        print("No active MGF files left to process after manifest filtering.")
        return

    print(f"Found {len(mgf_files)} MGF files to process. Filtering...")

    total_kept = 0
    total_skipped = 0

    for mgf_path in tqdm(mgf_files, desc="Processing MGF files"):
        run_id = os.path.basename(mgf_path).replace(".mgf", "")
        out_mgf_path = os.path.join(output_dir, f"{run_id}_unlabeled.mgf")
        
        # We can stream reading and writing to avoid loading the whole MGF in memory
        try:
            reader = mgf.read(mgf_path)
            
            # Use generator to stream output
            def filtered_spectra_generator():
                nonlocal total_kept, total_skipped
                count = 0
                for spectrum in reader:
                    params = spectrum.get("params", {})
                    scan_id = str(params.get("scans", ""))
                    if not scan_id:
                        title = str(params.get("title", ""))
                        if "scan=" in title:
                            scan_id = title.split("scan=")[-1].strip()
                        else:
                            scan_id = str(count)
                            
                    if (run_id, scan_id) in labelled_spectra:
                        total_skipped += 1
                    else:
                        total_kept += 1
                        yield spectrum
                    count += 1

            # Write to output MGF
            mgf.write(filtered_spectra_generator(), out_mgf_path)
            
        except Exception as e:
            print(f"Error processing {mgf_path}: {e}")

    print(f"\nSUCCESS: Extracted unlabeled spectra.")
    print(f"Total labeled (skipped): {total_skipped}")
    print(f"Total unlabeled (kept): {total_kept}")
    print(f"Unlabeled MGF files saved to: {output_dir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract unlabeled spectra into separate MGF files.")
    parser.add_argument("--mgf-dir", type=str, 
                        default=os.path.join(WORKSPACE_ROOT, "data", "mgf"),
                        help="Directory containing original MGF files")
    parser.add_argument("--psms", type=str, 
                        default=os.path.join(WORKSPACE_ROOT, "results", "immunopeptidome_psms.tsv"),
                        help="Path to PSM file with labelled spectra")
    parser.add_argument("--output-dir", type=str, 
                        default=os.path.join(WORKSPACE_ROOT, "data", "mgf_unlabeled"),
                        help="Directory to save unlabeled MGF files")
    parser.add_argument("--manifest", type=str,
                        default=None,
                        help="Path to sample manifest file to filter run_ids")
    parser.add_argument("--clean-output", action="store_true",
                        help="Remove stale *_unlabeled.mgf files whose run_id is not in the manifest")
    
    args = parser.parse_args()
    extract_unlabeled(args.mgf_dir, args.psms, args.output_dir, args.manifest, args.clean_output)
