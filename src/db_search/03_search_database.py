import os
import subprocess
import shutil
import xml.etree.ElementTree as ET

# =========================================================================
# CONFIGURATION (RELATIVE PATHS FOR PORTABILITY)
# =========================================================================
# Get the directory where THIS script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))

# Path definitions relative to Project Root
RAW_DIR = os.path.join(PROJECT_ROOT, "data", "raw")
RESULT_DIR = os.path.join(PROJECT_ROOT, "data", "psms")
BIN_DIR = os.path.join(PROJECT_ROOT, "bin")

# Execution Requirements
MQ_BIN = os.path.join(BIN_DIR, "dotnet8", "dotnet")
MQ_DLL = os.path.join(BIN_DIR, "MaxQuant_v2.8.0.0", "bin", "MaxQuantCmd.dll")
TEMPLATE_XML = os.path.join(PROJECT_ROOT, "configs", "mqpar_production.xml")
FASTA_PATH = os.path.join(PROJECT_ROOT, "data", "reference", "uniprot_human_reviewed.fasta")

# Ensure directories exist
os.makedirs(RESULT_DIR, exist_ok=True)

def run_sequential_search():
    # 1. Get all raw files recursively from data/raw
    import glob
    raw_files = glob.glob(os.path.join(RAW_DIR, "**", "*.raw"), recursive=True)
    raw_files.sort()
    
    if not raw_files:
        print(f"[ERROR] No .raw files found in {RAW_DIR}")
        return

    print(f"Found {len(raw_files)} raw files in pool.")

    for i, full_path in enumerate(raw_files):
        raw_filename = os.path.basename(full_path)
        dest_msms = os.path.join(RESULT_DIR, f"msms_{raw_filename}.txt")
        
        # SKIP LOGIC: Check if we already have the result
        if os.path.exists(dest_msms):
            print(f"--- Skipping [{i+1}/{len(raw_files)}]: {raw_filename} (Result already exists) ---")
            continue

        print(f"\n--- Processing [{i+1}/{len(raw_files)}]: {raw_filename} ---", flush=True)
        
        # 2. Patch XML for this specific file
        try:
            tree = ET.parse(TEMPLATE_XML)
            root = tree.getroot()
        except Exception as e:
            print(f"[ERROR] Could not parse template XML: {e}")
            return
        
        # Find and update filePaths
        file_paths = root.find("filePaths")
        # Clear existing ones
        for s in list(file_paths):
            file_paths.remove(s)
        # Add current file (Use realpath to bypass symlink metadata issues)
        new_string = ET.SubElement(file_paths, "string")
        new_string.text = os.path.realpath(full_path)

        # Update Fasta file path (MUST BE ABSOLUTE)
        fasta_info = root.find(".//FastaFileInfo")
        if fasta_info is not None:
            fasta_node = fasta_info.find("fastaFilePath")
            if fasta_node is not None:
                fasta_node.text = FASTA_PATH
        
        # Find and update experiment name
        experiments = root.find("experiments")
        for s in list(experiments):
            experiments.remove(s)
        new_exp = ET.SubElement(experiments, "string")
        new_exp.text = raw_filename.replace(".raw", "")
        
        # Save temp XML in temp directory (auto-cleaned after each run)
        temp_xml = os.path.join(PROJECT_ROOT, "temp", f"mqpar_temp_{i}.xml")
        os.makedirs(os.path.dirname(temp_xml), exist_ok=True)
        tree.write(temp_xml)
        
        # 3. Execute MaxQuant
        print(f"Starting MaxQuant search for {raw_filename}...", flush=True)
        
        try:
            # We run sequentially to prevent OOM crashes on the server RAM
            # Automated Timeout: 8 hours (28800 seconds). Prevents infinite hangs on corrupted files.
            subprocess.run([MQ_BIN, MQ_DLL, temp_xml], check=True, timeout=28800)
            print(f"Search successful for {raw_filename}")
            
            # 4. Extract results
            # IMPORTANT: Use realpath to find the combined folder next to the actual raw file
            raw_real_path = os.path.realpath(full_path)
            raw_parent = os.path.dirname(raw_real_path)
            search_results_dir = os.path.join(raw_parent, "combined", "txt")
            source_msms = os.path.join(search_results_dir, "msms.txt")
            
            if os.path.exists(source_msms):
                dest_msms = os.path.join(RESULT_DIR, f"msms_{raw_filename}.txt")
                shutil.copy(source_msms, dest_msms)
                print(f"Saved results to {dest_msms}", flush=True)
            else:
                print(f"WARNING: msms.txt not found for {raw_filename}")
                
            # 5. Cleanup combined folder to save disk space for the next run
            shutil.rmtree(os.path.join(raw_parent, "combined"), ignore_errors=True)
            if os.path.exists(temp_xml):
                os.remove(temp_xml)
            
        except subprocess.TimeoutExpired as e:
            print(f"ERROR: MaxQuant timed out (hung) for {raw_filename} after 8 hours. Skipping.")
            # Kill orphaned MaxQuant child processes just in case
            raw_parent = os.path.dirname(os.path.realpath(full_path))
            shutil.rmtree(os.path.join(raw_parent, "combined"), ignore_errors=True)
            continue
        except subprocess.CalledProcessError as e:
            raw_parent = os.path.dirname(os.path.realpath(full_path))
            shutil.rmtree(os.path.join(raw_parent, "combined"), ignore_errors=True)
            continue
        except Exception as e:
            raw_parent = os.path.dirname(os.path.realpath(full_path))
            shutil.rmtree(os.path.join(raw_parent, "combined"), ignore_errors=True)
            continue

if __name__ == "__main__":
    run_sequential_search()
