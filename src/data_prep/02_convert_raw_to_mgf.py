import os
import subprocess
import glob

# =========================================================================
# DYNAMIC PATH RESOLUTION
# =========================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))

# Path definitions relative to Project Root
RAW_DIR = os.path.join(PROJECT_ROOT, "data", "raw")
MGF_OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "mgf")
CONVERTER = os.path.join(PROJECT_ROOT, "bin", "ThermoRawFileParser", "ThermoRawFileParser")

def run_conversion():
    # 1. Create output directory
    os.makedirs(MGF_OUTPUT_DIR, exist_ok=True)
    
    # 2. Find all .raw files recursively in the localized data directory
    print(f"Scanning for RAW files in: {RAW_DIR}")
    raw_files = glob.glob(os.path.join(RAW_DIR, "**", "*.raw"), recursive=True)
    
    if not raw_files:
        print(f"[ERROR] No .raw files found in {RAW_DIR}")
        return

    print(f"Discovered {len(raw_files)} RAW files for conversion.")
    
    # 3. Verify Requirements
    DOTNET = os.path.join(PROJECT_ROOT, "bin", "dotnet8", "dotnet")
    CONVERTER_DLL = os.path.join(PROJECT_ROOT, "bin", "ThermoRawFileParser", "ThermoRawFileParser.dll")
    
    if not os.path.exists(DOTNET):
        print(f"[ERROR] dotnet binary not found at {DOTNET}")
        return
    if not os.path.exists(CONVERTER_DLL):
        print(f"[ERROR] Converter DLL not found at {CONVERTER_DLL}")
        return

    # 4. Conversion Process
    RESULTS_DIR = os.path.join(PROJECT_ROOT, "data", "psms")
    
    for i, raw_path in enumerate(raw_files):
        filename = os.path.basename(raw_path)
        # Check if result exists: msms_{filename}.txt
        # filename is e.g. sample.raw, result is msms_sample.raw.txt
        result_filename = f"msms_{filename}.txt"
        result_path = os.path.join(RESULTS_DIR, result_filename)
        
        if not os.path.exists(result_path):
            print(f"--- Skipping [{i+1}/{len(raw_files)}]: {filename} (No search results yet) ---")
            continue

        mgf_filename = filename.replace(".raw", ".mgf")
        mgf_path = os.path.join(MGF_OUTPUT_DIR, mgf_filename)
        
        if os.path.exists(mgf_path):
            print(f"--- Skipping [{i+1}/{len(raw_files)}]: {filename} (Already converted) ---")
            continue
            
        print(f"\n--- Converting [{i+1}/{len(raw_files)}]: {filename} ---")
        
        cmd = [
            DOTNET,
            "exec",
            CONVERTER_DLL,
            "-i", raw_path,
            "-o", MGF_OUTPUT_DIR,
            "-f", "0"
        ]
        
        try:
            subprocess.run(cmd, check=True)
            print(f"SUCCESS: Converted {filename} to MGF")
        except subprocess.CalledProcessError as e:
            print(f"FAILED: Error converting {filename}")
        except Exception as e:
            print(f"ERROR: Unexpected error: {e}")

    print(f"\n[COMPLETE] All converted files are in: {MGF_OUTPUT_DIR}")

if __name__ == "__main__":
    run_conversion()
