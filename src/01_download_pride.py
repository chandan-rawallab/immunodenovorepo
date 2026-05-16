import sys
import os
import requests
import subprocess
from pathlib import Path

# =========================================================================
# UNIVERSAL PROTEOMICS DATA RETRIEVAL TOOL (PRIDE & MASSIVE)
# =========================================================================

def get_pride_urls(project_id):
    """Fetches high-speed FTP URLs from the EBI PRIDE API."""
    url = f"https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{project_id}/files"
    print(f"[PRIDE] Fetching metadata for {project_id}...")
    
    try:
        resp = requests.get(url, headers={"Accept": "application/json"}, timeout=30)
        resp.raise_for_status()
        data = resp.json()
    except Exception as e:
        print(f"[ERROR] PRIDE API Error: {e}")
        return []

    urls = []
    for item in data:
        if item.get("fileName", "").lower().endswith(".raw"):
            for loc in item.get("publicFileLocations", []):
                if "FTP" in loc.get("name", ""):
                    urls.append(loc.get("value", ""))
                    break
    return urls

def download_massive(project_id, output_dir):
    """Uses recursive FTP scraping to retrieve raw data from UCSD MassIVE."""
    # Standard MassIVE FTP Pattern
    ftp_root = f"ftp://massive.ucsd.edu/{project_id}/"
    raw_folder = f"{ftp_root}raw/"
    
    print(f"[MASSIVE] Initializing recursive scan for: {raw_folder}")
    
    # We use wget --spider to first check if the /raw/ folder exists, 
    # fallback to root if /raw/ is missing (some older projects vary).
    try:
        download_dir = os.path.join(output_dir, project_id)
        os.makedirs(download_dir, exist_ok=True)
        
        # Wget command for MassIVE:
        # -r: recursive
        # -l 10: depth limit
        # -nd: no directories (keep files flat in the output)
        # -A: accept only .raw files
        # -np: dont go to parent
        cmd = [
            "wget", "-r", "-np", "-nd", "-l", "10",
            "-A", ".raw,.RAW",
            "-P", download_dir,
            "-c", # continue
            "--no-passive-ftp",
            raw_folder
        ]
        
        print(f"[MASSIVE] Launching FTP scraper...")
        subprocess.run(cmd, check=True)
        print(f"[SUCCESS] MassIVE data retrieved to {download_dir}")
        
    except Exception as e:
        print(f"[ERROR] MassIVE download failed: {e}")
        print("Tip: If 'raw/' folder failed, trying root folder...")
        # Fallback to root scan if /raw/ specifically doesn't exist
        cmd[-1] = ftp_root
        subprocess.run(cmd, check=False)

def download_pride(project_id, output_dir):
    """Standard PRIDE downloader using the gathered URLs."""
    urls = get_pride_urls(project_id)
    if not urls:
        print(f"[ERROR] No .raw files found for {project_id}")
        return

    download_dir = os.path.join(output_dir, project_id)
    os.makedirs(download_dir, exist_ok=True)
    
    # Save URLs to temp file for wget -i
    temp_list = "/tmp/pride_list.txt"
    with open(temp_list, "w") as f:
        for u in urls:
            f.write(f"{u}\n")
            
    print(f"[PRIDE] Downloading {len(urls)} files to {download_dir}...")
    cmd = [
        "wget", "-i", temp_list, 
        "-c", "-P", download_dir, 
        "-q", "--show-progress"
    ]
    subprocess.run(cmd, check=True)
    print(f"[SUCCESS] PRIDE data retrieved.")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Universal Proteomics Downloader (PRIDE/MassIVE)")
    parser.add_argument("-i", "--id", required=True, help="Dataset ID (e.g., PXD005231 or MSV000080620)")
    parser.add_argument("-o", "--output", default="data/raw", help="Root storage directory")
    
    args = parser.parse_args()
    project_id = args.id.upper()
    
    # Route based on ID prefix
    if project_id.startswith("MSV"):
        download_massive(project_id, args.output)
    elif project_id.startswith("PXD"):
        download_pride(project_id, args.output)
    else:
        print("[ERROR] Unrecognized ID format. Must start with PXD or MSV.")
