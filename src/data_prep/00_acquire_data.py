#!/usr/bin/env python3
"""
Task 2: Smart Data Acquisition Script
Acquires proteomics raw files and search results from PRIDE, MassIVE, or EGA.
"""

import sys
import os
import requests
import subprocess
import shutil
import json
import argparse
import zipfile
from pathlib import Path

def download_file(url, dest_path, resume=True):
    """Downloads a file using requests with progress logging, supporting resume/range requests."""
    dest_path = Path(dest_path)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Try converting ftp:// to https:// for PRIDE EBI URLs since HTTP is often faster and supports Range
    if url.startswith("ftp://ftp.pride.ebi.ac.uk/"):
        url = url.replace("ftp://ftp.pride.ebi.ac.uk/", "https://ftp.pride.ebi.ac.uk/")

    headers = {}
    existing_size = 0
    if resume and dest_path.exists():
        existing_size = dest_path.stat().st_size
        headers['Range'] = f"bytes={existing_size}-"

    print(f"Downloading {url} -> {dest_path}")
    try:
        # Stream response
        r = requests.get(url, headers=headers, stream=True, timeout=60)
        
        # If we requested Range and server returned 206, append to the file
        if resume and r.status_code == 206:
            mode = "ab"
            print(f"Resuming download from byte {existing_size}...")
        else:
            mode = "wb"
            if r.status_code != 200:
                # If server doesn't support Range or failed, fall back
                r = requests.get(url, stream=True, timeout=60)
                r.raise_for_status()

        total_size = int(r.headers.get('content-length', 0)) + existing_size
        block_size = 1024 * 1024  # 1MB blocks
        downloaded = existing_size

        with open(dest_path, mode) as f:
            for chunk in r.iter_content(chunk_size=block_size):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total_size > 0:
                        percent = (downloaded / total_size) * 100
                        sys.stdout.write(f"\rProgress: {percent:.2f}% ({downloaded}/{total_size} bytes)")
                        sys.stdout.flush()
        print("\nDownload finished.")
        return True
    except Exception as e:
        print(f"\n[WARNING] Requests download failed: {e}. Trying wget fallback...")
        try:
            # Fallback to wget
            cmd = ["wget", "-c", "-O", str(dest_path), url]
            subprocess.run(cmd, check=True)
            return True
        except Exception as wget_e:
            print(f"[ERROR] Both download methods failed. Wget error: {wget_e}")
            return False

def extract_zip(zip_path, extract_to):
    """Extracts a zip file to the target directory."""
    print(f"Extracting {zip_path} -> {extract_to}")
    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(extract_to)
        print("Extraction complete.")
        return True
    except Exception as e:
        print(f"[ERROR] Failed to extract zip: {e}")
        return False

def process_pride(accession, raw_dir, psm_dir, skip_raw):
    """Acquires data from PRIDE."""
    # 1. Fetch project information
    project_url = f"https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{accession}"
    print(f"[PRIDE] Fetching metadata for {accession}...")
    try:
        resp = requests.get(project_url, headers={"Accept": "application/json"}, timeout=30)
        resp.raise_for_status()
        metadata = resp.json()
    except Exception as e:
        print(f"[ERROR] Failed to fetch PRIDE metadata: {e}")
        return False

    # Save metadata
    metadata_path = Path("configs") / f"{accession}_metadata.json"
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=2)
    print(f"Saved project metadata to {metadata_path}")

    # Find cross-references
    cross_refs = metadata.get("crossReferences", [])
    other_links = metadata.get("otherOmicsLinks", [])
    print(f"[PRIDE] Found cross references: {cross_refs}")
    print(f"[PRIDE] Found other omics links: {other_links}")

    # 2. Fetch file list
    files_url = f"https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{accession}/files"
    print(f"[PRIDE] Fetching file list for {accession}...")
    try:
        resp = requests.get(files_url, headers={"Accept": "application/json"}, timeout=30)
        resp.raise_for_status()
        files_data = resp.json()
    except Exception as e:
        print(f"[ERROR] Failed to fetch file list: {e}")
        return False

    # Categorize files
    raw_files = []
    search_files = []
    sdrf_files = []

    for item in files_data:
        filename = item.get("fileName", "")
        file_category = item.get("fileCategory", {}).get("value", "")
        
        # Determine URL
        ftp_url = None
        for loc in item.get("publicFileLocations", []):
            if "FTP" in loc.get("name", ""):
                ftp_url = loc.get("value", "")
                break
        
        if not ftp_url:
            continue

        if filename.lower().endswith(".raw") or file_category == "RAW":
            raw_files.append((filename, ftp_url))
        elif filename.lower().endswith("sdrf.tsv") or filename.lower().endswith("sdrf.txt") or "sdrf" in filename.lower():
            sdrf_files.append((filename, ftp_url))
        elif file_category == "SEARCH" or filename.lower().endswith((".zip", ".rar", ".7z", "msms.txt", "peptides.txt")):
            search_files.append((filename, ftp_url))

    print(f"Categorized files: {len(raw_files)} RAW, {len(search_files)} SEARCH, {len(sdrf_files)} SDRF")

    # Download SDRF files
    for filename, ftp_url in sdrf_files:
        download_file(ftp_url, Path("configs") / filename)

    # Download SEARCH results
    for filename, ftp_url in search_files:
        dest = Path(psm_dir) / filename
        if download_file(ftp_url, dest):
            if filename.lower().endswith(".zip"):
                extract_zip(dest, Path(psm_dir))

    # Download RAW files
    if not skip_raw:
        dest_dir = Path(raw_dir) / accession
        dest_dir.mkdir(parents=True, exist_ok=True)
        for filename, ftp_url in raw_files:
            download_file(ftp_url, dest_dir / filename)
    else:
        print("Skipping raw file download as requested (--skip-raw).")

    return True

def process_massive(accession, raw_dir, psm_dir, skip_raw):
    """Acquires data from MassIVE."""
    # Standard MassIVE FTP Pattern
    ftp_root = f"ftp://massive.ucsd.edu/{accession}/"
    raw_folder = f"{ftp_root}raw/"
    
    print(f"[MASSIVE] Accession {accession} -> FTP root: {ftp_root}")
    
    # Download SEARCH results (by looking for zip files, peptides.txt, msms.txt, etc.)
    # We can use recursive wget with specific accepts
    dest_psm = Path(psm_dir) / accession
    dest_psm.mkdir(parents=True, exist_ok=True)
    
    print("[MASSIVE] Downloading search results/metadata...")
    cmd_search = [
        "wget", "-r", "-np", "-nd", "-l", "10",
        "-A", "*.zip,*.rar,*.txt,*.tsv",
        "-P", str(dest_psm),
        "-c", "--no-passive-ftp",
        ftp_root
    ]
    subprocess.run(cmd_search, check=False)
    
    # Unzip any downloaded zip files
    for zip_file in dest_psm.glob("*.zip"):
        extract_zip(zip_file, Path(psm_dir))

    # Download RAW files
    if not skip_raw:
        dest_raw = Path(raw_dir) / accession
        dest_raw.mkdir(parents=True, exist_ok=True)
        print("[MASSIVE] Downloading RAW files...")
        cmd_raw = [
            "wget", "-r", "-np", "-nd", "-l", "10",
            "-A", ".raw,.RAW",
            "-P", str(dest_raw),
            "-c", "--no-passive-ftp",
            raw_folder
        ]
        res = subprocess.run(cmd_raw, check=False)
        if res.returncode != 0:
            print("[WARNING] Failed to download from raw/ subdirectory. Trying root FTP folder...")
            cmd_raw[-1] = ftp_root
            subprocess.run(cmd_raw, check=False)
    else:
        print("Skipping raw file download as requested (--skip-raw).")

    return True

def process_ega(accession, raw_dir, credentials_path):
    """Acquires data from EGA (requires approved credentials)."""
    if not credentials_path or not Path(credentials_path).exists():
        print(f"[ERROR] EGA dataset {accession} requires approved credentials via --ega-credentials")
        return False
        
    pyega3_bin = shutil.which("pyega3")
    if not pyega3_bin:
        print("[ERROR] pyega3 is not installed or not in PATH. Please run: pip install pyega3")
        return False

    print(f"[EGA] Downloading {accession} using pyega3...")
    dest_dir = Path(raw_dir) / accession
    dest_dir.mkdir(parents=True, exist_ok=True)
    
    cmd = [pyega3_bin, "-cf", str(credentials_path), "fetch", accession, "--output-dir", str(dest_dir)]
    try:
        subprocess.run(cmd, check=True)
        print(f"[SUCCESS] EGA data retrieved.")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] pyega3 download failed: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Smart Data Acquisition Tool for Neoepitope Pipeline")
    parser.add_argument("-a", "--accession", help="Dataset Accession ID (e.g., PXD005231, MSV000085836, EGAD00001003XXX)")
    parser.add_argument("--skip-raw", action="store_true", help="Skip downloading raw files (download metadata/search files only)")
    parser.add_argument("--local", action="store_true", help="Local mode, bypass download checks")
    parser.add_argument("--raw-dir", default="data/raw", help="Path to raw data directory")
    parser.add_argument("--psm-dir", default="data/psms", help="Path to PSM/search results directory")
    parser.add_argument("--ega-credentials", help="Path to credentials JSON for EGA (pyega3)")
    
    args = parser.parse_args()
    
    if args.local:
        print("[LOCAL MODE] Bypassing online acquisition.")
        if args.accession:
            local_raw = Path(args.raw_dir) / args.accession
            local_psm = Path(args.psm_dir) / args.accession
        else:
            local_raw = Path(args.raw_dir)
            local_psm = Path(args.psm_dir)
            
        print(f"Checking raw directory: {local_raw}")
        print(f"Checking PSM directory: {local_psm}")
        
        # Verify if raw files exist
        raw_files = list(local_raw.glob("**/*.raw")) + list(local_raw.glob("**/*.RAW"))
        print(f"Found {len(raw_files)} local raw files.")
        return True

    if not args.accession:
        print("[ERROR] Accession ID (-a/--accession) is required unless running in --local mode.")
        sys.exit(1)

    acc = args.accession.upper()
    success = False
    
    if acc.startswith("PXD"):
        success = process_pride(acc, args.raw_dir, args.psm_dir, args.skip_raw)
    elif acc.startswith("MSV"):
        success = process_massive(acc, args.raw_dir, args.psm_dir, args.skip_raw)
    elif acc.startswith("EGAD"):
        success = process_ega(acc, args.raw_dir, args.ega_credentials)
    else:
        print(f"[ERROR] Unknown accession ID format: {args.accession}")
        sys.exit(1)
        
    if not success:
        print(f"[ERROR] Failed to acquire data for {args.accession}")
        sys.exit(1)
        
    print(f"[SUCCESS] Data acquisition complete for {args.accession}")

if __name__ == "__main__":
    main()
