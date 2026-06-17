#!/usr/bin/env python3
"""
Task 2: Smart Manifest Builder
Constructs configs/sample_manifest.tsv from SDRF metadata files, filename regex fallbacks, or manual overrides.
"""

import os
import re
import sys
import argparse
import pandas as pd
from pathlib import Path

def normalize_hla(value: str) -> str:
    """Normalize HLA allele to MHCflurry format (e.g., HLA-A*02:01)."""
    if not isinstance(value, str):
        return "TBD"
    value = value.strip()
    if not value or value.upper() in ["TBD", "UNKNOWN", "MISSING", "NAN", "NULL"]:
        return "TBD"
    
    # If starts with A*, B*, C*, prepend HLA-
    if re.match(r'^[A-C]\*\d{2}:\d{2}', value, re.IGNORECASE):
        return f"HLA-{value.upper()}"
        
    # If already HLA-A*02:01, just upper it
    if re.match(r'^HLA-[A-C]\*\d{2}:\d{2}', value, re.IGNORECASE):
        return value.upper()
        
    # Standardize basic format: remove HLA-, asterisks, colons
    clean = value.upper().replace("HLA-", "").replace("*", "").replace(":", "").strip()
    
    # Check if we have A0201 format (letter + 4 digits)
    match = re.match(r'^([A-C])(\d{2})(\d{2})$', clean)
    if match:
        return f"HLA-{match.group(1)}*{match.group(2)}:{match.group(3)}"
        
    # If it's just A*0201
    match = re.match(r'^([A-C])\*?(\d{2})(\d{2})$', clean)
    if match:
        return f"HLA-{match.group(1)}*{match.group(2)}:{match.group(3)}"
        
    return value

def parse_hla_list(hla_str: str) -> str:
    """Parse and normalize a list of HLA alleles (comma or semicolon separated)."""
    if not isinstance(hla_str, str) or hla_str.strip().upper() in ["TBD", "UNKNOWN", "", "NAN", "NULL"]:
        return "TBD"
    # Split by semicolon, comma, slash or whitespace
    alleles = re.split(r'[;,/\s]+', hla_str)
    normalized = [normalize_hla(a) for a in alleles if a.strip()]
    valid = [n for n in normalized if n != "TBD"]
    if not valid:
        return "TBD"
    return ",".join(sorted(list(set(valid))))

def find_column(df, patterns):
    """Find a column in a DataFrame matching a list of regex patterns."""
    for pattern in patterns:
        for col in df.columns:
            if re.search(pattern, str(col), re.IGNORECASE):
                return col
    return None

def extract_patient_from_filename(filename):
    """Extract patient ID from filename using common heuristics."""
    # Remove file extension and prefix/suffix structures like "msms_" or "combined"
    base = os.path.basename(filename)
    base_clean = re.sub(r'^(msms_|_|combined_)', '', base)
    base_clean = re.sub(r'\.(raw|txt|tsv|mgf|mgf\.txt|raw\.txt)$', '', base_clean, flags=re.IGNORECASE)
    
    # Try common patterns first
    match = re.search(r'\b(TIL\d+|GD\d+|PD\d+|Apher-\d+|CM\d+|RA\d+|MD\d+|pooledTIL\d+)\b', base_clean, re.IGNORECASE)
    if match:
        return match.group(1)
    
    # Try general Patient_A or PatientA patterns
    match = re.search(r'patient[_\s-]?([A-Z0-9]+)', base_clean, re.IGNORECASE)
    if match:
        return f"Patient_{match.group(1)}"
        
    # Try splitting by underscores/dashes and look for common identifier patterns
    parts = re.split(r'[_.-]', base_clean)
    for part in parts:
        if re.match(r'^(TIL|GD|PD|Apher|CM|RA|MD|Patient|pooledTIL)\d+$', part, re.IGNORECASE):
            return part
            
    # Default fallback: return the file base name or part of it
    if "_" in base_clean:
        return base_clean.split("_")[0]
    return base_clean

def build_manifest(accession, raw_dir, psm_dir, sdrf_path, hla_override, output_file):
    print(f"Building manifest for accession: {accession}")
    print(f"RAW directory: {raw_dir}")
    print(f"PSM directory: {psm_dir}")
    
    # 1. Discover all raw runs/files
    # Scan raw_dir and psm_dir for files
    raw_files = []
    run_ids = set()
    
    # Check RAW files (.raw, .RAW)
    if os.path.exists(raw_dir):
        for root, _, files in os.walk(raw_dir):
            for f in files:
                if f.lower().endswith(".raw"):
                    full_path = os.path.join(root, f)
                    run_id = f[:-4]
                    raw_files.append((run_id, f, full_path, "raw"))
                    run_ids.add(run_id)
                    
    # Check PSM files (msms_*.txt)
    if os.path.exists(psm_dir):
        for root, _, files in os.walk(psm_dir):
            for f in files:
                if f.startswith("msms_") and f.endswith(".txt"):
                    run_id = f.replace("msms_", "").replace(".raw.txt", "").replace(".txt", "")
                    if run_id not in run_ids:
                        full_path = os.path.join(root, f)
                        # Assume matching raw filename
                        raw_files.append((run_id, f"msms_{run_id}.raw.txt", full_path, "psm"))
                        run_ids.add(run_id)

    if not raw_files:
        print("WARNING: No raw or PSM files discovered. Checking configs/ for existing manifest context.")
    
    # 2. Try parsing SDRF
    sdrf_df = None
    sdrf_file_found = None
    
    # Search for SDRF file if not explicitly provided
    if not sdrf_path:
        search_dirs = [Path("configs"), Path(raw_dir), Path(psm_dir)]
        for d in search_dirs:
            if d.exists():
                sdrfs = list(d.glob("**/*sdrf*.tsv")) + list(d.glob("**/*sdrf*.txt"))
                if sdrfs:
                    sdrf_file_found = sdrfs[0]
                    print(f"Found SDRF file: {sdrf_file_found}")
                    break
    else:
        if os.path.exists(sdrf_path):
            sdrf_file_found = Path(sdrf_path)
            print(f"Using specified SDRF file: {sdrf_file_found}")
            
    if sdrf_file_found:
        try:
            sdrf_df = pd.read_csv(sdrf_file_found, sep='\t')
            print(f"Successfully loaded SDRF with {len(sdrf_df)} rows and columns: {list(sdrf_df.columns)}")
        except Exception as e:
            print(f"ERROR: Failed to load SDRF file {sdrf_file_found}: {e}")

    sdrf_mappings = {}
    if sdrf_df is not None:
        # Find columns
        file_col = find_column(sdrf_df, [r'comment\[file\s*name\]', r'comment\[data\s*file\]', r'file\s*name', r'filename'])
        patient_col = find_column(sdrf_df, [r'source\s*name', r'characteristics\[individual\]', r'sample\s*name', r'patient', r'characteristics\[organism\s*part\]'])
        hla_col = find_column(sdrf_df, [r'characteristics\[hla\s*allele\]', r'comment\[hla\s*allele\]', r'hla\s*allele', r'hla'])
        
        print(f"SDRF Columns identified: File={file_col}, Patient={patient_col}, HLA={hla_col}")
        
        if file_col and patient_col:
            for _, row in sdrf_df.iterrows():
                fname = str(row[file_col])
                # Clean filename
                fname_clean = os.path.basename(fname)
                run_id = fname_clean.replace(".raw", "").replace(".RAW", "").replace(".raw.txt", "")
                
                patient_id = str(row[patient_col])
                hla_val = str(row[hla_col]) if hla_col else "TBD"
                
                sdrf_mappings[run_id] = {
                    'patient_id': patient_id,
                    'hla_alleles': parse_hla_list(hla_val),
                    'hla_source': 'sdrf' if hla_val and hla_val.strip().upper() not in ["TBD", "UNKNOWN", ""] else 'missing'
                }
        else:
            print("WARNING: Could not identify mapping columns (File and Patient) in SDRF.")

    # 3. Handle HLA overrides
    hla_overrides = {}
    if hla_override and os.path.exists(hla_override):
        try:
            override_df = pd.read_csv(hla_override, sep='\t')
            # Expecting columns: patient_id, hla_alleles (or similar)
            p_col = find_column(override_df, [r'patient', r'sample'])
            h_col = find_column(override_df, [r'hla'])
            if p_col and h_col:
                for _, row in override_df.iterrows():
                    pat = str(row[p_col]).strip()
                    hla = str(row[h_col]).strip()
                    hla_overrides[pat] = parse_hla_list(hla)
                print(f"Loaded {len(hla_overrides)} HLA overrides from {hla_override}")
            else:
                print("WARNING: HLA override file missing 'patient_id' or 'hla_alleles' columns.")
        except Exception as e:
            print(f"ERROR: Failed to load HLA override file {hla_override}: {e}")

    # 4. Construct Manifest Rows
    manifest_rows = []
    
    # If we have discovered files
    if raw_files:
        for run_id, filename, full_path, file_type in raw_files:
            # Check SDRF mapping first
            if run_id in sdrf_mappings:
                patient_id = sdrf_mappings[run_id]['patient_id']
                hla_alleles = sdrf_mappings[run_id]['hla_alleles']
                hla_source = sdrf_mappings[run_id]['hla_source']
            else:
                # Heuristic pattern matching
                patient_id = extract_patient_from_filename(run_id)
                hla_alleles = "TBD"
                hla_source = "missing"
            
            # Apply HLA override if applicable
            if patient_id in hla_overrides:
                hla_alleles = hla_overrides[patient_id]
                hla_source = "manual"
                
            manifest_rows.append({
                'run_id': run_id,
                'patient_id': patient_id,
                'validation_id': patient_id,
                'filename': filename,
                'hla_alleles': hla_alleles,
                'hla_source': hla_source,
                'rna_expr_path': "",
                'rna_source': "missing",
                'study_id': accession,
                'cohort': accession,
                'sample_role': 'unknown',
                'raw_source': 'local_raw' if file_type == 'raw' else 'missing',
                'psm_source': 'local_psm' if file_type == 'psm' else 'missing',
                'include_in_pipeline': True,
                'notes': 'generated by smart manifest builder',
            })
    else:
        # Fallback to existing manifest if it exists and we're in local mode
        existing_manifest = Path("configs/sample_manifest.tsv")
        if existing_manifest.exists():
            print(f"Using existing manifest file as base: {existing_manifest}")
            try:
                m_df = pd.read_csv(existing_manifest, sep='\t')
                for _, row in m_df.iterrows():
                    run_id = row['run_id']
                    patient_id = row['patient_id']
                    filename = row.get('filename', f"msms_{run_id}.raw.txt")
                    validation_id = row.get('validation_id', patient_id)
                    sample_role = row.get('sample_role', 'unknown')
                    hla_alleles = row.get('hla_alleles', 'TBD')
                    hla_source = row.get('hla_source', 'manual' if hla_alleles != 'TBD' else 'missing')
                    rna_path = row.get('rna_expr_path', '')
                    rna_source = row.get('rna_source', 'missing')
                    cohort = row.get('cohort', accession)
                    study_id = row.get('study_id', cohort)
                    raw_source = row.get('raw_source', 'existing_manifest')
                    psm_source = row.get('psm_source', 'existing_manifest')
                    include_in_pipeline = row.get('include_in_pipeline', True)
                    notes = row.get('notes', 'loaded from existing manifest')
                    
                    if patient_id in hla_overrides:
                        hla_alleles = hla_overrides[patient_id]
                        hla_source = "manual"
                        
                    manifest_rows.append({
                        'run_id': run_id,
                        'patient_id': patient_id,
                        'validation_id': validation_id,
                        'filename': filename,
                        'hla_alleles': hla_alleles,
                        'hla_source': hla_source,
                        'rna_expr_path': rna_path,
                        'rna_source': rna_source,
                        'study_id': study_id,
                        'cohort': cohort,
                        'sample_role': sample_role,
                        'raw_source': raw_source,
                        'psm_source': psm_source,
                        'include_in_pipeline': include_in_pipeline,
                        'notes': notes,
                    })
            except Exception as e:
                print(f"ERROR reading existing manifest: {e}")

    # Build final DataFrame
    if not manifest_rows:
        print("ERROR: No runs found and no base manifest template could be loaded.")
        sys.exit(1)
        
    df_manifest = pd.DataFrame(manifest_rows)
    
    # Save manifest
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    df_manifest.to_csv(output_file, sep='\t', index=False)
    print(f"SUCCESS: Wrote sample manifest to {output_file} ({len(df_manifest)} runs).")
    
    # Print summary
    print("\nManifest Summary:")
    print(df_manifest[['patient_id', 'hla_source', 'hla_alleles']].drop_duplicates().to_string(index=False))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Smart Manifest Builder for Neoepitope Pipeline")
    parser.add_argument("-a", "--accession", default="PXD005231", help="Dataset accession ID")
    parser.add_argument("--raw-dir", default="data/raw", help="Directory where raw data is stored")
    parser.add_argument("--psm-dir", default="data/psms", help="Directory with search results")
    parser.add_argument("--sdrf", help="Optional path to specific SDRF file")
    parser.add_argument("--hla-override", help="Optional path to HLA override TSV (patient_id -> hla_alleles)")
    parser.add_argument("-o", "--output", default="configs/sample_manifest.tsv", help="Path to write final manifest")
    
    args = parser.parse_args()
    build_manifest(args.accession, args.raw_dir, args.psm_dir, args.sdrf, args.hla_override, args.output)
