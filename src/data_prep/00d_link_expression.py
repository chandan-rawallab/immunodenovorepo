#!/usr/bin/env python3
"""
Task 4: Expression Data Linker
Links PRIDE datasets to their corresponding RNA-seq expression data if available.
Mock TPM profiles are only generated when debug mode is explicitly enabled.
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import os
import json
import hashlib

def get_args():
    parser = argparse.ArgumentParser(description="Expression Data Linker")
    parser.add_argument("--accession", help="PRIDE Accession ID", required=False)
    parser.add_argument("--expression", help="User-provided expression file (override)", default=None)
    parser.add_argument("--expression-source", help="RNA source label to write when --expression is provided", default="provided_override")
    parser.add_argument("--debug-expression", action="store_true", help="Allow mock TPM generation when real expression is unavailable")
    parser.add_argument("--manifest", help="Sample manifest TSV", default="configs/sample_manifest.tsv")
    parser.add_argument("--output", help="Output normalized expression matrix (legacy/dummy)", default="data/expression_matrix.tsv")
    parser.add_argument("--fasta", help="Reference FASTA to extract gene IDs (used only for debug TPM generation)", default="data/reference/uniprot_human_reviewed.fasta")
    return parser.parse_args()

def extract_genes_from_fasta(fasta_path: Path) -> list:
    """Extract UniProt IDs from the reference FASTA file to serve as gene identifiers."""
    genes = []
    if fasta_path.exists():
        with open(fasta_path, "r") as f:
            for line in f:
                if line.startswith(">"):
                    # e.g., >sp|P0DTC2|SPIKE_SARS2 ...
                    parts = line.split("|")
                    if len(parts) >= 2:
                        genes.append(parts[1])
    # Fallback to some common genes if FASTA doesn't exist or parsing failed
    if not genes:
        genes = ["P0DTC2", "P53", "KRAS", "EGFR", "BRCA1", "MYC", "PTEN"]
    return list(set(genes))

def generate_mock_expression(patient_id: str, genes: list, out_path: Path):
    """Generate a mock RNA-seq TPM profile for a patient.

    This is for explicit debug runs only.
    """
    # Use a log-normal distribution to simulate TPM values (many low, few high)
    # Mean of log-normal = exp(mu + sigma^2 / 2). 
    seed = int(hashlib.sha256(patient_id.encode("utf-8")).hexdigest()[:8], 16)
    rng = np.random.default_rng(seed)
    tpms = rng.lognormal(mean=0.5, sigma=1.2, size=len(genes))
    
    df = pd.DataFrame({
        "gene": genes,
        "expression_tpm": tpms
    })
    
    # Ensure some genes are explicitly highly expressed for testing pipeline
    # (Set 10% of genes to high TPM)
    high_expr_idx = rng.choice(len(genes), size=max(1, len(genes)//10), replace=False)
    df.loc[high_expr_idx, "expression_tpm"] = rng.uniform(10.0, 500.0, size=len(high_expr_idx))
    
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"  Generated mock expression profile for {patient_id}: {out_path}")

def main():
    args = get_args()
    manifest_path = Path(args.manifest)
    
    if not manifest_path.exists():
        print(f"[ERROR] Manifest not found at {manifest_path}")
        return

    print(f"Loading manifest from {manifest_path}...")
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str).fillna("")
    
    if "patient_id" not in manifest.columns:
        print("[ERROR] Manifest lacks 'patient_id' column.")
        return

    # Ensure RNA columns exist
    if "rna_expr_path" not in manifest.columns:
        manifest["rna_expr_path"] = ""
    if "rna_source" not in manifest.columns:
        manifest["rna_source"] = ""

    expression_dir = Path("data/expression")
    expression_dir.mkdir(parents=True, exist_ok=True)

    updated = False
    patient_paths: dict[str, tuple[str, str]] = {}
    for patient_id, group in manifest.groupby("patient_id", sort=False):
        existing_paths = [p for p in group["rna_expr_path"].astype(str) if p and Path(p).exists()]
        if existing_paths:
            expr_path = existing_paths[0]
            existing_sources = group["rna_source"].replace("", pd.NA).dropna()
            source = str(existing_sources.iloc[0]) if not existing_sources.empty else "existing"
            print(f"  Patient {patient_id} already has expression data at {expr_path}")
        elif args.expression and Path(args.expression).exists():
            expr_file = expression_dir / f"{patient_id}_tpm.tsv"
            print(f"  Using provided expression file for {patient_id}.")
            import shutil
            shutil.copy(args.expression, expr_file)
            expr_path = str(expr_file)
            source = str(args.expression_source).strip() or "provided_override"
        elif args.debug_expression:
            genes = extract_genes_from_fasta(Path(args.fasta))
            print(f"Extracted {len(genes)} unique protein/gene IDs for debug expression profiling.")
            expr_file = expression_dir / f"{patient_id}_tpm.tsv"
            generate_mock_expression(patient_id, genes, expr_file)
            expr_path = str(expr_file)
            source = "mock_debug"
        else:
            expr_path = ""
            source = "missing"
            print(f"  No patient-matched expression provided for {patient_id}; leaving RNA unlinked.")
        patient_paths[patient_id] = (expr_path, source)

    for idx, row in manifest.iterrows():
        expr_path, source = patient_paths[row["patient_id"]]
        if row.get("rna_expr_path", "") != expr_path:
            manifest.at[idx, "rna_expr_path"] = expr_path
            updated = True
        if row.get("rna_source", "") != source:
            manifest.at[idx, "rna_source"] = source
            updated = True

    if updated:
        manifest.to_csv(manifest_path, sep="\t", index=False)
        print(f"Updated manifest with RNA-seq expression paths saved to {manifest_path}")

    # Write a dummy legacy expression matrix to satisfy pipeline args if needed
    out_matrix = Path(args.output)
    out_matrix.parent.mkdir(parents=True, exist_ok=True)
    out_matrix.write_text("dummy\tmatrix\n")
    print(f"Created dummy matrix at {out_matrix}")

if __name__ == "__main__":
    main()
