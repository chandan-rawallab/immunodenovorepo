#!/usr/bin/env python3
"""
Map CCLE (Cancer Cell Line Encyclopedia) RNA-seq expression to UniProt IDs,
producing per-patient expression TSVs that replace the mock_debug profiles.

This script uses surrogate LCL (lymphoblastoid) expression for B-cell line
patients and skin/melanoma expression for TIL patients.

Usage:
    python3 src/data_prep/00e_map_ccle_to_uniprot.py \
        --ccle data/expression/ccle_raw/CCLE_expression.csv \
        --sample-info data/expression/ccle_raw/sample_info.csv \
        --uniprot-map data/expression/ccle_raw/HUMAN_9606_idmapping_selected.tab.gz \
        --manifest configs/sample_manifest.tsv \
        --output-dir data/expression/

Download the required files first:
    wget -O data/expression/ccle_raw/CCLE_expression.csv \
      "https://depmap.org/portal/download/api/downloads?file_name=OmicsExpressionProteinCodingGenesTPMLogp1.csv"
    wget -O data/expression/ccle_raw/sample_info.csv \
      "https://depmap.org/portal/download/api/downloads?file_name=Model.csv"
    wget -O data/expression/ccle_raw/HUMAN_9606_idmapping_selected.tab.gz \
      "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz"
"""

import argparse
import gzip
import re
import sys
from pathlib import Path

import pandas as pd
import numpy as np


# Patient-type classification (must match sample_manifest.tsv patient_id values)
_BCELL_PATIENTS = {"GD149", "PD42", "RA957", "MD155", "CM467"}
_TIL_PATIENTS   = {"TIL1", "TIL3", "pooledTIL3"}
_CTRL_PATIENTS  = {"Apher-1", "Apher-6"}   # will use LCL surrogate too


def get_args():
    p = argparse.ArgumentParser(description="Map CCLE expression to UniProt IDs.")
    p.add_argument("--ccle", type=Path, required=True,
                   help="CCLE expression CSV (OmicsExpressionProteinCodingGenesTPMLogp1.csv)")
    p.add_argument("--sample-info", type=Path, required=True,
                   help="CCLE model metadata CSV (Model.csv)")
    p.add_argument("--uniprot-map", type=Path, required=True,
                   help="UniProt ID mapping file (HUMAN_9606_idmapping_selected.tab.gz)")
    p.add_argument("--manifest", type=Path, default=Path("configs/sample_manifest.tsv"))
    p.add_argument("--output-dir", type=Path, default=Path("data/expression"))
    p.add_argument("--dry-run", action="store_true",
                   help="Parse and map only; do not write files or update manifest")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Step 1: Build gene-symbol → UniProt ID mapping
# ---------------------------------------------------------------------------

def build_symbol_to_uniprot(mapping_path: Path) -> dict[str, str]:
    """
    Read HUMAN_9606_idmapping_selected.tab.gz.
    Columns: UniProtKB_AC, UniProtKB_ID, GeneID, RefSeq, GI, PDB, GO,
             UniRef100, UniRef90, UniRef50, UniParc, PIR, NCBI_taxon,
             MIM, UniGene, PubMed, EMBL, EMBL_CDS, Ensembl_G, Ensembl_T,
             Ensembl_P, Uniprot_isoforms
    We need columns 0 (UniProt AC) and 14 (Gene_Name / gene symbol).
    The selected tab file uses: UniProtAC, UniProtID(entry_name), GeneID, ...
    Actual column for gene symbol is index 2 in the *full* mapping file,
    but in the *selected* file it is not directly present. We parse from
    UniProtID column (format: GENENAME_HUMAN).
    """
    symbol_to_uniprot: dict[str, str] = {}
    open_fn = gzip.open if str(mapping_path).endswith(".gz") else open

    with open_fn(mapping_path, "rt", errors="replace") as fh:
        for line in fh:
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            uniprot_ac = parts[0].strip()
            entry_name = parts[1].strip()           # e.g. TP53_HUMAN
            if entry_name.endswith("_HUMAN"):
                symbol = entry_name[: -len("_HUMAN")]
                # Keep reviewed (Swiss-Prot) entry when multiple map to same symbol
                if symbol not in symbol_to_uniprot:
                    symbol_to_uniprot[symbol] = uniprot_ac

    print(f"  Built gene-symbol → UniProt map: {len(symbol_to_uniprot):,} entries")
    return symbol_to_uniprot


# ---------------------------------------------------------------------------
# Step 2: Load CCLE expression and compute median profiles
# ---------------------------------------------------------------------------

_GENE_COL_RE = re.compile(r"^([A-Z0-9_\-\.]+)\s*\(ENSG")  # "TP53 (ENSG00000141510)"


def load_ccle_and_compute_profiles(
    ccle_path: Path, sample_info_path: Path
) -> dict[str, pd.Series]:
    """
    Returns a dict with two keys:
        "lcl"      → median log2(TPM+1) profile across LCL/lymphoblast lines
        "melanoma" → median log2(TPM+1) profile across skin/melanoma lines

    Values are pandas Series with gene symbols as index, TPM (linear) as values.
    """
    print("Loading CCLE sample info ...")
    meta = pd.read_csv(sample_info_path, low_memory=False)
    # Column names vary between CCLE releases — try common ones
    lineage_col = next(
        (c for c in ["lineage", "Lineage", "primary_disease", "OncotreeLineage"] if c in meta.columns),
        None,
    )
    id_col = next(
        (c for c in ["DepMap_ID", "ModelID", "model_id"] if c in meta.columns),
        None,
    )
    if lineage_col is None or id_col is None:
        print(f"  Available columns in sample_info: {list(meta.columns)}")
        sys.exit("[ERROR] Cannot find lineage or ID column in CCLE sample info CSV.")

    lcl_ids = set(
        meta.loc[meta[lineage_col].str.lower().str.contains("lymph|leukemia|b_cell|lcl", na=False), id_col]
    )
    mel_ids = set(
        meta.loc[meta[lineage_col].str.lower().str.contains("skin|melanoma", na=False), id_col]
    )
    print(f"  LCL-type cell lines: {len(lcl_ids)} | Skin/melanoma: {len(mel_ids)}")

    print("Loading CCLE expression matrix (this may take ~1 min) ...")
    expr = pd.read_csv(ccle_path, index_col=0, low_memory=False)
    # rows = DepMap_IDs, cols = "SYMBOL (ENSGXXX)"
    print(f"  Expression matrix shape: {expr.shape}")

    # Parse gene symbol from column name
    new_cols = []
    for c in expr.columns:
        m = _GENE_COL_RE.match(str(c))
        new_cols.append(m.group(1) if m else str(c).split()[0])
    expr.columns = new_cols

    # Convert log2(TPM+1) → TPM
    expr_tpm = (2 ** expr) - 1

    profiles = {}
    for label, id_set in [("lcl", lcl_ids), ("melanoma", mel_ids)]:
        overlap = expr_tpm.index.intersection(list(id_set))
        if len(overlap) == 0:
            print(f"  WARNING: No {label} cell lines found in expression matrix.")
            profiles[label] = pd.Series(dtype=float)
        else:
            profiles[label] = expr_tpm.loc[overlap].median(axis=0)
            print(f"  {label} median profile: {len(overlap)} lines, {len(profiles[label])} genes")

    return profiles


# ---------------------------------------------------------------------------
# Step 3: Write per-patient expression TSVs
# ---------------------------------------------------------------------------

def write_patient_tpm(
    patient_id: str,
    profile: pd.Series,
    symbol_to_uniprot: dict,
    output_dir: Path,
    source_label: str,
) -> tuple[str, str]:
    """
    Write data/expression/{patient_id}_tpm.tsv.
    Returns (file_path, rna_source_label).
    """
    out_path = output_dir / f"{patient_id}_tpm.tsv"

    rows = []
    for symbol, tpm in profile.items():
        uniprot = symbol_to_uniprot.get(str(symbol))
        if uniprot:
            rows.append({"gene": uniprot, "expression_tpm": float(tpm)})

    df = pd.DataFrame(rows)
    if df.empty:
        print(f"  WARNING: No UniProt IDs mapped for {patient_id} — skipping")
        return "", "mapping_failed"

    df = df.dropna().groupby("gene", as_index=False)["expression_tpm"].max()
    df.to_csv(out_path, sep="\t", index=False)
    print(f"  Wrote {len(df):,} UniProt entries → {out_path}")
    return str(out_path), source_label


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = get_args()

    for path in [args.ccle, args.sample_info, args.uniprot_map, args.manifest]:
        if not path.exists():
            sys.exit(f"[ERROR] Required file not found: {path}\n"
                     "Run the download commands in this script's docstring first.")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("=== Step 1: Build gene-symbol → UniProt mapping ===")
    symbol_to_uniprot = build_symbol_to_uniprot(args.uniprot_map)

    print("\n=== Step 2: Load CCLE expression ===")
    profiles = load_ccle_and_compute_profiles(args.ccle, args.sample_info)

    print("\n=== Step 3: Read manifest ===")
    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str).fillna("")
    manifest.columns = manifest.columns.str.strip()
    manifest = manifest.apply(lambda col: col.str.strip() if col.dtype == "object" else col)
    patient_ids = manifest["patient_id"].unique()
    print(f"  Unique patients: {list(patient_ids)}")

    print("\n=== Step 4: Write per-patient expression files ===")
    path_map: dict[str, str] = {}
    source_map: dict[str, str] = {}

    for pid in patient_ids:
        if pid in _BCELL_PATIENTS or pid in _CTRL_PATIENTS:
            profile = profiles.get("lcl", pd.Series(dtype=float))
            source_label = "ccle_lcl_surrogate"
        elif pid in _TIL_PATIENTS:
            profile = profiles.get("melanoma", pd.Series(dtype=float))
            source_label = "ccle_melanoma_surrogate"
        else:
            print(f"  Unknown patient type for {pid} — using LCL profile")
            profile = profiles.get("lcl", pd.Series(dtype=float))
            source_label = "ccle_lcl_surrogate"

        if profile.empty:
            print(f"  Skipping {pid}: no profile data available")
            path_map[pid] = ""
            source_map[pid] = "no_data"
            continue

        if args.dry_run:
            print(f"  [DRY RUN] Would write {pid}_tpm.tsv ({source_label})")
            path_map[pid] = str(args.output_dir / f"{pid}_tpm.tsv")
            source_map[pid] = source_label
        else:
            file_path, src = write_patient_tpm(pid, profile, symbol_to_uniprot, args.output_dir, source_label)
            path_map[pid] = file_path
            source_map[pid] = src if file_path else "mapping_failed"

    if not args.dry_run:
        print("\n=== Step 5: Update manifest ===")
        for idx, row in manifest.iterrows():
            pid = row["patient_id"]
            if pid in path_map and path_map[pid]:
                manifest.at[idx, "rna_expr_path"] = path_map[pid]
                manifest.at[idx, "rna_source"] = source_map[pid]

        manifest.to_csv(args.manifest, sep="\t", index=False)
        print(f"  Manifest updated: {args.manifest}")

        # Verify no more mock_debug rows
        remaining_mock = (manifest["rna_source"] == "mock_debug").sum()
        if remaining_mock == 0:
            print("  ✅ All rna_source entries replaced (no more mock_debug)")
        else:
            print(f"  ⚠️  {remaining_mock} rows still have rna_source=mock_debug")

    print("\n=== Done ===")
    print("IMPORTANT: This data uses CCLE surrogate profiles, not patient-matched RNA-seq.")
    print("Disclose this in the paper Methods section as surrogate expression data.")


if __name__ == "__main__":
    main()
