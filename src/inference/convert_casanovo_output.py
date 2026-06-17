#!/usr/bin/env python3
"""
Convert Casanovo .mztab output → pipeline de_novo_candidates.tsv format.

Usage:
    python3 src/inference/convert_casanovo_output.py \
        --input  results/casanovo_all.mztab \
        --output results/de_novo_candidates_casanovo.tsv \
        --manifest configs/sample_manifest.tsv
"""

import argparse
from pathlib import Path
import pandas as pd


def parse_mztab(path: Path) -> pd.DataFrame:
    """Parse PSM section of a Casanovo mzTab file."""
    rows, header = [], None
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("PSH"):
                header = line.split("\t")[1:]
            elif line.startswith("PSM") and header:
                fields = line.split("\t")[1:]
                rows.append(dict(zip(header, fields)))
    if not rows:
        raise ValueError(f"No PSM rows found in {path}. Is this a valid Casanovo mzTab?")
    return pd.DataFrame(rows)


def _infer_sample_id(source_file: str, manifest: pd.DataFrame | None) -> str:
    if manifest is None or "filename" not in manifest.columns:
        return ""
    stem = Path(source_file).stem.lower()
    for _, row in manifest.iterrows():
        fn = str(row.get("filename", "")).lower().replace(".raw.txt", "")
        if stem in fn or fn in stem:
            return str(row.get("patient_id", ""))
    return ""


def convert(mztab_path: Path, manifest_path: Path | None, output_path: Path):
    print(f"Parsing {mztab_path} ...")
    raw = parse_mztab(mztab_path)
    print(f"  {len(raw):,} PSM rows | columns: {list(raw.columns)}")

    manifest = None
    if manifest_path and manifest_path.exists():
        manifest = pd.read_csv(manifest_path, sep="\t", dtype=str).fillna("")

    out = pd.DataFrame()

    # Peptide sequence
    seq_col = "sequence" if "sequence" in raw.columns else raw.columns[0]
    out["peptide"] = raw[seq_col].str.upper().str.replace(r"[^ACDEFGHIKLMNPQRSTVWY]", "", regex=True)

    # Score — pick best available column
    for sc in ["search_engine_score[1]", "opt_global_score", "opt_global_de_novo_score"]:
        if sc in raw.columns:
            out["score"] = pd.to_numeric(raw[sc], errors="coerce").fillna(0.0)
            break
    else:
        out["score"] = 0.0

    # Spectrum reference
    for rc in ["spectra_ref", "opt_global_spectrum_id"]:
        if rc in raw.columns:
            out["spectrum_id"] = raw[rc]
            break
    else:
        out["spectrum_id"] = ""

    # Source file
    for fc in ["opt_global_filename", "opt_global_source_file"]:
        if fc in raw.columns:
            out["source_file"] = raw[fc]
            break
    else:
        out["source_file"] = out["spectrum_id"]

    # Charge
    out["charge"] = pd.to_numeric(raw.get("charge", 2), errors="coerce").fillna(2).astype(int)

    # Map to sample_id
    out["sample_id"] = out["source_file"].apply(lambda f: _infer_sample_id(str(f), manifest))

    # Drop very short peptides
    before = len(out)
    out = out[out["peptide"].str.len() >= 8].reset_index(drop=True)
    print(f"  Kept {len(out):,}/{before:,} rows (length >= 8 AA)")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_path, sep="\t", index=False)
    print(f"Saved → {output_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, default=None)
    args = parser.parse_args()
    convert(args.input, args.manifest, args.output)


if __name__ == "__main__":
    main()
