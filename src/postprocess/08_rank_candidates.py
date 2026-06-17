#!/usr/bin/env python3
"""Activity 3: Ranking neoantigen candidates by binding affinity and expression."""

import argparse
import pandas as pd
from pathlib import Path
import subprocess
import tempfile
import shutil

SCRIPT_DIR = Path(__file__).resolve().parent
WORKSPACE_ROOT = SCRIPT_DIR.parent.parent
PATIENT_MATCHED_RNA_SOURCES = {
    "real_patient_matched",
    "patient_matched",
    "provided_patient_matched",
    "clinical_rna",
    "matched_rna",
}
NON_BIOLOGICAL_RNA_MARKERS = ("mock", "debug", "surrogate", "inferred", "missing", "provided_override")

def normalize_hla(value: str) -> str:
    """Normalize HLA allele to MHCflurry format (e.g., HLA-A*02:01)."""
    value = value.strip()
    if not value:
        return ""
    
    # If it's already in a good format, return it
    if "*" in value and ":" in value:
        return value
    
    # Handle A0201 or HLA-A0201 formats
    clean = value.replace("HLA-", "").replace("*", "").replace(":", "")
    if len(clean) >= 5 and clean[0].isalpha():
        # A0201 -> HLA-A*02:01
        return f"HLA-{clean[0]}*{clean[1:3]}:{clean[3:5]}"
    
    return value

def run_mhcflurry(peptides, alleles, output_path):
    """Run MHCflurry prediction for a list of peptides and alleles."""
    predictor = shutil.which("mhcflurry-predict")
    if not predictor:
        local_predictor = WORKSPACE_ROOT / ".venv" / "bin" / "mhcflurry-predict"
        if local_predictor.exists():
            predictor = str(local_predictor)
    if not predictor:
        print("Warning: mhcflurry-predict not found. Skipping binding prediction.")
        return False
    
    with tempfile.TemporaryDirectory() as tmpdir:
        input_csv = Path(tmpdir) / "input.csv"
        with open(input_csv, "w") as f:
            f.write("allele,peptide\n")
            for allele in alleles:
                for peptide in peptides:
                    f.write(f"{allele},{peptide}\n")
        
        try:
            subprocess.run([predictor, str(input_csv), "--out", str(output_path)], check=True, capture_output=True)
            return True
        except subprocess.CalledProcessError as e:
            print(f"Error running MHCflurry: {e.stderr.decode()}")
            return False

def normalize_protein_id(value: object) -> str:
    """Normalize accession-like IDs for expression lookup."""
    text = str(value or "").strip()
    if not text or text.lower() == "nan":
        return ""
    # UniProt isoforms often appear as Q16625-2 while expression may use Q16625.
    return text.split("-")[0]

def source_protein_ids(value: object) -> list[str]:
    """Split semicolon-delimited source protein annotations into lookup IDs."""
    text = str(value or "").strip()
    if not text or text.lower() == "nan":
        return []
    ids = []
    for part in text.split(";"):
        norm = normalize_protein_id(part)
        if norm and norm not in ids:
            ids.append(norm)
    return ids

def expression_for_source(source_value: object, rna_path: object, expression_lookup: dict[tuple[str, str], float]) -> float:
    """Return the max TPM among all source protein IDs for this patient RNA file."""
    path = str(rna_path or "")
    values = [
        expression_lookup[(path, protein_id)]
        for protein_id in source_protein_ids(source_value)
        if (path, protein_id) in expression_lookup
    ]
    return max(values) if values else 0.0

def collapse_duplicate_candidates(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse repeated sample+peptide rows while keeping the strongest evidence row."""
    if not {"sample_id", "peptide"}.issubset(df.columns):
        return df
    if not df.duplicated(["sample_id", "peptide"]).any():
        return df

    collapsed = []
    sort_cols = [c for c in ["binding_rank", "score"] if c in df.columns]
    for _, group in df.groupby(["sample_id", "peptide"], sort=False):
        ranked = group.copy()
        if "binding_rank" in ranked.columns:
            ranked["_binding_sort"] = pd.to_numeric(ranked["binding_rank"], errors="coerce").fillna(float("inf"))
            ranked = ranked.sort_values("_binding_sort", ascending=True)
        elif "score" in ranked.columns:
            ranked["_score_sort"] = pd.to_numeric(ranked["score"], errors="coerce").fillna(float("-inf"))
            ranked = ranked.sort_values("_score_sort", ascending=False)
        row = ranked.iloc[0].drop(labels=[c for c in ["_binding_sort", "_score_sort"] if c in ranked.columns]).copy()

        if "psm_count" in group.columns:
            row["psm_count"] = pd.to_numeric(group["psm_count"], errors="coerce").max()
        if "expression_tpm" in group.columns:
            row["expression_tpm"] = pd.to_numeric(group["expression_tpm"], errors="coerce").fillna(0.0).max()
        if "spectrum_id" in group.columns:
            spectra = [str(x) for x in group["spectrum_id"].dropna().unique()]
            row["spectrum_id"] = ";".join(spectra[:25])
        collapsed.append(row)

    return pd.DataFrame(collapsed).reset_index(drop=True)

def normalize_rna_source(value: object) -> str:
    """Normalize manifest RNA source labels for evidence gating."""
    text = str(value or "").strip().lower()
    if not text or text == "nan":
        return "missing"
    return text

def expression_evidence_status(value: object) -> str:
    """Classify whether RNA expression can support biological evidence classes."""
    source = normalize_rna_source(value)
    if source in PATIENT_MATCHED_RNA_SOURCES:
        return "patient_matched"
    if any(marker in source for marker in NON_BIOLOGICAL_RNA_MARKERS):
        return "debug_or_non_patient_matched"
    return "unverified"

def expression_can_support_biology(value: object) -> bool:
    return expression_evidence_status(value) == "patient_matched"

def main():
    parser = argparse.ArgumentParser(description="Rank neoantigen candidates.")
    parser.add_argument("--input", type=Path, required=True, help="Path to filtered_neoantigens.tsv")
    parser.add_argument("--manifest", type=Path, required=True, help="Path to sample_manifest.tsv")
    parser.add_argument("--output", type=Path, required=True, help="Output ranked TSV")
    parser.add_argument("--binding_rank_cutoff", type=float, default=2.0)
    parser.add_argument("--tpm_cutoff", type=float, default=1.0)
    parser.add_argument("--keep_duplicate_rows", action="store_true",
                        help="Keep repeated sample+peptide rows instead of collapsing them before evidence classification.")
    
    args = parser.parse_args()
    
    # 1. Load data
    df = pd.read_csv(args.input, sep="\t")
    manifest = pd.read_csv(args.manifest, sep="\t")
    
    # 2. Binding Prediction (MHCflurry)
    # Group by sample to get alleles
    results = []
    for sample_id, group in df.groupby("sample_id"):
        sample_meta = manifest[manifest['patient_id'] == sample_id] # Adjust if patient_id != sample_id
        if sample_meta.empty:
            sample_meta = manifest[manifest['run_id'] == sample_id]
            
        if sample_meta.empty:
            print(f"Warning: No metadata found for sample {sample_id}")
            continue
            
        # Get alleles from hla_alleles column or individual columns if they exist
        raw_alleles = str(sample_meta.iloc[0]['hla_alleles']).split(",")
        alleles = [normalize_hla(a) for a in raw_alleles if a and a.strip().lower() != "tbd"]
        
        if not alleles:
            print(f"Warning: No valid HLA alleles for sample {sample_id}")
            group_copy = group.copy()
            group_copy['best_hla'] = ""
            group_copy['binding_rank'] = None
            results.append(group_copy)
            continue
            
        peptides = group['peptide'].unique().tolist()
        with tempfile.TemporaryDirectory() as tmpdir:
            out_csv = Path(tmpdir) / "mhcflurry_out.csv"
            if run_mhcflurry(peptides, alleles, out_csv):
                mhc_df = pd.read_csv(out_csv)
                # Pick best allele per peptide
                score_candidates = [
                    'mhcflurry_presentation_percentile',
                    'presentation_percentile',
                    'mhcflurry_affinity_percentile',
                    'affinity_percentile',
                ]
                score_col = next((col for col in score_candidates if col in mhc_df.columns), None)
                if score_col is None:
                    print(f"Warning: MHCflurry output missing supported score columns: {score_candidates}")
                    group_copy = group.copy()
                    group_copy['best_hla'] = ""
                    group_copy['binding_rank'] = None
                    results.append(group_copy)
                    continue
                best_mhc = mhc_df.sort_values(score_col).groupby("peptide").first().reset_index()
                merged = group.merge(best_mhc[['peptide', 'allele', score_col]], on='peptide', how='left')
                merged.rename(columns={'allele': 'best_hla', score_col: 'binding_rank'}, inplace=True)
                results.append(merged)
            else:
                group_copy = group.copy()
                group_copy['best_hla'] = ""
                group_copy['binding_rank'] = None
                results.append(group_copy)
                
    if not results:
        print("No results to rank.")
        return
        
    df_ranked = pd.concat(results)
    
    # 3. Expression (linking by source_protein/gene if available)
    if 'expression_tpm' not in df_ranked.columns:
        df_ranked['expression_tpm'] = 0.0 # Default
        
    if 'rna_expr_path' in manifest.columns:
        # Build a manifest lookup: patient_id → rna_expr_path (one row per patient)
        rna_cols = ['patient_id', 'rna_expr_path']
        if 'rna_source' not in manifest.columns:
            manifest['rna_source'] = 'missing'
        rna_cols.append('rna_source')
        patient_rna = manifest[rna_cols].drop_duplicates('patient_id').dropna(subset=['rna_expr_path'])
        
        # Pre-load all unique RNA-seq files
        rna_dfs = {}
        for _, meta_row in patient_rna.iterrows():
            rna_path = meta_row['rna_expr_path']
            if rna_path and Path(rna_path).exists():
                try:
                    rdf = pd.read_csv(rna_path, sep='\t')
                    # Normalise: ensure we have gene and expression_tpm columns
                    if 'gene' not in rdf.columns and 'protein_id' in rdf.columns:
                        rdf = rdf.rename(columns={'protein_id': 'gene'})
                    rdf = rdf[['gene', 'expression_tpm']].copy()
                    rdf['gene'] = rdf['gene'].map(normalize_protein_id)
                    rdf['expression_tpm'] = pd.to_numeric(rdf['expression_tpm'], errors='coerce').fillna(0.0)
                    rna_dfs[rna_path] = rdf.groupby('gene', as_index=False)['expression_tpm'].max()
                except Exception as e:
                    print(f"Failed to load RNA-seq data from {rna_path}: {e}")
        
        if rna_dfs:
            # Tag each candidate row with the rna_expr_path for its patient
            df_ranked = df_ranked.merge(
                patient_rna,
                left_on='sample_id', right_on='patient_id', how='left'
            )
            source_col = 'source_protein' if 'source_protein' in df_ranked.columns else None
            if source_col:
                expression_lookup = {}
                for path, rdf in rna_dfs.items():
                    for _, expr_row in rdf.iterrows():
                        gene = str(expr_row['gene'])
                        if gene:
                            expression_lookup[(str(path), gene)] = float(expr_row['expression_tpm'])

                df_ranked['expression_tpm'] = df_ranked.apply(
                    lambda row: expression_for_source(row[source_col], row.get('rna_expr_path', ''), expression_lookup),
                    axis=1
                )
                df_ranked.rename(columns={'patient_id_x': 'patient_id'}, inplace=True)
    if 'rna_source' not in df_ranked.columns:
        if {'sample_id', 'patient_id', 'rna_source'}.issubset(manifest.columns):
            patient_sources = manifest[['patient_id', 'rna_source']].drop_duplicates('patient_id')
            df_ranked = df_ranked.merge(patient_sources, left_on='sample_id', right_on='patient_id', how='left')
            df_ranked.rename(columns={'patient_id_x': 'patient_id'}, inplace=True)
        else:
            df_ranked['rna_source'] = 'missing'
    df_ranked['rna_source'] = df_ranked['rna_source'].fillna('missing').map(normalize_rna_source)
    df_ranked['expression_evidence_status'] = df_ranked['rna_source'].map(expression_evidence_status)
    df_ranked['expression_supports_biology'] = df_ranked['rna_source'].map(expression_can_support_biology)

    if not args.keep_duplicate_rows:
        before = len(df_ranked)
        df_ranked = collapse_duplicate_candidates(df_ranked)
        if len(df_ranked) != before:
            print(f"Collapsed duplicate sample+peptide rows: {before} -> {len(df_ranked)}")

    # 4. Evidence Class
    def assign_class(row):
        has_binding = pd.notnull(row['binding_rank']) and row['binding_rank'] <= args.binding_rank_cutoff
        has_expression = row['expression_tpm'] >= args.tpm_cutoff and bool(row.get('expression_supports_biology', False))
        
        if row.get('mutation_type') == 'missense' and has_binding and has_expression:
            return 'A'
        elif has_binding and has_expression:
            return 'B'
        else:
            return 'C'
            
    df_ranked['evidence_class'] = df_ranked.apply(assign_class, axis=1)
    df_ranked['evidence_limitations'] = df_ranked.apply(
        lambda row: "" if row['expression_supports_biology'] else f"RNA source `{row['rna_source']}` cannot support biological expression evidence",
        axis=1,
    )
    
    # Sort by class and then by rank/score
    class_order = {'A': 0, 'B': 1, 'C': 2}
    df_ranked['class_sort'] = df_ranked['evidence_class'].map(class_order)
    df_ranked.sort_values(['class_sort', 'binding_rank'], ascending=[True, True], inplace=True)
    df_ranked.drop(columns=['class_sort'], inplace=True)
    
    df_ranked.to_csv(args.output, sep="\t", index=False)
    print(f"Saved {len(df_ranked)} ranked candidates to {args.output}")

if __name__ == "__main__":
    main()
