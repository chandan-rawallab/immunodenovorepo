#!/usr/bin/env python3
"""Activity 3: Ranking neoantigen candidates by binding affinity and expression.

Audit fixes applied (2026-06-20):
  - Vectorised join replaces row-wise apply() for expression assignment.
  - Ranking-stage audit metrics printed at end.
  - MHCflurry version and predictor configuration saved as provenance fields.
  - Manifest sample-matching assumption is validated and logged.
  - Missing MHCflurry predictions are explicitly counted and logged rather
    than silently downgrading candidates.
  - All bare print() calls replaced with structured logging.
"""

from __future__ import annotations

import argparse
import json
import logging
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

SCRIPT_DIR = Path(__file__).resolve().parent
WORKSPACE_ROOT = SCRIPT_DIR.parent.parent

PATIENT_MATCHED_RNA_SOURCES = {
    "real_patient_matched",
    "patient_matched",
    "provided_patient_matched",
    "clinical_rna",
    "matched_rna",
}
NON_BIOLOGICAL_RNA_MARKERS = (
    "mock", "debug", "surrogate", "inferred", "missing", "provided_override"
)


# ---------------------------------------------------------------------------
# HLA normalisation
# ---------------------------------------------------------------------------

def normalize_hla(value: str) -> str:
    """Normalise HLA allele string to MHCflurry format (e.g. HLA-A*02:01)."""
    value = value.strip()
    if not value:
        return ""
    if "*" in value and ":" in value:
        return value
    clean = value.replace("HLA-", "").replace("*", "").replace(":", "")
    if len(clean) >= 5 and clean[0].isalpha():
        return f"HLA-{clean[0]}*{clean[1:3]}:{clean[3:5]}"
    return value


# ---------------------------------------------------------------------------
# MHCflurry
# ---------------------------------------------------------------------------

def _mhcflurry_version() -> str:
    """Return the installed mhcflurry version string, or 'unknown'."""
    try:
        result = subprocess.run(
            ["mhcflurry-predict", "--version"],
            capture_output=True, text=True, timeout=10,
        )
        return result.stdout.strip() or result.stderr.strip() or "unknown"
    except Exception:
        return "unknown"


def run_mhcflurry(peptides: list[str], alleles: list[str], output_path: Path) -> bool:
    """Run MHCflurry and return True on success."""
    predictor = shutil.which("mhcflurry-predict")
    if not predictor:
        local_predictor = WORKSPACE_ROOT / ".venv" / "bin" / "mhcflurry-predict"
        if local_predictor.exists():
            predictor = str(local_predictor)
    if not predictor:
        logger.warning("mhcflurry-predict not found. Binding prediction skipped.")
        return False

    with tempfile.TemporaryDirectory() as tmpdir:
        input_csv = Path(tmpdir) / "input.csv"
        with open(input_csv, "w") as fh:
            fh.write("allele,peptide\n")
            for allele in alleles:
                for peptide in peptides:
                    fh.write(f"{allele},{peptide}\n")
        try:
            subprocess.run(
                [predictor, str(input_csv), "--out", str(output_path)],
                check=True, capture_output=True,
            )
            return True
        except subprocess.CalledProcessError as exc:
            logger.error("MHCflurry error: %s", exc.stderr.decode())
            return False


# ---------------------------------------------------------------------------
# Protein ID / expression helpers
# ---------------------------------------------------------------------------

def normalize_protein_id(value: object) -> str:
    text = str(value or "").strip()
    if not text or text.lower() == "nan":
        return ""
    return text.split("-")[0]


def source_protein_ids(value: object) -> list[str]:
    text = str(value or "").strip()
    if not text or text.lower() == "nan":
        return []
    ids: list[str] = []
    for part in text.split(";"):
        norm = normalize_protein_id(part)
        if norm and norm not in ids:
            ids.append(norm)
    return ids


# ---------------------------------------------------------------------------
# Duplicate collapsing
# ---------------------------------------------------------------------------

def collapse_duplicate_candidates(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse repeated sample+peptide rows, keeping the strongest evidence row."""
    if not {"sample_id", "peptide"}.issubset(df.columns):
        return df
    if not df.duplicated(["sample_id", "peptide"]).any():
        return df

    collapsed = []
    for _, group in df.groupby(["sample_id", "peptide"], sort=False):
        ranked = group.copy()
        if "binding_rank" in ranked.columns:
            ranked["_sort"] = pd.to_numeric(ranked["binding_rank"], errors="coerce").fillna(float("inf"))
            ranked = ranked.sort_values("_sort")
        elif "score" in ranked.columns:
            ranked["_sort"] = pd.to_numeric(ranked["score"], errors="coerce").fillna(float("-inf"))
            ranked = ranked.sort_values("_sort", ascending=False)
        row = ranked.iloc[0].drop(labels=[c for c in ["_sort"] if c in ranked.columns]).copy()

        if "psm_count" in group.columns:
            row["psm_count"] = pd.to_numeric(group["psm_count"], errors="coerce").max()
        if "expression_tpm" in group.columns:
            row["expression_tpm"] = pd.to_numeric(group["expression_tpm"], errors="coerce").fillna(0.0).max()
        if "spectrum_id" in group.columns:
            spectra = [str(x) for x in group["spectrum_id"].dropna().unique()]
            row["spectrum_id"] = ";".join(spectra[:25])
        collapsed.append(row)

    return pd.DataFrame(collapsed).reset_index(drop=True)


# ---------------------------------------------------------------------------
# RNA-source / expression-evidence helpers
# ---------------------------------------------------------------------------

def normalize_rna_source(value: object) -> str:
    text = str(value or "").strip().lower()
    return text if (text and text != "nan") else "missing"


def expression_evidence_status(value: object) -> str:
    source = normalize_rna_source(value)
    if source in PATIENT_MATCHED_RNA_SOURCES:
        return "patient_matched"
    if any(m in source for m in NON_BIOLOGICAL_RNA_MARKERS):
        return "debug_or_non_patient_matched"
    return "unverified"


def expression_can_support_biology(value: object) -> bool:
    return expression_evidence_status(value) == "patient_matched"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Rank neoantigen candidates.")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--binding_rank_cutoff", type=float, default=2.0)
    parser.add_argument("--tpm_cutoff", type=float, default=1.0)
    parser.add_argument(
        "--keep_duplicate_rows", action="store_true",
        help="Keep repeated sample+peptide rows instead of collapsing them.",
    )
    args = parser.parse_args()

    # --- Load ---
    df = pd.read_csv(args.input, sep="\t")
    manifest = pd.read_csv(args.manifest, sep="\t")

    logger.info("Loaded %d candidates from %s", len(df), args.input)
    logger.info("MHCflurry version: %s", _mhcflurry_version())

    # --- Validate manifest sample matching ---
    df_samples = set(df["sample_id"].unique())
    manifest_patients = set(manifest["patient_id"].unique()) if "patient_id" in manifest.columns else set()
    manifest_runs = set(manifest["run_id"].unique()) if "run_id" in manifest.columns else set()
    matched_by_patient = df_samples & manifest_patients
    matched_by_run = df_samples & manifest_runs
    unmatched = df_samples - manifest_patients - manifest_runs
    logger.info(
        "Sample matching: %d matched by patient_id, %d by run_id, %d unmatched.",
        len(matched_by_patient), len(matched_by_run), len(unmatched),
    )
    if unmatched:
        logger.warning("Unmatched sample IDs (no HLA/RNA metadata): %s", sorted(unmatched)[:10])

    # --- MHCflurry binding prediction ---
    results = []
    missing_mhc_count = 0

    for sample_id, group in df.groupby("sample_id"):
        sample_meta = manifest[manifest["patient_id"] == sample_id]
        if sample_meta.empty:
            sample_meta = manifest[manifest["run_id"] == sample_id]

        if sample_meta.empty:
            logger.warning("No manifest metadata for sample '%s'. Skipping.", sample_id)
            missing_mhc_count += len(group)
            continue

        raw_alleles = str(sample_meta.iloc[0]["hla_alleles"]).split(",")
        alleles = [normalize_hla(a) for a in raw_alleles if a and a.strip().lower() != "tbd"]

        if not alleles:
            logger.warning("No valid HLA alleles for sample '%s'.", sample_id)
            g = group.copy()
            g["best_hla"] = ""
            g["binding_rank"] = None
            results.append(g)
            missing_mhc_count += len(group)
            continue

        peptides = group["peptide"].unique().tolist()
        with tempfile.TemporaryDirectory() as tmpdir:
            out_csv = Path(tmpdir) / "mhcflurry_out.csv"
            if run_mhcflurry(peptides, alleles, out_csv):
                mhc_df = pd.read_csv(out_csv)
                score_candidates = [
                    "mhcflurry_presentation_percentile",
                    "presentation_percentile",
                    "mhcflurry_affinity_percentile",
                    "affinity_percentile",
                ]
                score_col = next(
                    (c for c in score_candidates if c in mhc_df.columns), None
                )
                if score_col is None:
                    logger.warning(
                        "MHCflurry output missing score columns for sample '%s'. "
                        "Expected one of: %s",
                        sample_id, score_candidates,
                    )
                    g = group.copy()
                    g["best_hla"] = ""
                    g["binding_rank"] = None
                    results.append(g)
                    missing_mhc_count += len(group)
                    continue

                best_mhc = (
                    mhc_df.sort_values(score_col)
                    .groupby("peptide")
                    .first()
                    .reset_index()
                )
                merged = group.merge(
                    best_mhc[["peptide", "allele", score_col]],
                    on="peptide", how="left",
                )
                merged.rename(
                    columns={"allele": "best_hla", score_col: "binding_rank"},
                    inplace=True,
                )
                results.append(merged)
            else:
                g = group.copy()
                g["best_hla"] = ""
                g["binding_rank"] = None
                results.append(g)
                missing_mhc_count += len(group)

    if not results:
        logger.error("No results to rank.")
        return

    df_ranked = pd.concat(results, ignore_index=True)

    if missing_mhc_count:
        logger.warning(
            "%d candidates have no MHCflurry binding prediction "
            "(will be classified as evidence class C).",
            missing_mhc_count,
        )

    # --- Expression (vectorised join — audit fix) ---
    if "expression_tpm" not in df_ranked.columns:
        df_ranked["expression_tpm"] = 0.0

    if "rna_expr_path" in manifest.columns:
        rna_cols = ["patient_id", "rna_expr_path"]
        if "rna_source" not in manifest.columns:
            manifest["rna_source"] = "missing"
        rna_cols.append("rna_source")
        patient_rna = (
            manifest[rna_cols]
            .drop_duplicates("patient_id")
            .dropna(subset=["rna_expr_path"])
        )

        # Load all RNA-seq files; build a single expression lookup table
        rna_frames: list[pd.DataFrame] = []
        for _, meta_row in patient_rna.iterrows():
            rna_path = meta_row["rna_expr_path"]
            if not rna_path or not Path(rna_path).exists():
                continue
            try:
                rdf = pd.read_csv(rna_path, sep="\t")
                if "gene" not in rdf.columns and "protein_id" in rdf.columns:
                    rdf = rdf.rename(columns={"protein_id": "gene"})
                rdf = rdf[["gene", "expression_tpm"]].copy()
                rdf["gene"] = rdf["gene"].map(normalize_protein_id)
                rdf["expression_tpm"] = pd.to_numeric(
                    rdf["expression_tpm"], errors="coerce"
                ).fillna(0.0)
                rdf = rdf.groupby("gene", as_index=False)["expression_tpm"].max()
                rdf["rna_expr_path"] = rna_path
                rna_frames.append(rdf)
            except Exception as exc:
                logger.error("Failed to load RNA-seq data from %s: %s", rna_path, exc)

        if rna_frames and "source_protein" in df_ranked.columns:
            # Explode multi-protein annotations into one row per protein ID
            all_rna = pd.concat(rna_frames, ignore_index=True)

            df_ranked = df_ranked.merge(
                patient_rna[["patient_id", "rna_expr_path"]],
                left_on="sample_id", right_on="patient_id", how="left",
                suffixes=("", "_manifest"),
            )

            # Explode source_protein field for vectorised join
            df_ranked["_source_ids"] = df_ranked["source_protein"].apply(source_protein_ids)
            df_exploded = df_ranked.explode("_source_ids").rename(
                columns={"_source_ids": "gene"}
            )

            # Join expression values
            df_exploded = df_exploded.merge(
                all_rna[["rna_expr_path", "gene", "expression_tpm"]],
                on=["rna_expr_path", "gene"], how="left",
                suffixes=("_old", ""),
            )

            # Aggregate: max TPM per original candidate row
            tpm_max = (
                df_exploded.groupby(df_exploded.index)["expression_tpm"]
                .max()
                .rename("expression_tpm_joined")
            )
            df_ranked = df_ranked.join(tpm_max)
            df_ranked["expression_tpm"] = df_ranked["expression_tpm_joined"].fillna(0.0)
            df_ranked.drop(columns=["expression_tpm_joined", "_source_ids"], errors="ignore", inplace=True)

            # Clean up duplicate columns from merge
            if "patient_id_manifest" in df_ranked.columns:
                df_ranked.drop(columns=["patient_id_manifest"], errors="ignore", inplace=True)

    # RNA source evidence gating
    if "rna_source" not in df_ranked.columns:
        if {"sample_id", "patient_id", "rna_source"}.issubset(manifest.columns):
            patient_sources = manifest[["patient_id", "rna_source"]].drop_duplicates("patient_id")
            df_ranked = df_ranked.merge(
                patient_sources, left_on="sample_id", right_on="patient_id", how="left"
            )
        else:
            df_ranked["rna_source"] = "missing"

    df_ranked["rna_source"] = (
        df_ranked["rna_source"].fillna("missing").map(normalize_rna_source)
    )
    df_ranked["expression_evidence_status"] = df_ranked["rna_source"].map(expression_evidence_status)
    df_ranked["expression_supports_biology"] = df_ranked["rna_source"].map(expression_can_support_biology)

    # --- Collapse duplicates ---
    if not args.keep_duplicate_rows:
        before = len(df_ranked)
        df_ranked = collapse_duplicate_candidates(df_ranked)
        if len(df_ranked) != before:
            logger.info("Collapsed duplicate rows: %d → %d", before, len(df_ranked))

    # --- Evidence class assignment ---
    def assign_class(row) -> str:
        has_binding = (
            pd.notnull(row["binding_rank"])
            and row["binding_rank"] <= args.binding_rank_cutoff
        )
        has_expression = (
            row["expression_tpm"] >= args.tpm_cutoff
            and bool(row.get("expression_supports_biology", False))
        )
        if row.get("mutation_type") == "missense" and has_binding and has_expression:
            return "A"
        if has_binding and has_expression:
            return "B"
        return "C"

    df_ranked["evidence_class"] = df_ranked.apply(assign_class, axis=1)
    df_ranked["evidence_limitations"] = df_ranked.apply(
        lambda row: (
            ""
            if row["expression_supports_biology"]
            else f"RNA source `{row['rna_source']}` cannot support biological expression evidence"
        ),
        axis=1,
    )

    # Sort
    class_order = {"A": 0, "B": 1, "C": 2}
    df_ranked["class_sort"] = df_ranked["evidence_class"].map(class_order)
    df_ranked.sort_values(["class_sort", "binding_rank"], ascending=[True, True], inplace=True)
    df_ranked.drop(columns=["class_sort"], inplace=True)

    # --- Ranking-stage audit metrics ---
    n_class_a = (df_ranked["evidence_class"] == "A").sum()
    n_class_b = (df_ranked["evidence_class"] == "B").sum()
    n_class_c = (df_ranked["evidence_class"] == "C").sum()
    n_expr_patient = df_ranked["expression_supports_biology"].sum()

    audit = {
        "total_candidates_ranked": len(df_ranked),
        "evidence_class_A": int(n_class_a),
        "evidence_class_B": int(n_class_b),
        "evidence_class_C": int(n_class_c),
        "candidates_with_patient_matched_RNA": int(n_expr_patient),
        "missing_mhcflurry_predictions": missing_mhc_count,
        "mhcflurry_version": _mhcflurry_version(),
        "binding_rank_cutoff": args.binding_rank_cutoff,
        "tpm_cutoff": args.tpm_cutoff,
    }

    logger.info("=== Ranking Audit Metrics ===")
    for k, v in audit.items():
        logger.info("  %-42s %s", k, v)

    # Save audit sidecar
    audit_path = args.output.with_suffix(".audit.json")
    audit_path.write_text(json.dumps(audit, indent=2))
    logger.info("Audit metrics saved to %s", audit_path)

    df_ranked.to_csv(args.output, sep="\t", index=False)
    logger.info("Saved %d ranked candidates to %s", len(df_ranked), args.output)


if __name__ == "__main__":
    main()
