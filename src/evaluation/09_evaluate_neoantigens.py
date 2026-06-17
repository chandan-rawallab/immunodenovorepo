#!/usr/bin/env python3
"""Activity 4: Evaluate de novo neoantigens against the Bassani-Sternberg 2016 validated list.

Usage:
    python 09_evaluate_neoantigens.py \
        --input results/ranked_neoantigens.tsv \
        --validated data/reference/s2_dataset_extracted/Dataset1/Dataset1.txt \
        --manifest configs/sample_manifest.tsv \
        --output results/evaluation_report.md \
        --top_n 10 25 50
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path

# Patients in PXD005231 — map from manifest patient_id/validation_id to Dataset1 intensity column names
# Dataset1.txt columns: "Intensity Apher1", "Intensity Apher6", "Intensity CD165",
#   "Intensity CM467", "Intensity GD149", "Intensity MD155", "Intensity PD42",
#   "Intensity RA957", "Intensity TIL1", "Intensity TIL3"
PATIENT_TO_COLUMN = {
    "Apher-1":    "Intensity Apher1",
    "Apher1":     "Intensity Apher1",
    "Apher-6":    "Intensity Apher6",
    "Apher6":     "Intensity Apher6",
    "CD165":      "Intensity CD165",
    "CM467":      "Intensity CM467",
    "GD149":      "Intensity GD149",
    "MD155":      "Intensity MD155",
    "PD42":       "Intensity PD42",
    "RA957":      "Intensity RA957",
    "TIL1":       "Intensity TIL1",
    "TIL3":       "Intensity TIL3",
    "pooledTIL3": "Intensity TIL3",
}

def load_validated_peptides(validated_path: Path) -> dict[str, set]:
    """
    Parse Dataset1.txt to produce a per-patient set of validated peptides.
    A peptide is considered 'detected' for a patient if the corresponding
    Intensity column is non-NaN and > 0.
    """
    print(f"Loading validated peptides from {validated_path}...")
    df = pd.read_csv(validated_path, sep="\t", low_memory=False)

    intensity_cols = [c for c in df.columns if c.startswith("Intensity")]
    validated: dict[str, set] = {}

    for col in intensity_cols:
        patient_key = col.replace("Intensity ", "").strip()
        mask = df[col].notna() & (df[col] > 0)
        validated[patient_key] = set(df.loc[mask, "Sequence"].dropna().str.strip())
        print(f"  {patient_key}: {len(validated[patient_key])} validated peptides")

    return validated


def precision_recall_at_n(
    predicted: list[str],
    validated: set[str],
    top_n: int
) -> dict:
    """Compute precision and recall for the top-N predicted peptides."""
    top = predicted[:top_n]
    tp = sum(1 for p in top if p in validated)
    precision = tp / max(len(top), 1)
    recall = tp / max(len(validated), 1)
    return {"top_n": top_n, "tp": tp, "precision": precision, "recall": recall}


def evaluate_per_patient(
    ranked: pd.DataFrame,
    validated: dict[str, set],
    manifest: pd.DataFrame,
    top_ns: list[int],
) -> list[dict]:
    rows = []
    # Get unique patients from the ranked output
    if "sample_id" not in ranked.columns:
        print("WARNING: 'sample_id' column missing in ranked output. Cannot evaluate per patient.")
        return rows

    validation_lookup = {}
    if {"patient_id", "validation_id"}.issubset(manifest.columns):
        validation_lookup = dict(zip(manifest["patient_id"].astype(str), manifest["validation_id"].astype(str)))

    for patient_id in ranked["sample_id"].unique():
        validation_id = validation_lookup.get(str(patient_id), str(patient_id))
        # Map validation_id → Dataset1 column key
        val_key = None
        for k, col in PATIENT_TO_COLUMN.items():
            if k == validation_id:
                val_key = col.replace("Intensity ", "").strip()
                break

        if val_key is None or val_key not in validated:
            print(f"  Skipping {patient_id} (validation_id={validation_id}): no matching validated set.")
            continue

        val_set = validated[val_key]
        # Peptides for this patient, ranked by evidence class then binding_rank
        patient_df = ranked[ranked["sample_id"] == patient_id].copy()
        peptide_list = patient_df["peptide"].tolist()

        for top_n in top_ns:
            metrics = precision_recall_at_n(peptide_list, val_set, top_n)
            metrics["patient_id"] = patient_id
            metrics["validation_id"] = validation_id
            metrics["n_predicted"] = len(peptide_list)
            metrics["n_validated"] = len(val_set)
            rows.append(metrics)

    return rows


def write_report(
    metrics_rows: list[dict],
    ranked: pd.DataFrame,
    output_path: Path,
    top_ns: list[int],
) -> None:
    lines = []
    lines.append("# Neoepitope Evaluation Report\n")
    lines.append(f"_Evaluated against Bassani-Sternberg 2016 Dataset1 (immunopeptidome MS)_\n\n")
    lines.append("> Evaluation is independent only when the ranked input was not extracted from Dataset1 helper/demo scripts.\n\n")
    if "rna_source" in ranked.columns and (ranked["rna_source"] == "mock_debug").any():
        lines.append("> RNA source warning: at least one candidate uses mock_debug RNA; TPM-supported evidence classes are pipeline-debug only.\n\n")

    # Overall statistics
    if not ranked.empty:
        lines.append("## Pipeline Output Summary\n")
        lines.append(f"- **Total ranked candidates:** {len(ranked)}\n")
        if "evidence_class" in ranked.columns:
            counts = ranked["evidence_class"].value_counts()
            for cls in ["A", "B", "C"]:
                lines.append(f"  - Class {cls}: {counts.get(cls, 0)}\n")
        lines.append("\n")

    if not metrics_rows:
        lines.append("_No per-patient metrics could be computed. Ensure sample_id values match validated patient IDs._\n")
    else:
        df_metrics = pd.DataFrame(metrics_rows)

        for top_n in top_ns:
            subset = df_metrics[df_metrics["top_n"] == top_n].copy()
            if subset.empty:
                continue
            lines.append(f"## Precision & Recall @ Top-{top_n}\n\n")
            lines.append("| Patient | Validation ID | Predicted | Validated | TP | Precision | Recall |\n")
            lines.append("|:--------|:--------------|----------:|----------:|---:|----------:|-------:|\n")
            for _, r in subset.sort_values("patient_id").iterrows():
                lines.append(
                    f"| {r['patient_id']} | {r.get('validation_id', r['patient_id'])} | {r['n_predicted']} | {r['n_validated']} "
                    f"| {r['tp']} | {r['precision']:.3f} | {r['recall']:.4f} |\n"
                )
            avg_prec = subset["precision"].mean()
            avg_rec  = subset["recall"].mean()
            lines.append(f"\n**Mean Precision@{top_n}:** {avg_prec:.3f} | **Mean Recall@{top_n}:** {avg_rec:.4f}\n\n")

    # Top-20 candidates table
    lines.append("## Top-20 Ranked Candidates (All Patients)\n\n")
    cols = [c for c in ["peptide", "sample_id", "evidence_class", "binding_rank", "expression_tpm",
                         "mutation_type", "mutation_pos", "wt_aa", "mut_aa", "source_protein"] if c in ranked.columns]
    lines.append(ranked[cols].head(20).to_markdown(index=False))
    lines.append("\n")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("".join(lines))
    print(f"\nEvaluation report saved to {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Evaluate predicted neoantigens against validated list.")
    parser.add_argument("--input",     type=Path, required=True,  help="Path to ranked_neoantigens.tsv")
    parser.add_argument("--validated", type=Path, required=True,  help="Path to Dataset1.txt (Bassani-Sternberg 2016)")
    parser.add_argument("--manifest",  type=Path, required=True,  help="Path to sample_manifest.tsv")
    parser.add_argument("--output",    type=Path, required=True,  help="Path for evaluation_report.md")
    parser.add_argument("--top_n",     type=int,  nargs="+",      default=[10, 25, 50],
                        help="Top-N cutoffs for precision/recall (default: 10 25 50)")
    args = parser.parse_args()

    # Load ranked candidates
    print(f"Loading ranked candidates from {args.input}...")
    ranked = pd.read_csv(args.input, sep="\t")
    print(f"  {len(ranked)} candidates loaded.")

    # Load validated peptides from Dataset1
    validated = load_validated_peptides(args.validated)

    # Load manifest
    manifest = pd.read_csv(args.manifest, sep="\t")

    # Compute per-patient metrics
    print("\nComputing per-patient precision/recall...")
    metrics_rows = evaluate_per_patient(ranked, validated, manifest, args.top_n)

    # Write report
    write_report(metrics_rows, ranked, args.output, args.top_n)


if __name__ == "__main__":
    main()
