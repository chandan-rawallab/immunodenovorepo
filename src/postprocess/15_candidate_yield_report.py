#!/usr/bin/env python3
"""
Candidate Yield Report (Phase 3, Priority 4)

Produces the mandatory waterfall/funnel table showing candidate attrition
from 1.44 million raw unlabeled spectra down to the final ranked neoantigen
candidates. This is a key documentation figure for the pipeline.

Usage:
    .venv/bin/python src/postprocess/15_candidate_yield_report.py \\
        --candidates    results/de_novo_candidates.tsv \\
        --filtered      results/filtered_neoantigens.tsv \\
        --ranked        results/ranked_neoantigens.tsv \\
        --mgf-unlabeled data/mgf_unlabeled \\
        --manifest      configs/sample_manifest.tsv \\
        --output        results/candidate_yield_report.md \\
        --output-json   results/candidate_yield_report.json
"""

import argparse
import json
from pathlib import Path
import pandas as pd


# ---------------------------------------------------------------------------
# Spectrum counting
# ---------------------------------------------------------------------------

def count_spectra_in_dir(mgf_dir: Path) -> dict[str, int]:
    """Count BEGIN IONS blocks per MGF file, return {run_id: count}."""
    counts: dict[str, int] = {}
    for mgf_file in sorted(mgf_dir.glob("*.mgf")):
        run_id = mgf_file.stem.replace("_unlabeled", "")
        n = 0
        with open(mgf_file, "rb") as fh:
            for line in fh:
                if line.strip() == b"BEGIN IONS":
                    n += 1
        counts[run_id] = n
    return counts


# ---------------------------------------------------------------------------
# Per-patient breakdown
# ---------------------------------------------------------------------------

def per_patient_breakdown(
    candidates_df: pd.DataFrame,
    filtered_df: pd.DataFrame,
    ranked_df: pd.DataFrame,
    spectrum_counts: dict[str, int],
    manifest: pd.DataFrame,
) -> pd.DataFrame:
    """Produce per-patient contribution table."""
    # Map run_id to patient_id from manifest
    run_to_patient = dict(
        zip(manifest["run_id"].str.strip(), manifest["patient_id"].str.strip())
    )
    
    rows = []
    patients = sorted(manifest["patient_id"].str.strip().unique())
    for patient in patients:
        runs = [r for r, p in run_to_patient.items() if p == patient]
        
        unlabeled = sum(spectrum_counts.get(r, 0) for r in runs)
        
        if "run_id" in candidates_df.columns:
            cand_count = candidates_df[candidates_df["run_id"].isin(runs)].shape[0]
        else:
            cand_count = 0
        
        if "sample_id" in filtered_df.columns:
            filt_count = filtered_df[filtered_df["sample_id"] == patient].shape[0]
            missense_count = filtered_df[
                (filtered_df["sample_id"] == patient) &
                (filtered_df.get("mutation_type", pd.Series(dtype=str)) == "missense")
            ].shape[0]
        else:
            filt_count = 0
            missense_count = 0

        if "sample_id" in ranked_df.columns:
            ranked_count = ranked_df[ranked_df["sample_id"] == patient].shape[0]
            class_a = ranked_df[
                (ranked_df["sample_id"] == patient) &
                (ranked_df.get("evidence_class", pd.Series(dtype=str)) == "A")
            ].shape[0]
        else:
            ranked_count = 0
            class_a = 0

        rows.append({
            "patient": patient,
            "runs": len(runs),
            "unlabeled_spectra": unlabeled,
            "denovo_predictions": cand_count,
            "after_filters": filt_count,
            "missense_only": missense_count,
            "ranked": ranked_count,
            "class_a": class_a,
        })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Markdown report builder
# ---------------------------------------------------------------------------

def build_report(
    spectrum_counts: dict[str, int],
    candidates_df: pd.DataFrame,
    filtered_df: pd.DataFrame,
    ranked_df: pd.DataFrame,
    manifest: pd.DataFrame,
) -> tuple[str, dict]:
    """Build markdown report and metrics dict."""

    total_unlabeled = sum(spectrum_counts.values())
    total_runs = len(spectrum_counts)

    n_denovo = len(candidates_df)

    # Score filter stage
    score_col = "score"
    if score_col in candidates_df.columns:
        n_score_pass = (candidates_df[score_col] >= -0.5).sum()
        n_score_pass_strict = (candidates_df[score_col] >= 0.0).sum()
    else:
        n_score_pass = n_denovo
        n_score_pass_strict = n_denovo

    # Filtered breakdown
    n_filtered = len(filtered_df)
    if "mutation_type" in filtered_df.columns:
        n_missense = (filtered_df["mutation_type"] == "missense").sum()
        n_other = (filtered_df["mutation_type"] == "other").sum()
    else:
        n_missense = 0
        n_other = n_filtered

    # Ranked breakdown
    n_ranked = len(ranked_df)
    if "evidence_class" in ranked_df.columns:
        n_class_a = (ranked_df["evidence_class"] == "A").sum()
        n_class_b = (ranked_df["evidence_class"] == "B").sum()
        n_class_c = (ranked_df["evidence_class"] == "C").sum()
    else:
        n_class_a = n_class_b = n_class_c = 0

    # Retention rates
    def pct(n, total):
        if total == 0:
            return "N/A"
        return f"{n/total*100:.3f}%"

    # Per-patient table
    patient_df = per_patient_breakdown(
        candidates_df, filtered_df, ranked_df, spectrum_counts, manifest
    )

    metrics = {
        "total_runs": total_runs,
        "total_unlabeled_spectra": total_unlabeled,
        "denovo_predictions": n_denovo,
        "score_pass_threshold_minus05": int(n_score_pass),
        "score_pass_threshold_0": int(n_score_pass_strict),
        "after_all_filters": n_filtered,
        "missense_only": int(n_missense),
        "other_type": int(n_other),
        "ranked_total": n_ranked,
        "class_a": int(n_class_a),
        "class_b": int(n_class_b),
        "class_c": int(n_class_c),
        "overall_retention_pct": float(n_ranked / total_unlabeled * 100) if total_unlabeled else 0,
    }

    # ---- Markdown ----
    md = []
    md.append("# Candidate Yield Report — Objective 3")
    md.append("")
    md.append("Pipeline: 31 patient runs → Personalized CNN-LSTM de novo sequencing → Tiered neoantigen candidates")
    md.append("")
    md.append("---")
    md.append("")
    md.append("## Candidate Attrition Waterfall")
    md.append("")
    md.append("| Stage | Count | Retention vs. Previous Stage |")
    md.append("|:------|------:|-----------------------------:|")
    md.append(f"| Unlabeled spectra (input) | {total_unlabeled:,} | — |")
    md.append(f"| De novo predictions (CNN-LSTM) | {n_denovo:,} | {pct(n_denovo, total_unlabeled)} |")
    md.append(f"| Score ≥ -0.5 | {n_score_pass:,} | {pct(n_score_pass, n_denovo)} |")
    md.append(f"| Score ≥ 0.0 (strict) | {n_score_pass_strict:,} | {pct(n_score_pass_strict, n_score_pass)} |")
    md.append(f"| After all filters (07_filter) | {n_filtered:,} | {pct(n_filtered, n_score_pass)} |")
    md.append(f"| — Missense substitution | {n_missense:,} | {pct(n_missense, n_filtered)} |")
    md.append(f"| — Other sequence novelty | {n_other:,} | {pct(n_other, n_filtered)} |")
    md.append(f"| Final ranked candidates | {n_ranked:,} | {pct(n_ranked, n_filtered)} |")
    md.append(f"| — Class A (strong binder + expressed) | {n_class_a:,} | {pct(n_class_a, n_ranked)} |")
    md.append(f"| — Class B (binder or expressed) | {n_class_b:,} | {pct(n_class_b, n_ranked)} |")
    md.append(f"| — Class C (low confidence) | {n_class_c:,} | {pct(n_class_c, n_ranked)} |")
    md.append("")
    md.append(f"> **Overall retention rate**: {pct(n_ranked, total_unlabeled)} of input spectra become ranked candidates.")
    md.append("")
    md.append("---")
    md.append("")
    md.append("## Per-Patient Breakdown")
    md.append("")
    if not patient_df.empty:
        md.append("| Patient | Runs | Unlabeled Spectra | De Novo Preds | After Filters | Missense | Ranked | Class A |")
        md.append("|:--------|-----:|------------------:|--------------:|--------------:|---------:|-------:|--------:|")
        for _, row in patient_df.iterrows():
            md.append(
                f"| {row['patient']} | {row['runs']} "
                f"| {row['unlabeled_spectra']:,} "
                f"| {row['denovo_predictions']:,} "
                f"| {row['after_filters']:,} "
                f"| {row['missense_only']:,} "
                f"| {row['ranked']:,} "
                f"| {row['class_a']:,} |"
            )
    else:
        md.append("_Per-patient data not available — check sample_id column in output files._")
    md.append("")
    md.append("---")
    md.append("")
    md.append("## Scientific Caveats")
    md.append("")
    md.append("1. **Train/Test Leakage**: The current train/test split operates at the PSM level, not the unique-peptide level.")
    md.append("   Audit shows 62.4% of test peptide sequences also appear in training. Exact match accuracy (25.1%) is inflated.")
    md.append("")
    md.append("2. **Expression values are CCLE surrogates**: No patient-matched RNA-seq is available.")
    md.append("   All `expression_tpm` values derive from tissue-matched CCLE cell lines.")
    md.append("   `rna_source = ccle_lcl_surrogate` for all 31 runs.")
    md.append("")
    md.append("3. **Mutation annotation is peptide-level**: The `predicted_protein_change` (e.g., G12V) is derived")
    md.append("   from Levenshtein-1 alignment against the reference proteome — not from a VCF/MAF somatic variant call.")
    md.append("   Interpret as a sequence-level inference, not a confirmed somatic mutation.")
    md.append("")
    md.append("4. **max(TPM) across multi-mapped proteins is biologically optimistic**: When a peptide maps to")
    md.append("   multiple UniProt accessions, expression is taken as the maximum across all matching proteins.")

    return "\n".join(md), metrics


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate the candidate yield waterfall report."
    )
    parser.add_argument(
        "--candidates",
        type=Path,
        default=Path("results/de_novo_candidates.tsv"),
    )
    parser.add_argument(
        "--filtered",
        type=Path,
        default=Path("results/filtered_neoantigens.tsv"),
    )
    parser.add_argument(
        "--ranked",
        type=Path,
        default=Path("results/ranked_neoantigens.tsv"),
    )
    parser.add_argument(
        "--mgf-unlabeled",
        type=Path,
        default=Path("data/mgf_unlabeled"),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("configs/sample_manifest.tsv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/candidate_yield_report.md"),
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path("results/candidate_yield_report.json"),
    )
    args = parser.parse_args()

    print("Counting unlabeled spectra per run (this takes ~60 seconds)...")
    spectrum_counts = count_spectra_in_dir(args.mgf_unlabeled)
    print(f"  Total: {sum(spectrum_counts.values()):,} spectra across {len(spectrum_counts)} runs")

    print("Loading pipeline output files...")
    candidates_df = pd.read_csv(args.candidates, sep="\t", low_memory=False)
    filtered_df = pd.read_csv(args.filtered, sep="\t", low_memory=False)
    ranked_df = pd.read_csv(args.ranked, sep="\t", low_memory=False)
    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str)

    print("Building report...")
    report_md, metrics = build_report(
        spectrum_counts, candidates_df, filtered_df, ranked_df, manifest
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(report_md)
    print(f"Markdown report saved to {args.output}")

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_json, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"JSON metrics saved to {args.output_json}")

    # Print waterfall to stdout
    print("\n=== Candidate Yield Summary ===")
    for line in report_md.split("\n"):
        if line.startswith("|"):
            print(line)


if __name__ == "__main__":
    main()
