#!/usr/bin/env python3
"""Preflight checks for Objective 3 pipeline data lineage and outputs."""

import argparse
from pathlib import Path

import pandas as pd


def read_tsv(path: Path, usecols=None) -> pd.DataFrame | None:
    if not path.exists():
        return None
    return pd.read_csv(path, sep="\t", usecols=usecols, low_memory=False)


def run_id_from_psm(path: Path) -> str:
    return path.name.replace("msms_", "").replace(".raw.txt", "").replace(".txt", "")


def run_id_from_mgf(path: Path) -> str:
    return path.name.replace("_unlabeled.mgf", "").replace(".mgf", "")


def discover_runs(path: Path, pattern: str, parser) -> set[str]:
    if not path.exists():
        return set()
    return {parser(p) for p in path.glob(pattern) if p.is_file()}


def fasta_stats(path: Path, sample_limit: int = 5) -> dict:
    total = 0
    human = 0
    examples = []
    if not path.exists():
        return {"exists": False, "total": 0, "human": 0, "examples": []}

    with path.open() as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            total += 1
            is_human = "OS=Homo sapiens" in line or "OX=9606" in line
            human += int(is_human)
            if not is_human and len(examples) < sample_limit:
                examples.append(line.strip())
    return {"exists": True, "total": total, "human": human, "examples": examples}


def report_set_delta(title: str, left_name: str, left: set[str], right_name: str, right: set[str]) -> list[str]:
    lines = [f"\n## {title}"]
    missing = sorted(left - right)
    extra = sorted(right - left)
    if not missing and not extra:
        lines.append(f"OK: `{left_name}` and `{right_name}` cover the same run IDs.")
        return lines

    if missing:
        lines.append(f"FAIL: {len(missing)} run IDs in `{left_name}` are missing from `{right_name}`.")
        lines.append("  " + ", ".join(missing[:20]) + (" ..." if len(missing) > 20 else ""))
    if extra:
        lines.append(f"WARN: {len(extra)} run IDs in `{right_name}` are absent from `{left_name}`.")
        lines.append("  " + ", ".join(extra[:20]) + (" ..." if len(extra) > 20 else ""))
    return lines


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate Objective 3 input/output consistency before reruns.")
    parser.add_argument("--manifest", type=Path, default=Path("configs/sample_manifest.tsv"))
    parser.add_argument("--psm-dir", type=Path, default=Path("data/psms"))
    parser.add_argument("--mgf-dir", type=Path, default=Path("data/mgf"))
    parser.add_argument("--unlabeled-mgf-dir", type=Path, default=Path("data/mgf_unlabeled"))
    parser.add_argument("--reference-fasta", type=Path, default=Path("data/reference/uniprot_human_reviewed.fasta"))
    parser.add_argument("--psms", type=Path, default=Path("results/immunopeptidome_psms.tsv"))
    parser.add_argument("--de-novo", type=Path, default=Path("results/de_novo_candidates.tsv"))
    parser.add_argument("--filtered", type=Path, default=Path("results/filtered_neoantigens.tsv"))
    parser.add_argument("--ranked", type=Path, default=Path("results/ranked_neoantigens.tsv"))
    parser.add_argument("--excluded-runs", type=Path, default=Path("configs/excluded_runs.tsv"))
    parser.add_argument("--expected-active-runs", type=int, default=None)
    parser.add_argument("--output", type=Path, default=None, help="Optional Markdown report path.")
    args = parser.parse_args()

    failures = 0
    warnings = 0
    lines = ["# Objective 3 Preflight Report"]

    manifest = read_tsv(args.manifest)
    if manifest is None:
        lines.append(f"\nFAIL: manifest not found: `{args.manifest}`")
        failures += 1
        manifest_runs = set()
    else:
        required = {"study_id", "run_id", "patient_id", "hla_alleles", "hla_source", "rna_source", "include_in_pipeline"}
        missing_cols = required - set(manifest.columns)
        if missing_cols:
            lines.append(f"\nFAIL: manifest missing columns: {', '.join(sorted(missing_cols))}")
            failures += 1
        manifest_runs = set(manifest.get("run_id", pd.Series(dtype=str)).dropna().astype(str))
        lines.append(f"\nManifest rows: {len(manifest)}; unique run IDs: {len(manifest_runs)}")
        if args.expected_active_runs is not None and len(manifest_runs) != args.expected_active_runs:
            lines.append(f"FAIL: expected {args.expected_active_runs} active manifest runs, observed {len(manifest_runs)}.")
            failures += 1
        if "rna_expr_path" not in manifest.columns:
            lines.append("WARN: manifest has no `rna_expr_path` column; expression linking will default to TPM 0.")
            warnings += 1
        else:
            missing_expr = []
            for raw_path in manifest["rna_expr_path"].dropna().astype(str):
                if raw_path and not Path(raw_path).exists():
                    missing_expr.append(raw_path)
            if missing_expr:
                lines.append(f"FAIL: {len(missing_expr)} RNA expression paths do not exist.")
                lines.append("  " + ", ".join(missing_expr[:10]))
                failures += 1
        if "validation_id" not in manifest.columns:
            lines.append("WARN: manifest has no `validation_id`; evaluation will fall back to `sample_id`.")
            warnings += 1
        if "rna_source" in manifest.columns:
            mock_rows = int((manifest["rna_source"] == "mock_debug").sum())
            if mock_rows:
                lines.append(f"WARN: {mock_rows} manifest rows use mock_debug RNA; TPM-supported evidence is debug-only.")
                warnings += 1
            surrogate_rows = int(manifest["rna_source"].astype(str).str.contains("surrogate", case=False, na=False).sum())
            if surrogate_rows:
                lines.append(f"WARN: {surrogate_rows} manifest rows use surrogate RNA; TPM support is not patient-matched biological evidence.")
                warnings += 1

    psm_file_runs = discover_runs(args.psm_dir, "msms_*.txt", run_id_from_psm)
    mgf_runs = discover_runs(args.mgf_dir, "*.mgf", run_id_from_mgf)
    unlabeled_runs = discover_runs(args.unlabeled_mgf_dir, "*_unlabeled.mgf", run_id_from_mgf)
    excluded_runs = set()
    excluded = read_tsv(args.excluded_runs)
    if excluded is not None and "run_id" in excluded.columns:
        excluded_runs = set(excluded["run_id"].dropna().astype(str))
        lines.append(f"\nExcluded run IDs: {len(excluded_runs)}")

    psm_file_runs_for_manifest = psm_file_runs - excluded_runs
    mgf_runs_for_manifest = mgf_runs - excluded_runs

    lines.extend(report_set_delta("Manifest vs PSM Files", "PSM files", psm_file_runs_for_manifest, "manifest", manifest_runs))
    lines.extend(report_set_delta("Manifest vs MGF Files", "MGF files", mgf_runs_for_manifest, "manifest", manifest_runs))
    for block in lines[-8:]:
        if block.startswith("FAIL:"):
            failures += 1
        elif block.startswith("WARN:"):
            warnings += 1

    stats = fasta_stats(args.reference_fasta)
    lines.append("\n## Reference FASTA")
    if not stats["exists"]:
        lines.append(f"FAIL: reference FASTA not found: `{args.reference_fasta}`")
        failures += 1
    else:
        lines.append(f"Entries: {stats['total']}; human-tagged entries: {stats['human']}")
        if stats["total"] and stats["human"] != stats["total"]:
            lines.append("FAIL: reference FASTA is not human-only. Non-human examples:")
            lines.extend(f"  - `{x}`" for x in stats["examples"])
            failures += 1

    psms = read_tsv(args.psms, usecols=lambda c: c in {"run_id", "sample_id", "peptide"})
    if psms is not None:
        unknown = int((psms.get("sample_id", "") == "Unknown").sum()) if "sample_id" in psms else 0
        lines.append("\n## PSM Output")
        lines.append(f"Rows: {len(psms)}; Unknown sample rows: {unknown}")
        if unknown:
            lines.append("FAIL: PSM output contains `sample_id=Unknown`; rerun Step 04 after fixing manifest.")
            failures += 1

    denovo = read_tsv(args.de_novo, usecols=lambda c: c in {"run_id", "peptide", "fdr", "score"})
    if denovo is not None:
        denovo_runs = set(denovo["run_id"].dropna().astype(str)) if "run_id" in denovo else set()
        lines.append("\n## De Novo Output")
        lines.append(f"Rows: {len(denovo)}; run IDs: {len(denovo_runs)}")
        if unlabeled_runs:
            lines.extend(report_set_delta("De Novo Output vs Current Unlabeled MGF", "de novo output", denovo_runs, "current unlabeled MGF", unlabeled_runs))
            if denovo_runs != unlabeled_runs:
                failures += 1
        if "fdr" in denovo and len(denovo):
            zero_fdr = int((denovo["fdr"] == 0).sum())
            if zero_fdr / len(denovo) > 0.25:
                lines.append(f"WARN: {zero_fdr}/{len(denovo)} de novo rows have FDR 0; target-decoy calibration needs review.")
                warnings += 1

    for label, path in [("Filtered", args.filtered), ("Ranked", args.ranked)]:
        df = read_tsv(path)
        if df is None:
            continue
        lines.append(f"\n## {label} Output")
        lines.append(f"Rows: {len(df)}")
        if "sample_id" in df:
            unknown = int((df["sample_id"] == "Unknown").sum())
            lines.append(f"Unknown sample rows: {unknown}")
            if unknown:
                lines.append(f"FAIL: {label.lower()} output contains unknown samples; upstream manifest mapping is broken.")
                failures += 1
        if {"sample_id", "peptide"}.issubset(df.columns):
            dupes = int(df.duplicated(["sample_id", "peptide"]).sum())
            lines.append(f"Duplicate sample+peptide rows: {dupes}")
            if dupes:
                lines.append("WARN: duplicate sample+peptide rows should be collapsed before evaluation.")
                warnings += 1

    lines.append("\n## Summary")
    lines.append(f"Failures: {failures}")
    lines.append(f"Warnings: {warnings}")

    text = "\n".join(lines) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(text)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
