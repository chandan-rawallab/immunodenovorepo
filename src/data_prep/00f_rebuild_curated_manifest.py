#!/usr/bin/env python3
"""Rebuild the active Objective 3 manifest from the curated publication cohort."""

import argparse
from pathlib import Path

import pandas as pd


MANIFEST_COLUMNS = [
    "study_id",
    "run_id",
    "patient_id",
    "validation_id",
    "filename",
    "cohort",
    "sample_role",
    "hla_alleles",
    "hla_source",
    "rna_expr_path",
    "rna_source",
    "raw_source",
    "psm_source",
    "include_in_pipeline",
    "notes",
]

EXCLUDED_REASON = "unverified_patient_hla_mapping"


def run_id_from_psm(path: Path) -> str:
    return path.name.replace("msms_", "").replace(".raw.txt", "").replace(".txt", "")


def discover_psm_runs(psm_dir: Path) -> dict[str, str]:
    return {run_id_from_psm(path): path.name for path in sorted(psm_dir.glob("msms_*.txt"))}


def normalize_patient_id(patient_id: str, run_id: str) -> tuple[str, str]:
    if patient_id == "CM647" or "CM647" in run_id:
        return "CM467", "CM467"
    return patient_id, patient_id


def sample_role_for(patient_id: str) -> str:
    if patient_id in {"Apher-1", "Apher-6"}:
        return "apheresis_control"
    if patient_id == "pooledTIL3":
        return "pooled_til"
    return "hla_peptidome"


def main() -> int:
    parser = argparse.ArgumentParser(description="Rebuild curated 31-run manifest and excluded-run audit table.")
    parser.add_argument("--base", type=Path, default=Path("configs/sample_manifest.tsv.bak"))
    parser.add_argument("--psm-dir", type=Path, default=Path("data/psms"))
    parser.add_argument("--output", type=Path, default=Path("configs/sample_manifest.tsv"))
    parser.add_argument("--excluded-output", type=Path, default=Path("configs/excluded_runs.tsv"))
    parser.add_argument("--cohort", default="PXD005231")
    args = parser.parse_args()

    if not args.base.exists():
        print(f"ERROR: curated base manifest not found: {args.base}")
        return 1

    base = pd.read_csv(args.base, sep="\t")
    required = {"run_id", "patient_id", "filename", "hla_alleles", "cohort"}
    missing = required - set(base.columns)
    if missing:
        print(f"ERROR: base manifest missing columns: {', '.join(sorted(missing))}")
        return 1

    psm_runs = discover_psm_runs(args.psm_dir)
    active_run_ids = set(base["run_id"].astype(str))
    excluded_run_ids = sorted(set(psm_runs) - active_run_ids)

    rows = []
    for _, row in base.iterrows():
        run_id = str(row["run_id"])
        patient_id, validation_id = normalize_patient_id(str(row["patient_id"]), run_id)
        notes = []
        if str(row["patient_id"]) != patient_id:
            notes.append(f"patient_id normalized from {row['patient_id']}")
        if run_id not in psm_runs:
            notes.append("WARNING: matching PSM file not currently present")

        rows.append(
            {
                "study_id": args.cohort,
                "run_id": run_id,
                "patient_id": patient_id,
                "validation_id": validation_id,
                "filename": psm_runs.get(run_id, str(row["filename"])),
                "cohort": str(row.get("cohort", args.cohort) or args.cohort),
                "sample_role": sample_role_for(patient_id),
                "hla_alleles": str(row["hla_alleles"]),
                "hla_source": "curated_publication",
                "rna_expr_path": "",
                "rna_source": "mock_debug",
                "raw_source": "pride_raw",
                "psm_source": "local_maxquant_psm",
                "include_in_pipeline": True,
                "notes": "; ".join(notes),
            }
        )

    manifest = pd.DataFrame(rows, columns=MANIFEST_COLUMNS)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(args.output, sep="\t", index=False)

    excluded_rows = []
    for run_id in excluded_run_ids:
        excluded_rows.append(
            {
                "run_id": run_id,
                "filename": psm_runs[run_id],
                "cohort": args.cohort,
                "reason": EXCLUDED_REASON,
                "notes": "Present in local PSM/MGF files but absent from curated publication manifest.",
            }
        )
    excluded = pd.DataFrame(excluded_rows, columns=["run_id", "filename", "cohort", "reason", "notes"])
    args.excluded_output.parent.mkdir(parents=True, exist_ok=True)
    excluded.to_csv(args.excluded_output, sep="\t", index=False)

    print(f"Wrote active manifest: {args.output} ({len(manifest)} rows)")
    print(f"Wrote excluded runs:   {args.excluded_output} ({len(excluded)} rows)")
    print("\nActive patients:")
    print(manifest.groupby(["patient_id", "validation_id"]).size().rename("runs").reset_index().to_string(index=False))
    if len(excluded):
        print("\nExcluded run IDs:")
        print("\n".join(excluded["run_id"].tolist()))

    if len(manifest) != 31:
        print(f"ERROR: expected 31 curated active rows, observed {len(manifest)}")
        return 1
    if len(excluded) != 9:
        print(f"ERROR: expected 9 excluded local rows, observed {len(excluded)}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
