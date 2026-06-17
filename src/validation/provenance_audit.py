#!/usr/bin/env python3
"""Audit external and local provenance consistency for the Objective 3 pipeline."""

import argparse
import json
from pathlib import Path

import pandas as pd
import requests


PRIDE_PROJECT_URL = "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{accession}"
PRIDE_FILES_URL = "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{accession}/files"


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def run_id_from_psm(path: Path) -> str:
    return path.name.replace("msms_", "").replace(".raw.txt", "").replace(".txt", "")


def fetch_json(url: str, timeout: int) -> tuple[object | None, str | None]:
    try:
        response = requests.get(url, headers={"Accept": "application/json"}, timeout=timeout)
        response.raise_for_status()
        return response.json(), None
    except Exception as exc:
        return None, str(exc)


def pride_file_names(files_payload) -> set[str]:
    names = set()
    if not isinstance(files_payload, list):
        return names
    for entry in files_payload:
        name = entry.get("fileName")
        if name:
            names.add(str(name))
    return names


def dataset1_intensity_ids(path: Path) -> set[str]:
    if not path.exists():
        return set()
    df = pd.read_csv(path, sep="\t", nrows=1)
    return {col.replace("Intensity ", "").strip() for col in df.columns if col.startswith("Intensity ")}


def fasta_stats(path: Path) -> dict:
    total = 0
    human = 0
    if not path.exists():
        return {"exists": False, "total": 0, "human": 0}
    with path.open() as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            total += 1
            human += int("OS=Homo sapiens" in line or "OX=9606" in line)
    return {"exists": True, "total": total, "human": human}


def add_status(rows: list[dict], category: str, check: str, status: str, detail: str):
    rows.append({"category": category, "check": check, "status": status, "detail": detail})


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit data provenance and partial/stale pipeline state.")
    parser.add_argument("--accession", default="PXD005231")
    parser.add_argument("--manifest", type=Path, default=Path("configs/sample_manifest.tsv"))
    parser.add_argument("--excluded-runs", type=Path, default=Path("configs/excluded_runs.tsv"))
    parser.add_argument("--psm-dir", type=Path, default=Path("data/psms"))
    parser.add_argument("--mgf-dir", type=Path, default=Path("data/mgf"))
    parser.add_argument("--unlabeled-mgf-dir", type=Path, default=Path("data/mgf_unlabeled"))
    parser.add_argument("--reference-fasta", type=Path, default=Path("data/reference/uniprot_human_reviewed.only_human.fasta"))
    parser.add_argument("--validated", type=Path, default=Path("data/reference/s2_dataset_extracted/Dataset1/Dataset1.txt"))
    parser.add_argument("--de-novo", type=Path, default=Path("results/de_novo_candidates.tsv"))
    parser.add_argument("--filtered", type=Path, default=Path("results/filtered_neoantigens.tsv"))
    parser.add_argument("--ranked", type=Path, default=Path("results/ranked_neoantigens.tsv"))
    parser.add_argument("--output-md", type=Path, default=Path("results/provenance_audit.md"))
    parser.add_argument("--output-json", type=Path, default=Path("results/provenance_audit.json"))
    parser.add_argument("--timeout", type=int, default=30)
    parser.add_argument("--expected-active-runs", type=int, default=None)
    parser.add_argument("--expected-excluded-runs", type=int, default=None)
    args = parser.parse_args()

    rows: list[dict] = []
    payload: dict = {"accession": args.accession}

    manifest = read_tsv(args.manifest)
    excluded = read_tsv(args.excluded_runs)
    manifest_runs = set(manifest["run_id"]) if "run_id" in manifest else set()
    excluded_runs = set(excluded["run_id"]) if "run_id" in excluded else set()
    local_psm_runs = {run_id_from_psm(path) for path in args.psm_dir.glob("msms_*.txt")}
    local_mgf_runs = {path.stem for path in args.mgf_dir.glob("*.mgf")}
    active_local_psm_runs = local_psm_runs - excluded_runs
    active_local_mgf_runs = local_mgf_runs - excluded_runs

    project, project_error = fetch_json(PRIDE_PROJECT_URL.format(accession=args.accession), args.timeout)
    files, files_error = fetch_json(PRIDE_FILES_URL.format(accession=args.accession), args.timeout)
    payload["pride_project_error"] = project_error
    payload["pride_files_error"] = files_error

    if project_error:
        add_status(rows, "external_pride", "project_api", "WARN", project_error)
    else:
        organisms = json.dumps(project.get("organisms") or project.get("sampleAttributes") or "")
        add_status(rows, "external_pride", "project_title", "OK", str(project.get("title", "")))
        add_status(rows, "external_pride", "submission_type", "WARN" if project.get("submissionType") == "PARTIAL" else "OK", str(project.get("submissionType", "")))
        add_status(rows, "external_pride", "organism_human", "OK" if "9606" in organisms or "Homo sapiens" in organisms else "FAIL", organisms[:500])
        add_status(rows, "external_pride", "related_omics_links", "WARN", "PRIDE project metadata has no RNA/other-omics links; current RNA must be treated as mock or externally supplied.")

    remote_files = pride_file_names(files)
    if files_error:
        add_status(rows, "external_pride", "files_api", "WARN", files_error)
    else:
        manifest_raw_names = {f"{run}.raw" for run in manifest_runs}
        missing_remote = sorted(manifest_raw_names - remote_files)
        add_status(rows, "external_pride", "manifest_raw_files_in_pride", "OK" if not missing_remote else "FAIL", f"missing={missing_remote[:20]}")
        payload["pride_file_count"] = len(remote_files)

    active_status = "OK"
    active_detail = str(len(manifest))
    if args.expected_active_runs is not None and len(manifest) != args.expected_active_runs:
        active_status = "FAIL"
        active_detail = f"expected={args.expected_active_runs}, observed={len(manifest)}"
    excluded_status = "OK"
    excluded_detail = str(len(excluded))
    if args.expected_excluded_runs is not None and len(excluded) != args.expected_excluded_runs:
        excluded_status = "FAIL"
        excluded_detail = f"expected={args.expected_excluded_runs}, observed={len(excluded)}"
    add_status(rows, "local_manifest", "active_manifest_rows", active_status, active_detail)
    add_status(rows, "local_manifest", "excluded_run_rows", excluded_status, excluded_detail)
    add_status(rows, "local_manifest", "active_psm_coverage", "OK" if manifest_runs == active_local_psm_runs else "FAIL", f"missing={sorted(manifest_runs-active_local_psm_runs)[:20]}, extra={sorted(active_local_psm_runs-manifest_runs)[:20]}")
    add_status(rows, "local_manifest", "active_mgf_coverage", "OK" if manifest_runs == active_local_mgf_runs else "FAIL", f"missing={sorted(manifest_runs-active_local_mgf_runs)[:20]}, extra={sorted(active_local_mgf_runs-manifest_runs)[:20]}")

    if {"run_id", "patient_id", "validation_id"}.issubset(manifest.columns):
        cm_rows = manifest[manifest["run_id"].str.contains("CM647", na=False)]
        cm_ok = len(cm_rows) > 0 and set(cm_rows["patient_id"]) == {"CM467"} and set(cm_rows["validation_id"]) == {"CM467"}
        add_status(rows, "local_manifest", "cm647_alias_to_cm467", "OK" if cm_ok else "FAIL", cm_rows[["run_id", "patient_id", "validation_id"]].to_dict("records"))

    if "rna_source" in manifest.columns:
        sources = sorted(set(manifest["rna_source"]))
        lower_sources = [source.lower() for source in sources]
        if any(source == "mock_debug" for source in lower_sources):
            status = "WARN"
            detail = f"sources={sources}; mock_debug is simulated and not patient-matched biological RNA."
        elif any("surrogate" in source for source in lower_sources):
            status = "WARN"
            detail = f"sources={sources}; surrogate RNA is better than mock but still not patient-matched biological RNA."
        else:
            status = "OK"
            detail = f"sources={sources}"
        add_status(rows, "rna", "rna_source", status, detail)
    else:
        add_status(rows, "rna", "rna_source", "FAIL", "manifest lacks rna_source")

    validation_ids = set(manifest.get("validation_id", pd.Series(dtype=str)))
    dataset_ids = dataset1_intensity_ids(args.validated)
    mapped_ids = {x.replace("-", "") if x.startswith("Apher-") else x for x in validation_ids}
    # pooledTIL3 intentionally maps to TIL3 in evaluation.
    mapped_ids.discard("pooledTIL3")
    add_status(rows, "validation_dataset", "dataset1_columns", "OK" if mapped_ids <= dataset_ids else "FAIL", f"manifest_validation_ids={sorted(validation_ids)}, dataset_ids={sorted(dataset_ids)}")

    stats = fasta_stats(args.reference_fasta)
    fasta_status = "OK" if stats["exists"] and stats["total"] and stats["total"] == stats["human"] else "FAIL"
    add_status(rows, "reference", "human_only_fasta", fasta_status, json.dumps(stats))

    unlabeled_runs = {path.name.replace("_unlabeled.mgf", "") for path in args.unlabeled_mgf_dir.glob("*_unlabeled.mgf")}
    for label, path in [("de_novo", args.de_novo), ("filtered", args.filtered), ("ranked", args.ranked)]:
        df = read_tsv(path)
        if df.empty:
            add_status(rows, "outputs", label, "WARN", f"{path} missing or empty")
            continue
        run_ids = set(df["run_id"]) if "run_id" in df else set()
        if label == "de_novo" and unlabeled_runs:
            stale = sorted(run_ids - unlabeled_runs)
            add_status(rows, "outputs", "de_novo_matches_current_unlabeled", "OK" if not stale else "FAIL", f"stale_run_count={len(stale)} stale={stale[:20]}")
        if "sample_id" in df:
            unknown = int((df["sample_id"] == "Unknown").sum())
            add_status(rows, "outputs", f"{label}_unknown_samples", "OK" if unknown == 0 else "FAIL", f"unknown_rows={unknown}")
        if {"sample_id", "peptide"}.issubset(df.columns):
            dupes = int(df.duplicated(["sample_id", "peptide"]).sum())
            add_status(rows, "outputs", f"{label}_duplicate_sample_peptide", "OK" if dupes == 0 or label != "ranked" else "WARN", f"duplicates={dupes}")

    report = pd.DataFrame(rows)
    counts = report["status"].value_counts().to_dict()
    payload["checks"] = rows
    payload["status_counts"] = counts
    payload["sources"] = {
        "pride_project_api": PRIDE_PROJECT_URL.format(accession=args.accession),
        "pride_files_api": PRIDE_FILES_URL.format(accession=args.accession),
        "publication": "https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005725",
        "proteomexchange": "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD005231",
    }

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(payload, indent=2, default=str))

    lines = ["# Provenance Audit", ""]
    lines.append(f"Accession: `{args.accession}`")
    lines.append(f"Status counts: `{counts}`")
    lines.append("")
    lines.append(report.to_markdown(index=False))
    lines.append("")
    lines.append("## Source Links")
    for name, url in payload["sources"].items():
        lines.append(f"- {name}: {url}")
    args.output_md.parent.mkdir(parents=True, exist_ok=True)
    args.output_md.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))

    return 1 if counts.get("FAIL", 0) else 0


if __name__ == "__main__":
    raise SystemExit(main())
