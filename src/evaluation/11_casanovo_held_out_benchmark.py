#!/usr/bin/env python3
"""
Casanovo Held-Out Benchmark (Phase 1, Priority 2)

Runs Casanovo on the same held-out test set used by 10_evaluate_denovo_model.py
so that both models are compared on identical spectra with identical metrics.

Approach:
  For each PSM in test_set_psms.tsv:
    1. Locate the spectrum in data/mgf/<run_id>.mgf by SCANS= value.
    2. Write a single-spectrum temporary MGF file.
    3. Call: conda run -n objective3-casanovo casanovo sequence <tmp.mgf>
    4. Parse the resulting .mztab for the top-1 prediction.
    5. Evaluate against ground truth peptide.

Metrics (identical to 10_evaluate_denovo_model.py):
  - exact_peptide_accuracy         (case-sensitive)
  - exact_peptide_accuracy_il      (I and L treated as equivalent)
  - token_accuracy                 (per-position, shorter of pred/truth)
  - edit_distance_le1_rate         (fraction with edit distance <= 1)
  - per_length_accuracy            (breakdown by 8, 9, 10, 11+ mers)

Usage:
    conda activate base
    .venv/bin/python src/evaluation/11_casanovo_held_out_benchmark.py \\
        --test-psms  results/checkpoints_curated31_v2/test_set_psms.tsv \\
        --mgf-dir    data/mgf \\
        --output     results/casanovo_accuracy_curated31_v2.json \\
        --casanovo-env  objective3-casanovo \\
        --tmp-dir    /tmp/casanovo_benchmark
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Amino acid utilities
# ---------------------------------------------------------------------------

def il_collapse(seq: str) -> str:
    """Treat I and L as equivalent (same nominal mass 113.08 Da)."""
    return seq.replace("I", "L")


def edit_distance(a: str, b: str) -> int:
    """Standard Levenshtein edit distance."""
    m, n = len(a), len(b)
    dp = list(range(n + 1))
    for i in range(1, m + 1):
        prev = dp[0]
        dp[0] = i
        for j in range(1, n + 1):
            temp = dp[j]
            if a[i - 1] == b[j - 1]:
                dp[j] = prev
            else:
                dp[j] = 1 + min(prev, dp[j], dp[j - 1])
            prev = temp
    return dp[n]


def token_accuracy(pred: str, truth: str) -> float:
    """Per-position accuracy over the shorter of the two sequences."""
    length = min(len(pred), len(truth))
    if length == 0:
        return 0.0
    matches = sum(p == t for p, t in zip(pred, truth))
    return matches / length


# ---------------------------------------------------------------------------
# MGF utilities
# ---------------------------------------------------------------------------

def build_mgf_index(mgf_path: Path) -> dict[str, int]:
    """
    Build a {scan_id: byte_offset} index for fast spectrum retrieval.
    Handles both SCANS= and TITLE=...scan=N formats.
    """
    index: dict[str, int] = {}
    with open(mgf_path, "rb") as fh:
        while True:
            offset = fh.tell()
            line = fh.readline()
            if not line:
                break
            decoded = line.decode("utf-8", errors="replace").strip()
            if decoded == "BEGIN IONS":
                # Scan ahead to find scan id
                scan_offset = fh.tell()
                scan_id = None
                for _ in range(20):
                    inner = fh.readline().decode("utf-8", errors="replace").strip()
                    if inner.startswith("SCANS="):
                        scan_id = inner.split("=", 1)[1].strip()
                        break
                    if inner.startswith("TITLE="):
                        title = inner.split("=", 1)[1]
                        # Try scan= substring
                        if "scan=" in title:
                            scan_id = title.split("scan=")[-1].split()[0].strip()
                        else:
                            scan_id = title.split()[-1].strip()
                        # Also try SCANS on next lines before giving up
                    if inner == "END IONS":
                        break
                if scan_id:
                    index[scan_id] = offset
                fh.seek(scan_offset)
    return index


def extract_spectrum_block(mgf_path: Path, byte_offset: int) -> str:
    """Return the raw MGF block (BEGIN IONS ... END IONS) at byte_offset."""
    lines = []
    with open(mgf_path, "rb") as fh:
        fh.seek(byte_offset)
        in_block = False
        for _ in range(5000):
            raw = fh.readline()
            if not raw:
                break
            line = raw.decode("utf-8", errors="replace")
            stripped = line.strip()
            if stripped == "BEGIN IONS":
                in_block = True
            if in_block:
                lines.append(line)
            if stripped == "END IONS":
                break
    return "".join(lines)


# ---------------------------------------------------------------------------
# Casanovo invocation
# ---------------------------------------------------------------------------

def run_casanovo_on_mgf(tmp_mgf: Path, output_dir: Path, casanovo_env: str) -> Path | None:
    """
    Run: conda run -n <env> casanovo sequence <tmp_mgf> --output_root <output_dir>
    Returns path to the resulting .mztab file, or None on failure.
    """
    cmd = [
        "conda", "run", "-n", casanovo_env,
        "casanovo", "sequence",
        str(tmp_mgf),
        "--output_root", str(output_dir),
    ]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=120,
    )
    # Find .mztab output
    mztab_files = list(output_dir.glob("*.mztab"))
    if not mztab_files:
        return None
    return mztab_files[0]


def parse_mztab_top1(mztab_path: Path) -> str | None:
    """Extract the top-1 peptide sequence from a Casanovo mzTab PSM section."""
    if not mztab_path or not mztab_path.exists():
        return None
    with open(mztab_path, "r", errors="replace") as fh:
        for line in fh:
            if line.startswith("PSM\t"):
                parts = line.strip().split("\t")
                # mzTab PSM columns: PSM, sequence, PSM_ID, accession, ...
                if len(parts) >= 2:
                    seq = parts[1].strip().upper()
                    # Strip any mod notation like [+16]
                    seq = "".join(c for c in seq if c.isalpha())
                    return seq if seq else None
    return None


# ---------------------------------------------------------------------------
# Main evaluation loop
# ---------------------------------------------------------------------------

def evaluate(
    test_psms_path: Path,
    mgf_dir: Path,
    output_path: Path,
    casanovo_env: str,
    tmp_dir: Path,
    max_spectra: int | None,
) -> dict:
    test_df = pd.read_csv(test_psms_path, sep="\t", dtype=str)
    print(f"Loaded {len(test_df)} test PSMs from {test_psms_path}")

    if max_spectra:
        test_df = test_df.head(max_spectra)
        print(f"  → Running on first {max_spectra} spectra (--max-spectra limit)")

    # Pre-build MGF indexes for all unique runs in the test set
    run_ids = test_df["run_id"].unique()
    mgf_indexes: dict[str, dict[str, int]] = {}
    mgf_paths: dict[str, Path] = {}
    for run_id in run_ids:
        mgf_path = mgf_dir / f"{run_id}.mgf"
        if not mgf_path.exists():
            print(f"  WARNING: MGF not found for run {run_id} at {mgf_path}")
            continue
        print(f"  Indexing {mgf_path.name} ...")
        mgf_indexes[run_id] = build_mgf_index(mgf_path)
        mgf_paths[run_id] = mgf_path
    print(f"Indexed {len(mgf_indexes)} MGF files.")

    tmp_dir.mkdir(parents=True, exist_ok=True)

    results = []
    skipped = 0

    for idx, row in test_df.iterrows():
        run_id = str(row["run_id"]).strip()
        spectrum_id = str(row["spectrum_id"]).strip()
        truth = str(row["peptide"]).strip().upper()

        if run_id not in mgf_indexes:
            skipped += 1
            continue

        index = mgf_indexes[run_id]
        if spectrum_id not in index:
            skipped += 1
            continue

        # Extract spectrum block
        block = extract_spectrum_block(mgf_paths[run_id], index[spectrum_id])
        if not block:
            skipped += 1
            continue

        # Write single-spectrum MGF
        spectrum_tmp = tmp_dir / f"spec_{idx}.mgf"
        spectrum_tmp.write_text(block)

        # Create a clean per-spectrum output dir
        out_subdir = tmp_dir / f"out_{idx}"
        out_subdir.mkdir(exist_ok=True)

        try:
            mztab_path = run_casanovo_on_mgf(spectrum_tmp, out_subdir, casanovo_env)
            prediction = parse_mztab_top1(mztab_path)
        except subprocess.TimeoutExpired:
            prediction = None
        except Exception as e:
            prediction = None

        # Clean up temp files
        spectrum_tmp.unlink(missing_ok=True)
        shutil.rmtree(out_subdir, ignore_errors=True)

        results.append({
            "run_id": run_id,
            "spectrum_id": spectrum_id,
            "truth": truth,
            "prediction": prediction or "",
            "truth_len": len(truth),
        })

        if (idx + 1) % 50 == 0:
            done = len(results)
            exact = sum(r["truth"] == r["prediction"] for r in results if r["prediction"])
            print(f"  [{done}/{len(test_df)}] Running exact accuracy: {exact/done*100:.1f}%")

    print(f"\nSkipped {skipped} spectra (MGF not found or scan ID mismatch).")
    print(f"Evaluated {len(results)} spectra.")

    # ---- Compute metrics ----
    if not results:
        print("ERROR: No results to evaluate.")
        return {}

    res_df = pd.DataFrame(results)
    valid = res_df[res_df["prediction"] != ""]

    n_total = len(res_df)
    n_valid = len(valid)

    exact = (valid["truth"] == valid["prediction"]).sum()
    exact_il = (valid["truth"].apply(il_collapse) == valid["prediction"].apply(il_collapse)).sum()

    token_accs = valid.apply(
        lambda r: token_accuracy(r["prediction"], r["truth"]), axis=1
    )

    edit_dists = valid.apply(
        lambda r: edit_distance(r["prediction"], r["truth"]), axis=1
    )
    edit_le1 = (edit_dists <= 1).sum()

    # Per-length breakdown
    per_length: dict[str, dict] = {}
    for length in [8, 9, 10, 11]:
        subset = valid[valid["truth_len"] == length]
        if len(subset) == 0:
            continue
        per_length[str(length)] = {
            "n": len(subset),
            "exact_accuracy": (subset["truth"] == subset["prediction"]).sum() / len(subset),
        }
    subset_12plus = valid[valid["truth_len"] >= 12]
    if len(subset_12plus) > 0:
        per_length["12+"] = {
            "n": len(subset_12plus),
            "exact_accuracy": (subset_12plus["truth"] == subset_12plus["prediction"]).sum() / len(subset_12plus),
        }

    metrics = {
        "model": "casanovo",
        "test_set": str(test_psms_path),
        "n_total_psms": n_total,
        "n_casanovo_predicted": n_valid,
        "n_skipped": skipped,
        "exact_peptide_accuracy": round(exact / n_valid, 4) if n_valid else 0.0,
        "exact_peptide_accuracy_il": round(exact_il / n_valid, 4) if n_valid else 0.0,
        "mean_token_accuracy": round(float(token_accs.mean()), 4) if n_valid else 0.0,
        "edit_distance_le1_rate": round(edit_le1 / n_valid, 4) if n_valid else 0.0,
        "per_length_accuracy": per_length,
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(metrics, f, indent=2)

    print(f"\n=== Casanovo Held-Out Benchmark Results ===")
    print(f"  Spectra evaluated:          {n_valid} / {n_total}")
    print(f"  Exact peptide accuracy:     {metrics['exact_peptide_accuracy']*100:.1f}%")
    print(f"  Exact accuracy (I=L):       {metrics['exact_peptide_accuracy_il']*100:.1f}%")
    print(f"  Mean token accuracy:        {metrics['mean_token_accuracy']*100:.1f}%")
    print(f"  Edit-dist ≤ 1 rate:         {metrics['edit_distance_le1_rate']*100:.1f}%")
    print(f"\nSaved to {output_path}")

    return metrics


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Benchmark Casanovo on the held-out test set from 05_train_denovo_model.py."
    )
    parser.add_argument(
        "--test-psms",
        type=Path,
        default=Path("results/checkpoints_curated31_v2/test_set_psms.tsv"),
        help="Path to test_set_psms.tsv generated during training.",
    )
    parser.add_argument(
        "--mgf-dir",
        type=Path,
        default=Path("data/mgf"),
        help="Directory containing labeled MGF files (one per run).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/casanovo_accuracy_curated31_v2.json"),
        help="Path for the output JSON metrics file.",
    )
    parser.add_argument(
        "--casanovo-env",
        default="objective3-casanovo",
        help="Name of the conda environment containing Casanovo.",
    )
    parser.add_argument(
        "--tmp-dir",
        type=Path,
        default=Path("/tmp/casanovo_benchmark"),
        help="Scratch directory for temporary single-spectrum MGF files.",
    )
    parser.add_argument(
        "--max-spectra",
        type=int,
        default=None,
        help="Limit evaluation to N spectra (for fast smoke testing).",
    )
    args = parser.parse_args()

    evaluate(
        test_psms_path=args.test_psms,
        mgf_dir=args.mgf_dir,
        output_path=args.output,
        casanovo_env=args.casanovo_env,
        tmp_dir=args.tmp_dir,
        max_spectra=args.max_spectra,
    )


if __name__ == "__main__":
    main()
