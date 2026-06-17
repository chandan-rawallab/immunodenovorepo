#!/usr/bin/env python3
"""Quarantine generated/demo clutter without deleting files."""

import argparse
import shutil
from pathlib import Path

import pandas as pd


QUARANTINE_RULES = {
    "generated_presentations": [
        "better_objective3_presentation.pptx",
        "objective3_presentation.pptx",
        "create_presentation.py",
        "make_beautiful_presentation.js",
        "src/create_premium_presentation.py",
        "temp_slides",
    ],
    "node_artifacts": [
        "node_modules",
        "package.json",
        "package-lock.json",
    ],
    "demo_cm467_circular": [
        "results/CM467_annotated_alleles.tsv",
        "results/CM467_extracted_data.tsv",
        "results/CM467_filtered_neoepitopes.tsv",
        "results/cm467_evaluation_report.md",
        "src/postprocess/cm467_annotate.py",
        "src/postprocess/cm467_extract.py",
        "src/postprocess/cm467_filter.py",
    ],
    "scratch": [
        "scratch",
        "dotnet-install.sh",
    ],
    "stale_results": [
        "results/20260516_180013_denovo_run",
        "results/20260519_163037_denovo_run",
        "results/20260519_181846_denovo_run",
        "results/20260520_093618_denovo_run",
        "results/20260520_105553_denovo_run",
        "results/20260520_105906_denovo_run",
    ],
}


def unique_destination(root: Path, source: Path) -> Path:
    dest = root / source.name
    if not dest.exists():
        return dest
    suffix = 1
    while True:
        candidate = root / f"{source.name}.{suffix}"
        if not candidate.exists():
            return candidate
        suffix += 1


def main() -> int:
    parser = argparse.ArgumentParser(description="Move generated/demo clutter into archive/quarantine.")
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--quarantine-root", type=Path, default=Path("archive/quarantine"))
    parser.add_argument("--apply", action="store_true", help="Actually move files. Default is dry-run.")
    args = parser.parse_args()

    rows = []
    for bucket, paths in QUARANTINE_RULES.items():
        bucket_dir = args.root / args.quarantine_root / bucket
        for raw_path in paths:
            source = args.root / raw_path
            if not source.exists():
                continue
            dest = unique_destination(bucket_dir, source)
            rows.append(
                {
                    "bucket": bucket,
                    "original_path": str(source),
                    "quarantine_path": str(dest),
                    "size_bytes": source.stat().st_size if source.is_file() else "",
                    "reason": bucket,
                    "applied": bool(args.apply),
                }
            )
            print(f"{'MOVE' if args.apply else 'DRY'} {source} -> {dest}")
            if args.apply:
                bucket_dir.mkdir(parents=True, exist_ok=True)
                shutil.move(str(source), str(dest))

    manifest_path = args.root / args.quarantine_root / "MANIFEST.tsv"
    if rows:
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        new_df = pd.DataFrame(rows)
        if args.apply and manifest_path.exists():
            old_df = pd.read_csv(manifest_path, sep="\t")
            new_df = pd.concat([old_df, new_df], ignore_index=True)
        if args.apply:
            new_df.to_csv(manifest_path, sep="\t", index=False)
            print(f"Wrote quarantine manifest: {manifest_path}")
        else:
            print(f"Dry-run only. Re-run with --apply to move files and write {manifest_path}.")
    else:
        print("No matching files found for quarantine rules.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
