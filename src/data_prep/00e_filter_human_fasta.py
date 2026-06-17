#!/usr/bin/env python3
"""Create a human-only FASTA from a mixed UniProt FASTA."""

import argparse
from pathlib import Path


def is_human_header(header: str) -> bool:
    return "OS=Homo sapiens" in header or "OX=9606" in header


def main() -> int:
    parser = argparse.ArgumentParser(description="Filter a FASTA to Homo sapiens entries only.")
    parser.add_argument("--input", type=Path, default=Path("data/reference/uniprot_human_reviewed.fasta"))
    parser.add_argument("--output", type=Path, default=Path("data/reference/uniprot_human_reviewed.only_human.fasta"))
    args = parser.parse_args()

    if not args.input.exists():
        print(f"ERROR: input FASTA not found: {args.input}")
        return 1

    args.output.parent.mkdir(parents=True, exist_ok=True)

    total = 0
    kept = 0
    write_entry = False
    with args.input.open() as src, args.output.open("w") as dst:
        for line in src:
            if line.startswith(">"):
                total += 1
                write_entry = is_human_header(line)
                kept += int(write_entry)
            if write_entry:
                dst.write(line)

    print(f"Read {total} FASTA entries.")
    print(f"Wrote {kept} human entries to {args.output}.")
    if kept == 0:
        print("ERROR: no human entries found; check FASTA header format.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
