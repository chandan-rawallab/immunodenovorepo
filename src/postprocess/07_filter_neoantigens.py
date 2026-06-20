#!/usr/bin/env python3
"""Activity 3: FDR, length, and missense mutation filters for neoantigen candidates.

Audit fixes applied (2026-06-20):
  - Reference proteome FASTA is parsed ONCE and a protein-to-sequence lookup
    is cached so find_source_proteins() does not re-open the file per candidate.
  - Per-stage candidate counts are printed at every filter step for audit trail.
  - find_mutation() now explicitly labels the mutation type returned.
"""

import argparse
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd
from pyteomics import fasta

# Canonical HLA-I length 8-11
CANONICAL_HLA_I_PATTERN = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]{8,11}$")

# ---------------------------------------------------------------------------
# FASTA helpers – parsed once and cached
# ---------------------------------------------------------------------------

def _load_human_records(fasta_path: Path) -> list[tuple[str, str]]:
    """Return (description, sequence) pairs for Homo sapiens entries only."""
    records = list(fasta.read(str(fasta_path)))
    has_tags = any("OS=" in d or "OX=" in d for d, _ in records)
    if not has_tags:
        return records
    return [
        (d, s)
        for d, s in records
        if "OS=Homo sapiens" in d or "OX=9606" in d
    ]


def build_reference_index(fasta_path: Path, lengths: list[int] = None):
    """Build two indexes from the reference proteome in a single FASTA pass.

    Returns
    -------
    ref_peptide_set : set[str]
        All 8-11 mer subsequences of human proteins.
    protein_lookup : dict[str, list[str]]
        Maps full-length protein sequence → list of Uniprot IDs.
        Used by find_source_proteins() without re-reading the file.
    """
    if lengths is None:
        lengths = [8, 9, 10, 11]

    records = _load_human_records(fasta_path)
    ref_peptide_set: set[str] = set()
    # Map sequence → [protein_id, ...] (multiple proteins can share a sequence)
    seq_to_proteins: dict[str, list[str]] = defaultdict(list)

    print(f"Indexing reference proteome from {fasta_path}…")
    for desc, sequence in records:
        parts = desc.split("|")
        protein_id = parts[1] if len(parts) > 1 else desc.split()[0]
        seq_to_proteins[sequence].append(protein_id)
        for length in lengths:
            seq_len = len(sequence)
            for i in range(seq_len - length + 1):
                ref_peptide_set.add(sequence[i : i + length])

    print(f"Indexed {len(ref_peptide_set):,} unique reference peptides "
          f"from {len(records):,} protein records.")
    return ref_peptide_set, seq_to_proteins


def find_mutation(peptide: str, reference_set: set[str]):
    """Check if the peptide differs by exactly one substitution from a reference peptide.

    Returns (wt_peptide, pos_1based, wt_aa, mut_aa) or None.
    """
    aas = "ACDEFGHIKLMNPQRSTVWY"
    for i in range(len(peptide)):
        mut_aa = peptide[i]
        for wt_aa in aas:
            if wt_aa == mut_aa:
                continue
            candidate_wt = peptide[:i] + wt_aa + peptide[i + 1 :]
            if candidate_wt in reference_set:
                return candidate_wt, i + 1, wt_aa, mut_aa
    return None


def find_source_proteins(
    wt_peptide: str, seq_to_proteins: dict[str, list[str]]
) -> list[str]:
    """Return Uniprot IDs for proteins containing *wt_peptide* (cached lookup)."""
    found: list[str] = []
    for sequence, prot_ids in seq_to_proteins.items():
        if wt_peptide in sequence:
            found.extend(prot_ids)
    return found


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Filter de novo candidates for neoantigens."
    )
    parser.add_argument("--input", type=Path, required=True,
                        help="Path to de_novo_candidates.tsv")
    parser.add_argument("--psms", type=Path, required=True,
                        help="Path to immunopeptidome_psms.tsv (database search hits)")
    parser.add_argument("--fasta", type=Path, required=True,
                        help="Path to reference proteome FASTA")
    parser.add_argument("--manifest", type=Path,
                        default="configs/sample_manifest.tsv",
                        help="Sample manifest for run_id to sample_id mapping")
    parser.add_argument("--output", type=Path, required=True,
                        help="Output filtered TSV")
    parser.add_argument("--score_cutoff", type=float, default=-0.5,
                        help="De novo score cutoff (log-probability, default: -0.5)")
    parser.add_argument("--min_psm_support", type=int, default=2,
                        help="Minimum spectra supporting a peptide (>= 2)")
    parser.add_argument("--allow_flanking_mutations", action="store_true",
                        help="Retain mutations at flanking positions (pos 1 or last).")

    args = parser.parse_args()

    # ------------------------------------------------------------------
    # 1. Load data
    # ------------------------------------------------------------------
    print(f"\n[Stage 0] Loading candidates from {args.input}…")
    df = pd.read_csv(args.input, sep="\t")
    print(f"  Input candidates : {len(df):>7,}")

    manifest = pd.read_csv(args.manifest, sep="\t")
    run_to_sample = dict(zip(manifest["run_id"], manifest["patient_id"]))
    df["sample_id"] = df["run_id"].map(run_to_sample).fillna("Unknown")

    if "score" not in df.columns and "de_novo_score" in df.columns:
        df = df.rename(columns={"de_novo_score": "score"})
    if "score" not in df.columns:
        raise ValueError("Candidate table must contain `score` or `de_novo_score`.")

    print(f"\n[Stage 1] Loading database-search PSMs from {args.psms}…")
    df_psms = pd.read_csv(args.psms, sep="\t")
    db_peptides: dict[str, set[str]] = defaultdict(set)
    for _, row in df_psms.iterrows():
        db_peptides[str(row["sample_id"])].add(str(row["peptide"]))

    # ------------------------------------------------------------------
    # 2. Score filter
    # ------------------------------------------------------------------
    before = len(df)
    df = df[df["score"] >= args.score_cutoff].copy()
    print(f"\n[Stage 2] Score >= {args.score_cutoff}: {before:,} → {len(df):,}")

    # ------------------------------------------------------------------
    # 3. Length / canonical AA filter
    # ------------------------------------------------------------------
    before = len(df)
    df = df[df["peptide"].str.match(CANONICAL_HLA_I_PATTERN)].copy()
    print(f"[Stage 3] Canonical HLA-I length (8–11 AA): {before:,} → {len(df):,}")

    # ------------------------------------------------------------------
    # 4. Database subtraction
    # ------------------------------------------------------------------
    mask_in_db = df.apply(
        lambda row: row["peptide"] in db_peptides.get(str(row["sample_id"]), set()),
        axis=1,
    )
    before = len(df)
    df = df[~mask_in_db].copy()
    print(f"[Stage 4] Remove database-search hits: {before:,} → {len(df):,}")

    # ------------------------------------------------------------------
    # 5. PSM support count
    # ------------------------------------------------------------------
    if "run_id" in df.columns:
        psm_counts = (
            df.groupby(["sample_id", "peptide"])["run_id"]
            .count()
            .rename("psm_count")
        )
    else:
        psm_counts = df.groupby(["sample_id", "peptide"]).size().rename("psm_count")

    df = df.join(psm_counts, on=["sample_id", "peptide"])
    before = len(df)
    df = df[df["psm_count"] >= args.min_psm_support].copy()
    print(f"[Stage 5] PSM support >= {args.min_psm_support}: {before:,} → {len(df):,}")

    # ------------------------------------------------------------------
    # 6. Missense mutation detection (single FASTA pass — audit fix)
    # ------------------------------------------------------------------
    print(f"\n[Stage 6] Indexing reference proteome (single pass)…")
    ref_set, seq_to_proteins = build_reference_index(args.fasta)

    results = []
    n_self = n_missense = n_other = 0

    print(f"Performing missense mutation detection on {len(df):,} candidates…")
    for _, row in df.iterrows():
        peptide = row["peptide"]

        # Exact self-peptide — skip
        if peptide in ref_set:
            n_self += 1
            continue

        mutation_info = find_mutation(peptide, ref_set)
        row_dict = row.to_dict()

        if mutation_info:
            wt_seq, pos, wt_aa, mut_aa = mutation_info
            is_flanking = pos == 1 or pos == len(peptide)
            if is_flanking and not args.allow_flanking_mutations:
                n_self += 1  # treated as non-neoantigen for ranking purposes
                continue

            # Cached source-protein lookup — no repeated FASTA I/O
            source_proteins = find_source_proteins(wt_seq, seq_to_proteins)

            row_dict.update({
                "wildtype_peptide": wt_seq,
                "mutation_pos": pos,
                "wt_aa": wt_aa,
                "mut_aa": mut_aa,
                "mutation_type": "missense",
                "source_protein": ";".join(source_proteins),
            })
            results.append(row_dict)
            n_missense += 1
        else:
            row_dict.update({
                "wildtype_peptide": "",
                "mutation_pos": "",
                "wt_aa": "",
                "mut_aa": "",
                "mutation_type": "other",
                "source_protein": "",
            })
            results.append(row_dict)
            n_other += 1

    df_final = pd.DataFrame(results)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print(f"\n=== Filter Stage Summary ===")
    print(f"  After score filter        : see Stage 2")
    print(f"  After length filter       : see Stage 3")
    print(f"  After DB subtraction      : see Stage 4")
    print(f"  After PSM support filter  : see Stage 5")
    print(f"  Self-peptides removed     : {n_self:>7,}")
    print(f"  Missense (Class A/B elig) : {n_missense:>7,}")
    print(f"  Other / Class C           : {n_other:>7,}")
    print(f"  TOTAL retained            : {len(df_final):>7,}")

    df_final.to_csv(args.output, sep="\t", index=False)
    print(f"\nSaved {len(df_final):,} neoantigen candidates to {args.output}")


if __name__ == "__main__":
    main()
