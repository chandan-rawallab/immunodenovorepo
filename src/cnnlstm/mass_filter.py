#!/usr/bin/env python3
"""
Phase 2 accuracy improvement: mass-constrained filtering of de novo predictions.

For each predicted peptide sequence, computes the theoretical monoisotopic neutral
mass and compares it against the precursor neutral mass from the MGF spectrum header.
Predictions whose mass deviates by more than `tolerance_ppm` are rejected.

Usage (standalone filter after 06_predict_denovo.py):
    python3 src/cnnlstm/mass_filter.py \
        --candidates results/de_novo_candidates.tsv \
        --mgf-dir data/mgf \
        --output results/de_novo_candidates_massfiltered.tsv \
        --tolerance-ppm 20

Or import and call filter_dataframe() directly.
"""

import argparse
import re
from pathlib import Path
from typing import Optional

import pandas as pd

# Monoisotopic residue masses (Da) — standard values from NIST
_AA_RESIDUE_MASS = {
    "A": 71.03711,
    "R": 156.10111,
    "N": 114.04293,
    "D": 115.02694,
    "C": 103.00919,
    "E": 129.04259,
    "Q": 128.05858,
    "G": 57.02146,
    "H": 137.05891,
    "I": 113.08406,
    "L": 113.08406,
    "K": 128.09496,
    "M": 131.04049,
    "F": 147.06841,
    "P": 97.05276,
    "S": 87.03203,
    "T": 101.04768,
    "W": 186.07931,
    "Y": 163.06333,
    "V": 99.06841,
}

_WATER = 18.01056     # H2O added to both termini for a full peptide
_PROTON = 1.00728     # proton mass for charge-state conversion


def peptide_neutral_mass(sequence: str) -> Optional[float]:
    """
    Return the monoisotopic neutral mass of the peptide.
    Returns None if sequence contains unrecognised amino acids.
    """
    total = _WATER
    for aa in sequence.upper():
        m = _AA_RESIDUE_MASS.get(aa)
        if m is None:
            return None
        total += m
    return total


def precursor_neutral_mass(mz: float, charge: int) -> float:
    """Convert observed precursor m/z and charge to neutral mass."""
    return (mz * charge) - (charge * _PROTON)


def ppm_error(observed: float, theoretical: float) -> float:
    """Return absolute PPM mass error."""
    if theoretical == 0:
        return float("inf")
    return abs(observed - theoretical) / theoretical * 1_000_000


# ---------------------------------------------------------------------------
# MGF precursor index builder
# ---------------------------------------------------------------------------

_TITLE_RE = re.compile(r"TITLE=(.+)", re.IGNORECASE)
_PEPMASS_RE = re.compile(r"PEPMASS=([\d.]+)\s*([\d.]+)?", re.IGNORECASE)
_CHARGE_RE = re.compile(r"CHARGE=(\d+)", re.IGNORECASE)


def build_precursor_index(mgf_dir: Path) -> dict:
    """
    Parse all MGF files in mgf_dir and return a dict:
        {spectrum_title: {"mz": float, "charge": int, "neutral_mass": float}}

    Only the TITLE, PEPMASS, and CHARGE fields are read — no peak data.
    """
    index = {}
    mgf_files = list(mgf_dir.glob("*.mgf")) + list(mgf_dir.glob("*.MGF"))
    if not mgf_files:
        print(f"[mass_filter] WARNING: No MGF files found in {mgf_dir}")
        return index

    for mgf_path in mgf_files:
        title = None
        mz = None
        charge = 2  # default if not specified
        in_scan = False
        with open(mgf_path, "r", errors="replace") as fh:
            for line in fh:
                line = line.strip()
                if line == "BEGIN IONS":
                    in_scan = True
                    title = mz = None
                    charge = 2
                    continue
                if line == "END IONS":
                    if in_scan and title is not None and mz is not None:
                        nm = precursor_neutral_mass(mz, charge)
                        index[title.strip()] = {
                            "mz": mz,
                            "charge": charge,
                            "neutral_mass": nm,
                        }
                    in_scan = False
                    continue
                if not in_scan:
                    continue
                m = _TITLE_RE.match(line)
                if m:
                    title = m.group(1)
                    continue
                m = _PEPMASS_RE.match(line)
                if m:
                    mz = float(m.group(1))
                    continue
                m = _CHARGE_RE.match(line)
                if m:
                    charge = int(m.group(1).rstrip("+"))
                    continue

    print(f"[mass_filter] Indexed {len(index):,} spectra from {len(mgf_files)} MGF files")
    return index


# ---------------------------------------------------------------------------
# Main filtering function
# ---------------------------------------------------------------------------

def filter_dataframe(
    df: pd.DataFrame,
    precursor_index: dict,
    tolerance_ppm: float = 20.0,
    spectrum_col: str = "spectrum_id",
    peptide_col: str = "peptide",
) -> pd.DataFrame:
    """
    Filter a de novo candidates DataFrame by precursor mass consistency.

    Rows where the peptide mass deviates from the spectrum precursor by more
    than `tolerance_ppm` are dropped.  Rows where the spectrum title is not
    found in the index (missing MGF data) are kept with a warning.

    Adds columns:
        mass_ppm_error  — absolute PPM deviation (NaN if spectrum not found)
        mass_valid      — bool
    """
    if spectrum_col not in df.columns:
        print(f"[mass_filter] WARNING: column '{spectrum_col}' not in DataFrame — skipping mass filter")
        df["mass_ppm_error"] = float("nan")
        df["mass_valid"] = True
        return df

    ppm_errors = []
    mass_valid = []

    for _, row in df.iterrows():
        seq = str(row.get(peptide_col, ""))
        spec_id = str(row.get(spectrum_col, ""))

        theo_mass = peptide_neutral_mass(seq)
        meta = precursor_index.get(spec_id)

        if theo_mass is None:
            # Sequence has unknown AA — reject
            ppm_errors.append(float("nan"))
            mass_valid.append(False)
        elif meta is None:
            # Spectrum not found in index — keep conservatively
            ppm_errors.append(float("nan"))
            mass_valid.append(True)
        else:
            err = ppm_error(meta["neutral_mass"], theo_mass)
            ppm_errors.append(err)
            mass_valid.append(err <= tolerance_ppm)

    df = df.copy()
    df["mass_ppm_error"] = ppm_errors
    df["mass_valid"] = mass_valid

    before = len(df)
    df = df[df["mass_valid"]].reset_index(drop=True)
    after = len(df)
    print(
        f"[mass_filter] Mass filter (±{tolerance_ppm} ppm): "
        f"{before:,} → {after:,} candidates ({before - after:,} removed)"
    )
    return df


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Mass-constraint filter for de novo candidates.")
    parser.add_argument("--candidates", type=Path, required=True, help="de_novo_candidates.tsv")
    parser.add_argument("--mgf-dir", type=Path, required=True, help="Directory of MGF files")
    parser.add_argument("--output", type=Path, required=True, help="Filtered output TSV")
    parser.add_argument("--tolerance-ppm", type=float, default=20.0,
                        help="PPM tolerance for precursor mass match (default: 20)")
    parser.add_argument("--spectrum-col", default="spectrum_id",
                        help="Column name for spectrum title in candidates TSV")
    parser.add_argument("--peptide-col", default="peptide",
                        help="Column name for predicted sequence in candidates TSV")
    args = parser.parse_args()

    print(f"Loading candidates from {args.candidates} ...")
    df = pd.read_csv(args.candidates, sep="\t", low_memory=False)
    print(f"  Loaded {len(df):,} rows")

    print(f"Building precursor index from {args.mgf_dir} ...")
    index = build_precursor_index(args.mgf_dir)

    df_filtered = filter_dataframe(
        df,
        index,
        tolerance_ppm=args.tolerance_ppm,
        spectrum_col=args.spectrum_col,
        peptide_col=args.peptide_col,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    df_filtered.to_csv(args.output, sep="\t", index=False)
    print(f"Saved {len(df_filtered):,} mass-valid candidates to {args.output}")


if __name__ == "__main__":
    main()
