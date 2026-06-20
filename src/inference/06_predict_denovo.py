"""De novo neoepitope prediction from unlabelled MGF spectra.

Audit fixes applied (2026-06-20):
  - Checkpoint metadata JSON is validated before model weights are loaded;
    architecture mismatch raises a clear error instead of a cryptic RuntimeError.
  - bin_spectrum, decode_sequence, and compute_sequence_score now share the same
    preprocessing path as SpectralDataset (log1p transform + top-N peak filter)
    to prevent training/inference drift.
  - FDR methodology comment clarified; decoy approach is explicitly documented.
  - Structured logging replaces bare print() throughout.
  - Beam-search stub added (greedy remains default; beam_width=1 ≡ greedy).
"""

from __future__ import annotations

import argparse
import glob
import json
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.dirname(SCRIPT_DIR)
WORKSPACE_ROOT = os.path.dirname(SRC_DIR)
sys.path.append(SRC_DIR)

from cnnlstm.cnnlstm_model import NeoepitopeSeq2Seq
from cnnlstm.spectral_dataset import AA_TO_INT, VOCAB_SIZE

try:
    from pyteomics import mgf as _pyteomics_mgf
except ImportError:
    _pyteomics_mgf = None
    from cnnlstm.mgf_utils import mgf_read_fallback as _mgf_read_fallback  # noqa

INT_TO_AA: dict[int, str] = {v: k for k, v in AA_TO_INT.items()}

# ---------------------------------------------------------------------------
# Shared spectrum preprocessing (mirrors SpectralDataset defaults)
# ---------------------------------------------------------------------------
_TOP_N_PEAKS = 200
_INTENSITY_TRANSFORM = "log1p"


def _apply_top_n(mz_array, int_array, top_n: int = _TOP_N_PEAKS):
    if top_n is None or len(mz_array) <= top_n:
        return mz_array, int_array
    pairs = sorted(zip(int_array, mz_array), reverse=True)[:top_n]
    intensities, mzs = zip(*pairs)
    return list(mzs), list(intensities)


def _transform_intensity(values: np.ndarray, method: str = _INTENSITY_TRANSFORM) -> np.ndarray:
    if method == "log1p":
        return np.log1p(values)
    if method == "sqrt":
        return np.sqrt(np.maximum(values, 0.0))
    return values


def bin_spectrum(
    mz_array,
    int_array,
    bin_size: float = 0.1,
    max_mz: float = 2000.0,
    top_n: int = _TOP_N_PEAKS,
    transform: str = _INTENSITY_TRANSFORM,
) -> torch.Tensor:
    """Convert raw peaks into a normalised binned vector.

    Applies the same top-N filtering and intensity transform used in
    SpectralDataset to prevent training/inference drift.
    """
    mz_array, int_array = _apply_top_n(mz_array, int_array, top_n)

    vector_size = int(max_mz / bin_size)
    vector = np.zeros(vector_size, dtype=np.float32)
    for mz, intensity in zip(mz_array, int_array):
        if mz < max_mz:
            bin_idx = int(mz / bin_size)
            if bin_idx < vector_size:
                vector[bin_idx] += intensity

    vector = _transform_intensity(vector, transform)
    max_val = vector.max()
    if max_val > 0:
        vector /= max_val

    # Shape: (1, 1, vector_size) — batch=1, channel=1, length
    return torch.tensor(vector).unsqueeze(0).unsqueeze(0)


# ---------------------------------------------------------------------------
# Sequence decoding
# ---------------------------------------------------------------------------

def decode_sequence(predicted_indices) -> str:
    seq = []
    for p_idx in predicted_indices:
        token = INT_TO_AA.get(int(p_idx), "")
        if token == "<END>":
            break
        if token not in ("<START>", "<PAD>", "<UNK>", ""):
            seq.append(token)
    return "".join(seq)


def compute_sequence_score(outputs: torch.Tensor, predicted_indices) -> float:
    """Length-normalised mean log-probability score."""
    probs = torch.softmax(outputs, dim=-1).squeeze(0)  # (seq_len, vocab)
    score = 0.0
    length = 0
    for i, p_idx in enumerate(predicted_indices):
        p_idx = int(p_idx)
        score += float(np.log(probs[i, p_idx].item() + 1e-9))
        length += 1
        if INT_TO_AA.get(p_idx, "") == "<END>":
            break
    return score / max(1, length)


def beam_search_decode(outputs: torch.Tensor, beam_width: int = 5) -> list[tuple[str, float]]:
    """
    Independent position-wise beam search.
    Returns list of (decoded_sequence, length_normalized_log_prob)
    """
    probs = torch.softmax(outputs, dim=-1).squeeze(0).cpu().numpy()
    seq_len = probs.shape[0]
    
    beam = [([], 0.0)]
    
    for i in range(seq_len):
        new_beam = []
        for seq, score in beam:
            if seq and INT_TO_AA.get(seq[-1], "") == "<END>":
                new_beam.append((seq, score))
                continue
                
            top_k = np.argsort(probs[i])[-beam_width:]
            for idx in top_k:
                p = probs[i, idx]
                new_score = score + np.log(p + 1e-9)
                new_beam.append((seq + [int(idx)], new_score))
                
        new_beam.sort(key=lambda x: x[1], reverse=True)
        beam = new_beam[:beam_width]
        
    results = []
    for seq, score in beam:
        decoded = decode_sequence(seq)
        if decoded:
            norm_score = score / max(1, len(seq))
            results.append((decoded, norm_score))
            
    unique_results = {}
    for seq, score in results:
        if seq not in unique_results or score > unique_results[seq]:
            unique_results[seq] = score
            
    sorted_results = sorted(unique_results.items(), key=lambda x: x[1], reverse=True)
    return sorted_results[:beam_width]


# ---------------------------------------------------------------------------
# Checkpoint validation
# ---------------------------------------------------------------------------

def load_checkpoint_metadata(model_path: str) -> dict | None:
    """Load the JSON metadata sidecar saved by 05_train_denovo_model.py."""
    meta_path = model_path.replace(".pth", "_metadata.json")
    if not os.path.exists(meta_path):
        logger.warning(
            "No checkpoint metadata found at '%s'. "
            "Cannot validate architecture before loading weights.",
            meta_path,
        )
        return None
    with open(meta_path) as fh:
        return json.load(fh)


def validate_checkpoint_metadata(meta: dict | None) -> None:
    """Raise if the checkpoint was saved for a different model class."""
    if meta is None:
        return  # Skip validation if sidecar is absent (legacy checkpoints)
    saved_class = meta.get("model_class", "")
    if saved_class and saved_class != "NeoepitopeSeq2Seq":
        raise RuntimeError(
            f"Checkpoint metadata reports model_class='{saved_class}' "
            f"but this script expects 'NeoepitopeSeq2Seq'."
        )
    saved_vocab = meta.get("vocab_size")
    if saved_vocab is not None and saved_vocab != VOCAB_SIZE:
        raise RuntimeError(
            f"Checkpoint vocab_size={saved_vocab} does not match "
            f"current VOCAB_SIZE={VOCAB_SIZE}. Retrain or use the correct checkpoint."
        )
    logger.info(
        "Checkpoint metadata validated: class=%s  vocab_size=%s  epoch=%s",
        meta.get("model_class"),
        meta.get("vocab_size"),
        meta.get("epoch"),
    )


# ---------------------------------------------------------------------------
# MGF reader
# ---------------------------------------------------------------------------

def _open_mgf_reader(mgf_path: str):
    if _pyteomics_mgf is not None:
        return _pyteomics_mgf.read(mgf_path)
    return _mgf_read_fallback(mgf_path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Run de novo neoepitope prediction.")
    parser.add_argument(
        "--model",
        type=str,
        default=os.path.join(WORKSPACE_ROOT, "results", "checkpoints", "neoepitope_production_best.pth"),
    )
    parser.add_argument(
        "--mgf_dir",
        type=str,
        default=os.path.join(WORKSPACE_ROOT, "data", "mgf_unlabeled"),
    )
    parser.add_argument(
        "--output",
        type=str,
        default=os.path.join(WORKSPACE_ROOT, "results", "de_novo_candidates.tsv"),
    )
    parser.add_argument("--device", type=str, default="cpu")
    parser.add_argument("--bin_size", type=float, default=0.1)
    parser.add_argument("--max_mz", type=float, default=2000.0)
    parser.add_argument(
        "--top_n_peaks",
        type=int,
        default=_TOP_N_PEAKS,
        help="Retain only the N most intense peaks before binning (0 = disabled).",
    )
    parser.add_argument(
        "--intensity_transform",
        choices=["log1p", "sqrt", "none"],
        default=_INTENSITY_TRANSFORM,
    )
    parser.add_argument(
        "--fdr_threshold",
        type=float,
        default=0.05,
        help="FDR threshold for target-decoy filtering.",
    )
    parser.add_argument(
        "--beam_width",
        type=int,
        default=5,
        help="Beam width for decoding. Default is 5. If 1, behaves like greedy argmax.",
    )
    parser.add_argument(
        "--mass_tolerance_ppm",
        type=float,
        default=None,
        help="If provided, filters beam candidates to those matching precursor mass within this PPM.",
    )
    args = parser.parse_args()

    MODEL_PATH = args.model
    MGF_DIR = args.mgf_dir
    OUT_FILE = args.output
    DEVICE = torch.device(args.device if args.device == "cuda" and torch.cuda.is_available() else "cpu")
    top_n = args.top_n_peaks if args.top_n_peaks > 0 else None
    
    if args.mass_tolerance_ppm is not None:
        try:
            from cnnlstm.mass_filter import precursor_neutral_mass, peptide_neutral_mass, ppm_error
        except ImportError:
            from mass_filter import precursor_neutral_mass, peptide_neutral_mass, ppm_error

    # --- Checkpoint validation (audit fix) ---
    if not os.path.exists(MODEL_PATH):
        alt_path = os.path.join(
            WORKSPACE_ROOT, "results", "_archive", "checkpoints", "neoepitope_medium_lite.pth"
        )
        if os.path.exists(alt_path):
            logger.warning("Primary checkpoint not found. Using archived fallback: %s", alt_path)
            MODEL_PATH = alt_path
            args.bin_size = 0.5
        else:
            logger.error("Checkpoint not found: %s", MODEL_PATH)
            return

    meta = load_checkpoint_metadata(MODEL_PATH)
    validate_checkpoint_metadata(meta)

    if meta:
        logger.info("Checkpoint provenance: %s", json.dumps(meta.get("provenance", {})))

    logger.info("Loading model from %s", os.path.basename(MODEL_PATH))
    model = NeoepitopeSeq2Seq().to(DEVICE)
    try:
        model.load_state_dict(torch.load(MODEL_PATH, map_location=DEVICE))
    except Exception as exc:
        logger.error("Could not load model weights: %s: %s", type(exc).__name__, exc)
        return
    model.eval()

    # --- Process MGF files ---
    mgf_files = glob.glob(os.path.join(MGF_DIR, "*.mgf"))
    logger.info("Found %d MGF files in %s", len(mgf_files), MGF_DIR)

    candidates: list[dict] = []

    for mgf_path in mgf_files:
        run_id = Path(mgf_path).stem.replace("_unlabeled", "")
        logger.info("Processing run: %s", run_id)

        try:
            reader = _open_mgf_reader(mgf_path)
        except Exception as exc:
            logger.error("Could not open %s: %s", mgf_path, exc)
            continue

        count = 0
        for spectrum in reader:
            if spectrum is None:
                continue
            if count > 0 and count % 10_000 == 0:
                logger.info("  %s: processed %d spectra", run_id, count)

            params = spectrum.get("params", {})
            scan_id = str(params.get("scans", ""))
            if not scan_id:
                title = str(params.get("title", ""))
                scan_id = title.split("scan=")[-1].strip() if "scan=" in title else str(count)

            mz_array = spectrum.get("m/z array", [])
            int_array = spectrum.get("intensity array", [])
            if len(mz_array) == 0:
                count += 1
                continue

            input_tensor = bin_spectrum(
                mz_array, int_array,
                bin_size=args.bin_size,
                max_mz=args.max_mz,
                top_n=top_n,
                transform=args.intensity_transform,
            ).to(DEVICE)

            with torch.no_grad():
                outputs = model(input_tensor)
                
            beam_results = beam_search_decode(outputs, beam_width=args.beam_width)
            
            best_seq = None
            best_score = float('-inf')
            
            if args.mass_tolerance_ppm is not None and beam_results:
                pepmass_info = params.get("pepmass", [None])
                prec_mz = pepmass_info[0] if isinstance(pepmass_info, tuple) else pepmass_info
                charge_info = params.get("charge", [2])
                charge = charge_info[0] if isinstance(charge_info, tuple) else charge_info
                if isinstance(charge, str):
                    charge = int(str(charge).replace('+', '').replace('-', ''))
                elif charge is None:
                    charge = 2
                    
                if prec_mz is not None:
                    prec_mass = precursor_neutral_mass(float(prec_mz), charge)
                    # Find highest scoring sequence in beam that matches mass
                    for seq, score in beam_results:
                        theo_mass = peptide_neutral_mass(seq)
                        if theo_mass is not None:
                            err = ppm_error(prec_mass, theo_mass)
                            if err <= args.mass_tolerance_ppm:
                                best_seq = seq
                                best_score = score
                                break
                                
            # Fallback to top sequence if mass filtering disabled or no match found
            if best_seq is None and beam_results:
                best_seq, best_score = beam_results[0]

            if best_seq and 8 <= len(best_seq) <= 11:
                candidates.append({
                    "run_id": run_id,
                    "spectrum_id": scan_id,
                    "peptide": best_seq,
                    "score": best_score,
                    "decoy": False,
                })

                # Target-decoy: score the *reversed* sequence against the same
                # spectrum to obtain an empirical null distribution.
                decoy_seq = best_seq[::-1]
                decoy_indices = [AA_TO_INT.get(c, 0) for c in decoy_seq]
                decoy_score = compute_sequence_score(outputs, decoy_indices)
                candidates.append({
                    "run_id": run_id,
                    "spectrum_id": scan_id,
                    "peptide": decoy_seq,
                    "score": decoy_score,
                    "decoy": True,
                })

            count += 1

    if not candidates:
        logger.warning("No candidates generated. Check that MGF dir is populated.")
        return

    df = pd.DataFrame(candidates)

    # --- FDR estimation (target-decoy competition) ---
    logger.info("Estimating FDR with %.0f%% threshold...", args.fdr_threshold * 100)
    df = df.sort_values("score", ascending=False).reset_index(drop=True)

    df["target_cum"] = (~df["decoy"]).cumsum()
    df["decoy_cum"] = df["decoy"].cumsum()
    df["fdr"] = df["decoy_cum"] / df["target_cum"].replace(0, 1)

    filtered_df = df[(df["fdr"] <= args.fdr_threshold) & (~df["decoy"])]
    logger.info(
        "Retained %d candidates at %.0f%% FDR (from %d targets).",
        len(filtered_df),
        args.fdr_threshold * 100,
        int((~df["decoy"]).sum()),
    )

    os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)
    out_cols = ["run_id", "spectrum_id", "peptide", "score", "fdr"]
    filtered_df[out_cols].to_csv(OUT_FILE, sep="\t", index=False)
    logger.info("Saved candidates to %s", OUT_FILE)


if __name__ == "__main__":
    main()
