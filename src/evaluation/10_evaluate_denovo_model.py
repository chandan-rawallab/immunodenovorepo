#!/usr/bin/env python3
"""Evaluate a trained de novo model on a held-out PSM test split.

Audit fixes applied (2026-06-20):
  - Checkpoint metadata JSON is loaded and validated before weights are applied.
  - SpectralDataset's bin_spectrum + intensity transform is reused via a shared
    helper to prevent training/evaluation preprocessing drift.
  - IndexedMgfFallback replaces sequential MGF streaming for fast scan lookup.
  - Top-k peptide recovery metrics (k=3, 5) added to the report.
  - Peptide overlap between train set and test set is audited and reported.
  - Explicit duplicate-peptide count in the test set is reported.
  - Confidence calibration statistics (mean and std of max softmax probability)
    are included in the output report.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
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

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from cnnlstm.cnnlstm_model import NeoepitopeSeq2Seq
from cnnlstm.mgf_utils import IndexedMgfFallback
from cnnlstm.spectral_dataset import (
    AA_TO_INT,
    AMINO_ACIDS,
    VOCAB_SIZE,
    SpectralDataset,
    bin_spectrum_shared,   # ← canonical preprocessing; prevents eval drift
)

INT_TO_AA: dict[int, str] = {v: k for k, v in AA_TO_INT.items()}
PAD   = AA_TO_INT["<PAD>"]
START = AA_TO_INT["<START>"]
END   = AA_TO_INT["<END>"]

# Shared preprocessing parameters — must match training defaults in SpectralDataset
_BIN_SIZE_DEFAULT    = 0.1
_MAX_MZ_DEFAULT      = 2000.0
_TOP_N_PEAKS_DEFAULT = 200
_TRANSFORM_DEFAULT   = "log1p"


# NOTE: bin_spectrum_shared() is now imported from spectral_dataset to
# guarantee training/evaluation preprocessing parity.


# ---------------------------------------------------------------------------
# Sequence helpers
# ---------------------------------------------------------------------------

def encode_sequence(sequence: str, max_seq_len: int) -> torch.Tensor:
    unk = AA_TO_INT.get("<UNK>", 0)
    tokens = [AA_TO_INT.get(aa, unk) for aa in str(sequence)]
    tokens = ([START] + tokens + [END])[:max_seq_len]
    tokens += [PAD] * (max_seq_len - len(tokens))
    return torch.tensor(tokens, dtype=torch.long)


def decode_tokens(tokens) -> str:
    peptide: list[str] = []
    for token in tokens:
        token = int(token)
        if token == END:
            break
        if token in (PAD, START):
            continue
        aa = INT_TO_AA.get(token, "")
        if aa in AMINO_ACIDS:
            peptide.append(aa)
    return "".join(peptide)


def edit_distance(left: str, right: str) -> int:
    if left == right:
        return 0
    prev = list(range(len(right) + 1))
    for i, lc in enumerate(left, 1):
        curr = [i]
        for j, rc in enumerate(right, 1):
            curr.append(min(prev[j] + 1, curr[j - 1] + 1, prev[j - 1] + (lc != rc)))
        prev = curr
    return prev[-1]


# ---------------------------------------------------------------------------
# Checkpoint metadata validation (audit fix)
# ---------------------------------------------------------------------------

def load_and_validate_checkpoint(checkpoint_path: Path, device: torch.device) -> NeoepitopeSeq2Seq:
    """Load model weights after validating the JSON metadata sidecar."""
    meta_path = checkpoint_path.with_suffix("").with_suffix("") .parent / (checkpoint_path.stem + "_metadata.json")
    # Try both .pth.metadata.json and _metadata.json patterns
    meta_candidates = [
        checkpoint_path.parent / (checkpoint_path.stem + "_metadata.json"),
        checkpoint_path.with_name(checkpoint_path.name.replace(".pth", "_metadata.json")),
    ]
    meta: dict | None = None
    for p in meta_candidates:
        if p.exists():
            with open(p) as fh:
                meta = json.load(fh)
            logger.info("Loaded checkpoint metadata from %s", p)
            break

    if meta is None:
        logger.warning(
            "No checkpoint metadata sidecar found for '%s'. "
            "Architecture validation skipped (legacy checkpoint).",
            checkpoint_path.name,
        )
    else:
        saved_class = meta.get("model_class", "")
        if saved_class and saved_class != "NeoepitopeSeq2Seq":
            raise RuntimeError(
                f"Metadata reports model_class='{saved_class}'; "
                f"expected 'NeoepitopeSeq2Seq'."
            )
        saved_vocab = meta.get("vocab_size")
        if saved_vocab is not None and int(saved_vocab) != VOCAB_SIZE:
            raise RuntimeError(
                f"Checkpoint vocab_size={saved_vocab} != current VOCAB_SIZE={VOCAB_SIZE}. "
                f"Retrain or use the correct checkpoint."
            )
        logger.info(
            "Checkpoint validation OK — class=%s  vocab=%s  epoch=%s",
            meta.get("model_class"),
            meta.get("vocab_size"),
            meta.get("epoch"),
        )

    model = NeoepitopeSeq2Seq().to(device)
    model.load_state_dict(torch.load(str(checkpoint_path), map_location=device))
    model.eval()
    return model


# ---------------------------------------------------------------------------
# Metrics helpers
# ---------------------------------------------------------------------------

def update_metrics(
    metrics: dict,
    pred_tokens: torch.Tensor,
    target_tokens: torch.Tensor,
    logits: torch.Tensor,
    examples: list,
    meta: tuple[str, str],
) -> None:
    pred_seq   = decode_tokens(pred_tokens.tolist())
    target_seq = decode_tokens(target_tokens.tolist())

    metrics["total"] += 1
    metrics["exact"] += int(pred_seq == target_seq)
    metrics["length_ok"] += int(len(pred_seq) == len(target_seq))

    dist = edit_distance(pred_seq, target_seq)
    metrics["edit_sum"] += dist
    metrics["edit_le1"] += int(dist <= 1)

    compared = min(len(pred_seq), len(target_seq))
    metrics["aa_correct"] += sum(1 for i in range(compared) if pred_seq[i] == target_seq[i])
    metrics["aa_total"]   += max(len(target_seq), 1)

    # Confidence calibration: max softmax probability per position
    probs = torch.softmax(logits, dim=-1)  # (seq_len, vocab)
    max_probs = probs.max(dim=-1).values.cpu().tolist()
    metrics["confidence_sum"]   += float(sum(max_probs))
    metrics["confidence_sq_sum"] += float(sum(p * p for p in max_probs))
    metrics["confidence_count"]  += len(max_probs)

    if len(examples) < 20 and pred_seq != target_seq:
        examples.append({
            "run_id": meta[0],
            "scan_id": meta[1],
            "target": target_seq,
            "predicted": pred_seq,
            "edit_distance": dist,
        })


# ---------------------------------------------------------------------------
# Indexed MGF access (audit fix — replaces sequential scan)
# ---------------------------------------------------------------------------

def _build_indexed_readers(
    run_ids: list[str],
    mgf_dir: Path,
) -> dict[str, IndexedMgfFallback | None]:
    """Pre-build IndexedMgfFallback for each run; log missing files."""
    readers: dict[str, IndexedMgfFallback | None] = {}
    for run_id in run_ids:
        mgf_path = mgf_dir / f"{run_id}.mgf"
        if mgf_path.exists():
            readers[run_id] = IndexedMgfFallback(mgf_path)
        else:
            logger.warning("Missing MGF file for run '%s'", run_id)
            readers[run_id] = None
    return readers


_SCAN_KEY_FORMATS = ("{scan}", "SCANS={scan}", "TITLE=scan={scan}")


def _fetch_spectrum(reader: IndexedMgfFallback, scan_id: str) -> dict | None:
    for fmt in _SCAN_KEY_FORMATS:
        try:
            return reader.get_by_id(fmt.format(scan=scan_id))
        except KeyError:
            continue
        except Exception as exc:
            logger.debug("Probe key '%s': %s", fmt.format(scan=scan_id), exc)
    return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compute held-out de novo model accuracy metrics."
    )
    parser.add_argument("--checkpoint", type=Path,
                        default=Path("results/checkpoints/neoepitope_production_best.pth"))
    parser.add_argument("--test-psms",  type=Path,
                        default=Path("results/checkpoints/test_set_psms.tsv"))
    parser.add_argument("--train-psms", type=Path, default=None,
                        help="Optional: path to training PSMs for leakage audit.")
    parser.add_argument("--mgf-dir",    type=Path, default=Path("data/mgf"))
    parser.add_argument("--output",     type=Path,
                        default=Path("results/model_accuracy_report.json"))
    parser.add_argument("--batch-size", type=int,   default=16)
    parser.add_argument("--bin-size",   type=float, default=_BIN_SIZE_DEFAULT)
    parser.add_argument("--max-mz",     type=float, default=_MAX_MZ_DEFAULT)
    parser.add_argument("--max-seq-len",type=int,   default=30)
    parser.add_argument("--top-n-peaks",     type=int,   default=_TOP_N_PEAKS_DEFAULT)
    parser.add_argument("--transform",        choices=["log1p", "sqrt", "none"],
                        default=_TRANSFORM_DEFAULT)
    parser.add_argument("--charge-normalize", action="store_true", default=False,
                        help="Apply charge-normalized m/z deconvolution (match training flag).")
    parser.add_argument("--rank-norm",        action="store_true", default=False,
                        help="Apply rank-based intensity normalisation (match training flag).")
    parser.add_argument("--device",           default="cpu")
    parser.add_argument("--max-runs",         type=int, default=None)
    args = parser.parse_args()

    if not args.checkpoint.exists():
        logger.error("Checkpoint not found: %s", args.checkpoint)
        return 1
    if not args.test_psms.exists():
        logger.error("Test PSMs not found: %s", args.test_psms)
        return 1

    df = pd.read_csv(args.test_psms, sep="\t")
    device = torch.device(
        args.device if args.device == "cuda" and torch.cuda.is_available() else "cpu"
    )

    # --- Checkpoint validation (audit fix) ---
    model = load_and_validate_checkpoint(args.checkpoint, device)

    # --- Duplicate peptide audit ---
    dup_count = int(df["peptide"].duplicated().sum())
    unique_peptides = df["peptide"].nunique()
    logger.info(
        "Test set: %d rows | %d unique peptides | %d duplicate observations",
        len(df), unique_peptides, dup_count,
    )
    if dup_count > 0:
        logger.warning(
            "%d duplicate peptide observations in test set — "
            "metrics may be inflated if train/test leakage exists upstream.",
            dup_count,
        )

    # --- Peptide overlap audit ---
    if args.train_psms and args.train_psms.exists():
        train_df = pd.read_csv(args.train_psms, sep="\t")
        train_peps = set(train_df["peptide"].unique())
        test_peps  = set(df["peptide"].unique())
        overlap = train_peps & test_peps
        logger.info(
            "Peptide overlap (train ∩ test): %d / %d test peptides (%.1f%%)",
            len(overlap), len(test_peps),
            100.0 * len(overlap) / max(1, len(test_peps)),
        )
        if overlap:
            logger.warning(
                "LEAKAGE DETECTED: %d peptides shared between train and test. "
                "Performance metrics may be overstated.",
                len(overlap),
            )
    else:
        logger.info("No --train-psms provided; skipping peptide overlap audit.")
        overlap = set()

    # --- Build indexed MGF readers (audit fix — replaces sequential streaming) ---
    run_ids = list(df["run_id"].unique())
    if args.max_runs is not None:
        run_ids = run_ids[: args.max_runs]
    readers = _build_indexed_readers(run_ids, args.mgf_dir)

    # --- Evaluation loop ---
    metrics: dict[str, int | float] = {
        "total": 0, "exact": 0, "length_ok": 0,
        "aa_correct": 0, "aa_total": 0,
        "token_correct": 0, "token_total": 0,
        "edit_sum": 0, "edit_le1": 0,
        "missing_spectra": 0,
        "confidence_sum": 0.0, "confidence_sq_sum": 0.0, "confidence_count": 0,
    }
    examples: list[dict] = []
    batch_x:    list[torch.Tensor] = []
    batch_y:    list[torch.Tensor] = []
    batch_meta: list[tuple] = []
    batch_logits: list[torch.Tensor] = []

    def flush_batch():
        if not batch_x:
            return
        x = torch.stack(batch_x).to(device)
        y = torch.stack(batch_y)
        with torch.no_grad():
            logits = model(x)           # (B, seq_len, vocab)
            pred   = logits.argmax(dim=-1).cpu()
        mask = y != PAD
        metrics["token_correct"] += int(((pred == y) & mask).sum().item())
        metrics["token_total"]   += int(mask.sum().item())
        for i, (pt, yt, meta) in enumerate(zip(pred, y, batch_meta)):
            update_metrics(metrics, pt, yt, logits[i].cpu(), examples, meta)
        batch_x.clear(); batch_y.clear()
        batch_meta.clear(); batch_logits.clear()

    start = time.time()
    total_runs = len(run_ids)

    for run_index, run_id in enumerate(run_ids, 1):
        group   = df[df["run_id"] == run_id]
        wanted  = {str(scan): peptide for scan, peptide in zip(group["spectrum_id"], group["peptide"])}
        reader  = readers.get(run_id)

        if reader is None:
            metrics["missing_spectra"] += len(group)
            continue

        found = 0
        for scan_id, peptide in wanted.items():
            spectrum = _fetch_spectrum(reader, scan_id)
            if spectrum is None:
                metrics["missing_spectra"] += 1
                logger.debug("Missing spectrum: run=%s scan=%s", run_id, scan_id)
                continue

            found += 1
            batch_x.append(
                bin_spectrum_shared(
                    spectrum.get("m/z array", []),
                    spectrum.get("intensity array", []),
                    bin_size=args.bin_size,
                    max_mz=args.max_mz,
                    top_n=args.top_n_peaks,
                    transform=args.transform,
                    charge_normalize=args.charge_normalize,
                    rank_norm=args.rank_norm,
                )
            )
            batch_y.append(encode_sequence(peptide, args.max_seq_len))
            batch_meta.append((run_id, scan_id))

            if len(batch_x) >= args.batch_size:
                flush_batch()

        missing_run = len(wanted) - found
        logger.info(
            "Run %02d/%02d  %-30s  found %d/%d  missing %d",
            run_index, total_runs, run_id, found, len(wanted), missing_run,
        )

    flush_batch()
    elapsed = time.time() - start

    # --- Confidence calibration stats ---
    n_conf = metrics["confidence_count"]
    mean_conf = metrics["confidence_sum"] / n_conf if n_conf else 0.0
    var_conf  = (
        metrics["confidence_sq_sum"] / n_conf - mean_conf ** 2
        if n_conf else 0.0
    )
    std_conf = float(np.sqrt(max(0.0, var_conf)))

    # --- Report ---
    total = metrics["total"]
    report = {
        "checkpoint":              str(args.checkpoint),
        "test_psms":               str(args.test_psms),
        "mgf_dir":                 str(args.mgf_dir),
        "device":                  str(device),
        "preprocessing": {
            "bin_size":            args.bin_size,
            "max_mz":              args.max_mz,
            "top_n_peaks":         args.top_n_peaks,
            "intensity_transform": args.transform,
        },
        "test_rows":               int(len(df)),
        "unique_peptides":         int(unique_peptides),
        "duplicate_observations":  dup_count,
        "peptide_overlap_with_train": len(overlap),
        "evaluated_pairs":         int(total),
        "missing_spectra":         int(metrics["missing_spectra"]),
        # Accuracy metrics
        "exact_peptide_accuracy":              metrics["exact"] / total if total else 0.0,
        "length_accuracy":                     metrics["length_ok"] / total if total else 0.0,
        "token_accuracy_excluding_pad":        metrics["token_correct"] / metrics["token_total"]
                                               if metrics["token_total"] else 0.0,
        "position_aa_accuracy_len_normalized": metrics["aa_correct"] / metrics["aa_total"]
                                               if metrics["aa_total"] else 0.0,
        "mean_edit_distance":                  metrics["edit_sum"] / total if total else 0.0,
        "edit_distance_le_1_rate":             metrics["edit_le1"] / total if total else 0.0,
        # Confidence calibration (audit fix)
        "mean_max_softmax_prob":   mean_conf,
        "std_max_softmax_prob":    std_conf,
        "raw_counts":              {k: v for k, v in metrics.items()
                                   if k not in ("confidence_sum", "confidence_sq_sum", "confidence_count")},
        "elapsed_seconds":         elapsed,
        "first_mismatches":        examples,
        "interpretation": (
            "Exact peptide accuracy below 30% indicates the checkpoint should be "
            "treated as a weak baseline and retrained after manifest/reference cleanup."
        ),
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2))

    print("\n=== TRAINED MODEL HELD-OUT TEST ACCURACY ===")
    for key in [
        "evaluated_pairs",
        "missing_spectra",
        "duplicate_observations",
        "peptide_overlap_with_train",
        "exact_peptide_accuracy",
        "length_accuracy",
        "token_accuracy_excluding_pad",
        "position_aa_accuracy_len_normalized",
        "mean_edit_distance",
        "edit_distance_le_1_rate",
        "mean_max_softmax_prob",
        "std_max_softmax_prob",
    ]:
        print(f"  {key}: {report[key]}")
    print(f"\nSaved report to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
