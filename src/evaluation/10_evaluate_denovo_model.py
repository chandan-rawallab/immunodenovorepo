#!/usr/bin/env python3
"""Evaluate a trained de novo model on a held-out PSM test split."""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from pyteomics import mgf

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from cnnlstm.cnnlstm_model import NeoepitopeSeq2Seq
from cnnlstm.spectral_dataset import AA_TO_INT, AMINO_ACIDS


INT_TO_AA = {value: key for key, value in AA_TO_INT.items()}
PAD = AA_TO_INT["<PAD>"]
START = AA_TO_INT["<START>"]
END = AA_TO_INT["<END>"]


def encode_sequence(sequence: str, max_seq_len: int) -> torch.Tensor:
    tokens = [AA_TO_INT.get(aa, 0) for aa in list(str(sequence))]
    tokens = ([START] + tokens + [END])[:max_seq_len]
    tokens += [PAD] * (max_seq_len - len(tokens))
    return torch.tensor(tokens, dtype=torch.long)


def decode_tokens(tokens: list[int]) -> str:
    peptide = []
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
    previous = list(range(len(right) + 1))
    for i, left_char in enumerate(left, 1):
        current = [i]
        for j, right_char in enumerate(right, 1):
            current.append(min(previous[j] + 1, current[j - 1] + 1, previous[j - 1] + (left_char != right_char)))
        previous = current
    return previous[-1]


def bin_spectrum(mz_array, intensity_array, bin_size: float, max_mz: float) -> torch.Tensor:
    vector_size = int(max_mz / bin_size)
    vector = np.zeros(vector_size, dtype=np.float32)
    for mz, intensity in zip(mz_array, intensity_array):
        if mz < max_mz:
            idx = int(mz / bin_size)
            if idx < vector_size:
                vector[idx] += intensity
    max_intensity = vector.max()
    if max_intensity > 0:
        vector /= max_intensity
    return torch.tensor(vector).unsqueeze(0)


def scan_id_from_spectrum(spectrum, fallback: int) -> str:
    params = spectrum.get("params", {})
    scan = str(params.get("scans", "")).strip()
    if scan:
        return scan
    title = str(params.get("title", ""))
    if "scan=" in title:
        return title.split("scan=")[-1].strip().split()[0]
    return str(fallback)


def update_metrics(metrics: dict, pred_tokens: torch.Tensor, target_tokens: torch.Tensor, examples: list, meta: tuple[str, str]):
    pred_seq = decode_tokens(pred_tokens.tolist())
    target_seq = decode_tokens(target_tokens.tolist())
    metrics["total"] += 1
    metrics["exact"] += int(pred_seq == target_seq)
    metrics["length_ok"] += int(len(pred_seq) == len(target_seq))

    distance = edit_distance(pred_seq, target_seq)
    metrics["edit_sum"] += distance
    metrics["edit_le1"] += int(distance <= 1)

    compared = min(len(pred_seq), len(target_seq))
    metrics["aa_correct"] += sum(1 for i in range(compared) if pred_seq[i] == target_seq[i])
    metrics["aa_total"] += max(len(target_seq), 1)

    if len(examples) < 20 and pred_seq != target_seq:
        examples.append(
            {
                "run_id": meta[0],
                "scan_id": meta[1],
                "target": target_seq,
                "predicted": pred_seq,
                "edit_distance": distance,
            }
        )


def main() -> int:
    parser = argparse.ArgumentParser(description="Compute held-out de novo model accuracy metrics.")
    parser.add_argument("--checkpoint", type=Path, default=Path("results/checkpoints/neoepitope_production_best.pth"))
    parser.add_argument("--test-psms", type=Path, default=Path("results/checkpoints/test_set_psms.tsv"))
    parser.add_argument("--mgf-dir", type=Path, default=Path("data/mgf"))
    parser.add_argument("--output", type=Path, default=Path("results/model_accuracy_report.json"))
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--bin-size", type=float, default=0.1)
    parser.add_argument("--max-mz", type=float, default=2000.0)
    parser.add_argument("--max-seq-len", type=int, default=30)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--max-runs", type=int, default=None, help="Optional smoke-test limit on number of runs to evaluate.")
    args = parser.parse_args()

    if not args.checkpoint.exists():
        print(f"ERROR: checkpoint not found: {args.checkpoint}")
        return 1
    if not args.test_psms.exists():
        print(f"ERROR: test split not found: {args.test_psms}")
        return 1

    df = pd.read_csv(args.test_psms, sep="\t")
    device = torch.device(args.device if args.device == "cuda" and torch.cuda.is_available() else "cpu")
    model = NeoepitopeSeq2Seq().to(device)
    model.load_state_dict(torch.load(args.checkpoint, map_location=device))
    model.eval()

    metrics = {
        "total": 0,
        "exact": 0,
        "length_ok": 0,
        "aa_correct": 0,
        "aa_total": 0,
        "token_correct": 0,
        "token_total": 0,
        "edit_sum": 0,
        "edit_le1": 0,
        "missing_spectra": 0,
    }
    examples = []
    batch_x = []
    batch_y = []
    batch_meta = []

    def flush_batch():
        if not batch_x:
            return
        x = torch.stack(batch_x).to(device)
        y = torch.stack(batch_y)
        with torch.no_grad():
            pred = model(x).argmax(dim=-1).cpu()
        mask = y != PAD
        metrics["token_correct"] += int(((pred == y) & mask).sum().item())
        metrics["token_total"] += int(mask.sum().item())
        for pred_tokens, target_tokens, meta in zip(pred, y, batch_meta):
            update_metrics(metrics, pred_tokens, target_tokens, examples, meta)
        batch_x.clear()
        batch_y.clear()
        batch_meta.clear()

    start = time.time()
    run_count = df["run_id"].nunique()
    grouped_runs = list(df.groupby("run_id", sort=True))
    if args.max_runs is not None:
        grouped_runs = grouped_runs[: args.max_runs]

    for run_index, (run_id, group) in enumerate(grouped_runs, 1):
        wanted = {str(scan): peptide for scan, peptide in zip(group["spectrum_id"], group["peptide"])}
        wanted_scans = set(wanted)
        found = set()
        mgf_path = args.mgf_dir / f"{run_id}.mgf"
        if not mgf_path.exists():
            metrics["missing_spectra"] += len(group)
            print(f"missing MGF: {run_id} ({len(group)} rows)")
            continue

        for idx, spectrum in enumerate(mgf.read(str(mgf_path))):
            scan = scan_id_from_spectrum(spectrum, idx)
            if scan not in wanted_scans:
                continue
            found.add(scan)
            batch_x.append(bin_spectrum(spectrum.get("m/z array", []), spectrum.get("intensity array", []), args.bin_size, args.max_mz))
            batch_y.append(encode_sequence(wanted[scan], args.max_seq_len))
            batch_meta.append((run_id, scan))
            if len(batch_x) >= args.batch_size:
                flush_batch()
            if len(found) == len(wanted_scans):
                break

        missing = len(wanted_scans - found)
        metrics["missing_spectra"] += missing
        print(f"run {run_index:02d}/{len(grouped_runs)} {run_id}: found {len(found)}/{len(wanted_scans)} missing {missing}")

    flush_batch()
    elapsed = time.time() - start

    total = metrics["total"]
    report = {
        "checkpoint": str(args.checkpoint),
        "test_psms": str(args.test_psms),
        "mgf_dir": str(args.mgf_dir),
        "device": str(device),
        "test_rows": int(len(df)),
        "evaluated_pairs": int(total),
        "missing_spectra": int(metrics["missing_spectra"]),
        "exact_peptide_accuracy": metrics["exact"] / total if total else 0.0,
        "length_accuracy": metrics["length_ok"] / total if total else 0.0,
        "token_accuracy_excluding_pad": metrics["token_correct"] / metrics["token_total"] if metrics["token_total"] else 0.0,
        "position_aa_accuracy_len_normalized": metrics["aa_correct"] / metrics["aa_total"] if metrics["aa_total"] else 0.0,
        "mean_edit_distance": metrics["edit_sum"] / total if total else 0.0,
        "edit_distance_le_1_rate": metrics["edit_le1"] / total if total else 0.0,
        "counts": metrics,
        "elapsed_seconds": elapsed,
        "first_mismatches": examples,
        "interpretation": (
            "Exact peptide accuracy below 30% indicates the checkpoint should be treated as a weak baseline and retrained "
            "after manifest/reference cleanup."
        ),
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2))

    print("\n=== TRAINED MODEL HELD-OUT TEST ACCURACY ===")
    for key in [
        "evaluated_pairs",
        "missing_spectra",
        "exact_peptide_accuracy",
        "length_accuracy",
        "token_accuracy_excluding_pad",
        "position_aa_accuracy_len_normalized",
        "mean_edit_distance",
        "edit_distance_le_1_rate",
    ]:
        print(f"{key}: {report[key]}")
    print(f"Saved report to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
