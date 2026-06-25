"""Train the Objective 3 CNN-LSTM de novo sequencing model.

Improvements applied (2026-06-25) to reach ≥ 30% exact-peptide accuracy:
  - CosineAnnealingWarmRestarts replaces the step-decay LR schedule.
    This allows the optimiser to escape flat regions and converges more reliably
    on sequence-to-sequence tasks (T_0=10, T_mult=2).
  - Label smoothing (ε = 0.10) added to CrossEntropyLoss.
    Reduces overconfidence and improves calibration for the output softmax.
  - Patience increased from 10 → 25 epochs to let CosineAnnealing complete
    at least two restart cycles before early-stopping.
  - Mixed-precision (AMP) training enabled when CUDA is available via
    torch.amp.autocast — roughly 2× faster and uses less VRAM.
  - Optional precursor mass/charge tensors are now assembled per-batch from
    the PSM file and forwarded to the model when the column is present.
  - Default LR lowered to 3 e-4 (better starting point for AdamW with
    CosineAnnealing on this architecture).
  - Gradient-accumulation parameter added (default: 2 steps) to simulate a
    larger effective batch size without increasing VRAM usage.

Audit fixes carried forward (2026-06-20):
  - Train/test split is peptide-grouped (GroupShuffleSplit) to prevent
    identical peptide sequences from appearing in both train and test.
  - Training configuration and manifest provenance are saved alongside
    model checkpoints as JSON metadata files.
  - audit_summary() from SpectralDataset is reported before training starts.
"""

import json
import logging
import os
import sys
from datetime import datetime, timezone

import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Subset

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# =========================================================================
# DYNAMIC PATH RESOLUTION (FOR ACADEMIC STRUCTURE)
# =========================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.append(PROJECT_ROOT)

from cnnlstm.cnnlstm_model import NeoepitopeSeq2Seq
from cnnlstm.spectral_dataset import SpectralDataset, VOCAB_SIZE, AMINO_ACIDS


# =========================================================================
# Peptide-grouped split — prevents train/test leakage
# =========================================================================

def _peptide_grouped_split(
    psm_df, train_frac: float = 0.70, val_frac: float = 0.15, random_state: int = 42
) -> tuple[list[int], list[int], list[int]]:
    """Return train/val/test index lists split by unique peptide identity.

    This prevents spectra of the same peptide sequence from leaking across
    train/val/test partitions, which would inflate held-out accuracy.
    """
    import numpy as np

    rng = np.random.default_rng(random_state)
    unique_peptides = psm_df["peptide"].unique()
    rng.shuffle(unique_peptides)

    n = len(unique_peptides)
    n_train = int(n * train_frac)
    n_val = int(n * val_frac)

    train_peps = set(unique_peptides[:n_train])
    val_peps = set(unique_peptides[n_train : n_train + n_val])

    train_idx, val_idx, test_idx = [], [], []
    for idx, row in psm_df.iterrows():
        pep = row["peptide"]
        if pep in train_peps:
            train_idx.append(idx)
        elif pep in val_peps:
            val_idx.append(idx)
        else:
            test_idx.append(idx)

    logger.info(
        "Peptide-grouped split → train %d | val %d | test %d PSMs "
        "(%d / %d / %d unique peptides)",
        len(train_idx), len(val_idx), len(test_idx),
        len(train_peps), len(val_peps), n - n_train - n_val,
    )
    return train_idx, val_idx, test_idx


# =========================================================================
# Checkpoint metadata sidecar
# =========================================================================

def _save_checkpoint_metadata(
    path: str,
    config: dict,
    psm_file: str,
    mgf_dir: str,
    epoch: int,
    val_loss: float,
    dataset_summary: dict,
) -> None:
    """Save JSON sidecar file alongside every model checkpoint."""
    meta = {
        "saved_at": datetime.now(timezone.utc).isoformat(),
        "epoch": epoch,
        "val_loss": val_loss,
        "model_class": "NeoepitopeSeq2Seq",
        "vocab_size": VOCAB_SIZE,
        "training_config": config,
        "provenance": {
            "psm_file": psm_file,
            "mgf_dir": mgf_dir,
        },
        "dataset_audit": dataset_summary,
    }
    meta_path = path.replace(".pth", "_metadata.json")
    with open(meta_path, "w") as fh:
        json.dump(meta, fh, indent=2)
    logger.info("Saved checkpoint metadata to %s", meta_path)


# =========================================================================
# Precursor mass/charge helpers
# =========================================================================

def _has_precursor_columns(psm_df) -> bool:
    return "precursor_mz" in psm_df.columns and "charge" in psm_df.columns


def _build_precursor_tensor(rows, device: torch.device) -> torch.Tensor | None:
    """Assemble a (B, 2) [neutral_mass, charge] tensor for a batch of PSM rows.

    neutral_mass = (precursor_mz * charge) - (charge * 1.00728)
    Falls back to None if columns are absent or values are non-numeric.
    """
    try:
        masses = (rows["precursor_mz"].values * rows["charge"].values
                  - rows["charge"].values * 1.00728).tolist()
        charges = rows["charge"].values.tolist()
        t = torch.tensor(
            [[m, c] for m, c in zip(masses, charges)], dtype=torch.float32
        )
        return t.to(device)
    except Exception:
        return None


# =========================================================================
# Training entry point
# =========================================================================

def train_production(
    mgf_dir=None,
    psm_file=None,
    checkpoint_dir=None,
    batch_size: int = 32,
    epochs: int = 150,
    lr: float = 3e-4,
    num_workers: int = 4,
    t0: int = 10,
    t_mult: int = 2,
    label_smoothing: float = 0.10,
    patience: int = 25,
    grad_accum_steps: int = 2,
    dry_run: bool = False,          # CPU sanity-check mode: 2 epochs, 64 samples
):
    # ── 1. Configuration ─────────────────────────────────────────────────────
    WORKSPACE_ROOT = os.path.dirname(PROJECT_ROOT)
    MGF_DIR = mgf_dir or os.path.join(WORKSPACE_ROOT, "data", "mgf")
    PSM_FILE = psm_file or os.path.join(
        WORKSPACE_ROOT, "results", "immunopeptidome_psms.tsv"
    )
    CHECKPOINT_DIR = checkpoint_dir or os.path.join(
        WORKSPACE_ROOT, "results", "checkpoints"
    )

    os.makedirs(CHECKPOINT_DIR, exist_ok=True)

    DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    USE_AMP = DEVICE.type == "cuda"
    scaler = torch.amp.GradScaler(enabled=USE_AMP)

    if dry_run:
        logger.warning(
            "DRY-RUN mode: capping to 2 epochs, 64 samples, batch_size=8, "
            "num_workers=0. This is for pipeline verification only, NOT real training."
        )
        epochs = 2
        batch_size = 8
        num_workers = 0
    elif DEVICE.type == "cpu":
        logger.warning(
            "No CUDA GPU detected — training on CPU. "
            "Full training will be extremely slow. "
            "Use --dry-run to do a quick pipeline check, or run on Google Colab "
            "(see notebooks/colab_train.ipynb) for GPU-accelerated training."
        )
        num_workers = 0  # multiprocessing on CPU can deadlock

    config = {
        "batch_size": batch_size,
        "epochs": epochs,
        "lr": lr,
        "t0": t0,
        "t_mult": t_mult,
        "label_smoothing": label_smoothing,
        "patience": patience,
        "grad_accum_steps": grad_accum_steps,
        "num_workers": num_workers,
        "device": str(DEVICE),
        "amp": USE_AMP,
    }

    logger.info("Initialising production training on %s  (AMP=%s)", DEVICE, USE_AMP)
    logger.info("Spectra : %s", MGF_DIR)
    logger.info("PSMs    : %s", PSM_FILE)
    logger.info("Config  : %s", json.dumps(config))

    # ── 2. Load dataset ───────────────────────────────────────────────────────
    full_dataset = SpectralDataset(PSM_FILE, MGF_DIR, bin_size=0.1)

    if len(full_dataset) == 0:
        logger.error(
            "No labelled spectra found. Has 04_extract_psms.py finished yet? "
            "MGF dir: %s  PSM file: %s",
            MGF_DIR, PSM_FILE,
        )
        return

    audit = full_dataset.audit_summary()
    logger.info("Dataset audit: %s", json.dumps(audit))

    # ── 3. Peptide-grouped train/val/test split ───────────────────────────────
    train_idx, val_idx, test_idx = _peptide_grouped_split(full_dataset.psms)

    train_dataset = Subset(full_dataset, train_idx)
    val_dataset   = Subset(full_dataset, val_idx)
    test_dataset  = Subset(full_dataset, test_idx)  # noqa: F841

    # Dry-run: cap to first 64 samples only
    if dry_run:
        from torch.utils.data import Subset as _Sub
        cap = min(64, len(train_idx))
        train_dataset = _Sub(full_dataset, train_idx[:cap])
        val_dataset   = _Sub(full_dataset, val_idx[:min(16, len(val_idx))])

    # Persist test set for locked evaluation
    test_psms = full_dataset.psms.iloc[test_idx]
    test_psms_path = os.path.join(CHECKPOINT_DIR, "test_set_psms.tsv")
    test_psms.to_csv(test_psms_path, sep="\t", index=False)
    logger.info("Test-set PSMs saved to %s (%d rows)", test_psms_path, len(test_psms))

    # Leakage audit
    train_peps = set(full_dataset.psms.iloc[train_idx]["peptide"])
    test_peps  = set(full_dataset.psms.iloc[test_idx]["peptide"])
    overlap = train_peps & test_peps
    if overlap:
        logger.warning(
            "LEAKAGE WARNING: %d peptides appear in both train and test splits.",
            len(overlap),
        )
    else:
        logger.info("Peptide overlap between train and test: 0 (clean split).")

    # Check for precursor columns
    has_prec = _has_precursor_columns(full_dataset.psms)
    logger.info("Precursor mass/charge conditioning: %s", has_prec)

    train_loader = DataLoader(
        train_dataset, batch_size=batch_size, shuffle=True,
        num_workers=num_workers, pin_memory=(DEVICE.type == "cuda"),
    )
    val_loader = DataLoader(
        val_dataset, batch_size=batch_size, shuffle=False,
        num_workers=num_workers, pin_memory=(DEVICE.type == "cuda"),
    )

    # ── 4. Model, loss, optimiser, scheduler ─────────────────────────────────
    model = NeoepitopeSeq2Seq().to(DEVICE)

    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    logger.info("Model parameters: %d (%.1f M)", n_params, n_params / 1e6)

    # Label smoothing reduces overconfidence on ambiguous spectra
    criterion = nn.CrossEntropyLoss(
        ignore_index=0,
        label_smoothing=label_smoothing,
    )

    optimizer = torch.optim.AdamW(
        model.parameters(), lr=lr, weight_decay=1e-4, eps=1e-8
    )

    # CosineAnnealingWarmRestarts: period doubles after each restart
    scheduler = torch.optim.lr_scheduler.CosineAnnealingWarmRestarts(
        optimizer, T_0=t0, T_mult=t_mult, eta_min=lr * 1e-3
    )

    logger.info(
        "Training: %d train | %d val | %d test PSMs",
        len(train_dataset), len(val_dataset), len(test_dataset),
    )

    best_val_loss = float("inf")
    trigger_times = 0

    # ── 5. Training loop ──────────────────────────────────────────────────────
    for epoch in range(1, epochs + 1):
        model.train()
        total_loss = 0.0
        optimizer.zero_grad()

        for batch_idx, (spectra, sequences) in enumerate(train_loader):
            spectra, sequences = spectra.to(DEVICE), sequences.to(DEVICE)

            with torch.amp.autocast(device_type=DEVICE.type, enabled=USE_AMP):
                # Precursor conditioning: pull from the PSM subset for this batch
                precursor = None
                if has_prec:
                    start = batch_idx * batch_size
                    end   = start + spectra.size(0)
                    rows  = full_dataset.psms.iloc[
                        [train_idx[i] for i in range(
                            min(start, len(train_idx)),
                            min(end,   len(train_idx))
                        )]
                    ]
                    precursor = _build_precursor_tensor(rows, DEVICE)

                outputs = model(spectra, precursor=precursor)
                outputs = outputs.view(-1, VOCAB_SIZE)
                sequences_flat = sequences.view(-1)
                loss = criterion(outputs, sequences_flat) / grad_accum_steps

            scaler.scale(loss).backward()

            if (batch_idx + 1) % grad_accum_steps == 0 or (
                batch_idx + 1 == len(train_loader)
            ):
                scaler.unscale_(optimizer)
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                scaler.step(optimizer)
                scaler.update()
                optimizer.zero_grad()

            total_loss += loss.item() * grad_accum_steps  # undo the division

            if batch_idx % 100 == 0:
                logger.info(
                    "Epoch [%d/%d] | Batch %d/%d | Loss: %.4f | LR: %.2e",
                    epoch, epochs, batch_idx, len(train_loader),
                    loss.item() * grad_accum_steps,
                    optimizer.param_groups[0]["lr"],
                )

        avg_train_loss = total_loss / len(train_loader)

        # ── Validation ────────────────────────────────────────────────────────
        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for spectra, sequences in val_loader:
                spectra, sequences = spectra.to(DEVICE), sequences.to(DEVICE)
                with torch.amp.autocast(device_type=DEVICE.type, enabled=USE_AMP):
                    outputs = model(spectra).view(-1, VOCAB_SIZE)
                    val_loss += criterion(outputs, sequences.view(-1)).item()

        avg_val_loss = val_loss / len(val_loader)

        # Step scheduler each epoch (CosineAnnealingWarmRestarts uses epoch)
        scheduler.step(epoch)

        # ── Sample prediction for visual sanity check ─────────────────────────
        with torch.no_grad():
            sample_spectra, _ = next(iter(val_loader))
            sample_spectra = sample_spectra[0].unsqueeze(0).to(DEVICE)
            pred = model(sample_spectra).argmax(dim=-1).squeeze().tolist()
            pred_chars = []
            for t in pred:
                if t == 0:
                    continue
                if t == 21:
                    pred_chars.append("<S>")
                elif t == 22:
                    pred_chars.append("<E>")
                elif 0 < t <= len(AMINO_ACIDS):
                    pred_chars.append(AMINO_ACIDS[t - 1])
                else:
                    pred_chars.append("?")
            logger.info("Sample prediction (epoch %d): %s", epoch, "".join(pred_chars))

        logger.info(
            "==> Epoch %d | Train %.4f | Val %.4f | LR %.2e",
            epoch, avg_train_loss, avg_val_loss,
            optimizer.param_groups[0]["lr"],
        )

        # ── Checkpointing and early stopping ─────────────────────────────────
        if avg_val_loss < best_val_loss:
            best_val_loss = avg_val_loss
            trigger_times = 0
            best_path = os.path.join(CHECKPOINT_DIR, "neoepitope_production_best.pth")
            torch.save(model.state_dict(), best_path)
            _save_checkpoint_metadata(
                best_path, config, PSM_FILE, MGF_DIR, epoch, avg_val_loss, audit
            )
            logger.info("Saved BEST model (val_loss=%.4f) → %s", avg_val_loss, best_path)
        else:
            trigger_times += 1
            if trigger_times >= patience:
                logger.info(
                    "Early stopping triggered at epoch %d "
                    "(no improvement for %d epochs).",
                    epoch, patience,
                )
                break

        if epoch % 10 == 0:
            ckpt_path = os.path.join(
                CHECKPOINT_DIR, f"neoepitope_production_epoch_{epoch}.pth"
            )
            torch.save(model.state_dict(), ckpt_path)
            _save_checkpoint_metadata(
                ckpt_path, config, PSM_FILE, MGF_DIR, epoch, avg_val_loss, audit
            )
            logger.info("Saved periodic checkpoint → %s", ckpt_path)

    # ── Final save ────────────────────────────────────────────────────────────
    final_path = os.path.join(CHECKPOINT_DIR, "neoepitope_production_final.pth")
    torch.save(model.state_dict(), final_path)
    _save_checkpoint_metadata(
        final_path, config, PSM_FILE, MGF_DIR, epoch, avg_val_loss, audit
    )
    logger.info("FINAL model saved → %s", final_path)


# =========================================================================
# CLI
# =========================================================================
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Train the Objective 3 CNN/LSTM model on production MGF + MaxQuant labels."
    )
    parser.add_argument("--mgf-dir", default="/home/amity/hla_data_mgf")
    parser.add_argument(
        "--psm-file",
        default=None,
        help="Path to immunopeptidome_psms.tsv. Defaults to results/immunopeptidome_psms.tsv.",
    )
    parser.add_argument(
        "--checkpoint-dir",
        default=os.path.join(PROJECT_ROOT, "results", "checkpoints_v3"),
    )
    parser.add_argument("--batch-size",      type=int,   default=32)
    parser.add_argument("--epochs",          type=int,   default=200)
    parser.add_argument("--lr",              type=float, default=3e-4)
    parser.add_argument("--num-workers",     type=int,   default=0)
    parser.add_argument("--t0",              type=int,   default=10,
                        help="CosineAnnealingWarmRestarts T_0 (first restart epoch).")
    parser.add_argument("--t-mult",          type=int,   default=2,
                        help="CosineAnnealingWarmRestarts T_mult.")
    parser.add_argument("--label-smoothing", type=float, default=0.10)
    parser.add_argument("--patience",        type=int,   default=25)
    parser.add_argument("--grad-accum",      type=int,   default=2,
                        help="Gradient accumulation steps (effective batch = batch-size × grad-accum).")
    parser.add_argument("--dry-run",         action="store_true", default=False,
                        help="CPU pipeline check: 2 epochs, 64 samples. No GPU needed.")
    args = parser.parse_args()

    train_production(
        mgf_dir=args.mgf_dir,
        psm_file=args.psm_file,
        checkpoint_dir=args.checkpoint_dir,
        batch_size=args.batch_size,
        epochs=args.epochs,
        lr=args.lr,
        num_workers=args.num_workers,
        t0=args.t0,
        t_mult=args.t_mult,
        label_smoothing=args.label_smoothing,
        patience=args.patience,
        grad_accum_steps=args.grad_accum,
        dry_run=args.dry_run,
    )
