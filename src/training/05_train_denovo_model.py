"""Train the Objective 3 CNN-LSTM de novo sequencing model.

Audit fixes applied (2026-06-20):
  - Train/test split is now peptide-grouped (GroupShuffleSplit) to prevent
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
    # Remaining peptides go to test

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
        len(train_idx),
        len(val_idx),
        len(test_idx),
        len(train_peps),
        len(val_peps),
        n - n_train - n_val,
    )
    return train_idx, val_idx, test_idx


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


def train_production(
    mgf_dir=None,
    psm_file=None,
    checkpoint_dir=None,
    batch_size=32,
    epochs=150,
    lr=0.001,
    num_workers=4,
    warmup_epochs=5,
):
    # 1. Configuration (relative to project structure)
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

    config = {
        "batch_size": batch_size,
        "epochs": epochs,
        "lr": lr,
        "num_workers": num_workers,
        "warmup_epochs": warmup_epochs,
        "device": str(DEVICE),
    }

    logger.info("Initialising production training on %s", DEVICE)
    logger.info("Spectra : %s", MGF_DIR)
    logger.info("PSMs    : %s", PSM_FILE)

    # 2. Load dataset
    full_dataset = SpectralDataset(PSM_FILE, MGF_DIR, bin_size=0.1)

    if len(full_dataset) == 0:
        logger.error(
            "No labelled spectra found. Has 04_extract_psms.py finished yet? "
            "MGF dir: %s  PSM file: %s",
            MGF_DIR,
            PSM_FILE,
        )
        return

    # Audit summary before training
    audit = full_dataset.audit_summary()
    logger.info("Dataset audit: %s", json.dumps(audit))

    # 3. Peptide-grouped train/val/test split (audit fix)
    train_idx, val_idx, test_idx = _peptide_grouped_split(full_dataset.psms)

    train_dataset = Subset(full_dataset, train_idx)
    val_dataset = Subset(full_dataset, val_idx)
    test_dataset = Subset(full_dataset, test_idx)  # noqa: F841  (kept for reference)

    # Persist test set for locked evaluation
    test_psms = full_dataset.psms.iloc[test_idx]
    test_psms_path = os.path.join(CHECKPOINT_DIR, "test_set_psms.tsv")
    test_psms.to_csv(test_psms_path, sep="\t", index=False)
    logger.info("Test-set PSMs saved to %s (%d rows)", test_psms_path, len(test_psms))

    # Audit: report peptide overlap between train and test
    train_peps = set(full_dataset.psms.iloc[train_idx]["peptide"])
    test_peps = set(full_dataset.psms.iloc[test_idx]["peptide"])
    overlap = train_peps & test_peps
    if overlap:
        logger.warning(
            "LEAKAGE WARNING: %d peptides appear in both train and test splits.",
            len(overlap),
        )
    else:
        logger.info("Peptide overlap between train and test: 0 (clean split).")

    train_loader = DataLoader(
        train_dataset, batch_size=batch_size, shuffle=True, num_workers=num_workers
    )
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)

    # 4. Model, loss, optimiser, scheduler
    model = NeoepitopeSeq2Seq().to(DEVICE)
    criterion = nn.CrossEntropyLoss(ignore_index=0)
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)

    WARMUP = warmup_epochs

    def lr_lambda(epoch):
        if epoch < WARMUP:
            return float(epoch + 1) / WARMUP
        return 0.5 ** ((epoch - WARMUP) // 20)

    scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda)

    logger.info(
        "Training: %d train | %d val | %d test PSMs",
        len(train_dataset),
        len(val_dataset),
        len(test_dataset),
    )

    best_val_loss = float("inf")
    patience = 10
    trigger_times = 0

    for epoch in range(1, epochs + 1):
        model.train()
        total_loss = 0.0

        for batch_idx, (spectra, sequences) in enumerate(train_loader):
            spectra, sequences = spectra.to(DEVICE), sequences.to(DEVICE)

            optimizer.zero_grad()
            outputs = model(spectra)
            outputs = outputs.view(-1, VOCAB_SIZE)
            sequences = sequences.view(-1)
            loss = criterion(outputs, sequences)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()
            total_loss += loss.item()

            if batch_idx % 100 == 0:
                logger.info(
                    "Epoch [%d/%d] | Batch %d/%d | Loss: %.4f",
                    epoch,
                    epochs,
                    batch_idx,
                    len(train_loader),
                    loss.item(),
                )

        avg_train_loss = total_loss / len(train_loader)
        scheduler.step()

        # Validation
        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for spectra, sequences in val_loader:
                spectra, sequences = spectra.to(DEVICE), sequences.to(DEVICE)
                outputs = model(spectra).view(-1, VOCAB_SIZE)
                sequences = sequences.view(-1)
                val_loss += criterion(outputs, sequences).item()

        avg_val_loss = val_loss / len(val_loader)

        # Sample prediction for visual sanity check
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
            epoch,
            avg_train_loss,
            avg_val_loss,
            scheduler.get_last_lr()[0],
        )

        # Checkpointing and early stopping
        if avg_val_loss < best_val_loss:
            best_val_loss = avg_val_loss
            trigger_times = 0
            best_path = os.path.join(
                CHECKPOINT_DIR, "neoepitope_production_best.pth"
            )
            torch.save(model.state_dict(), best_path)
            _save_checkpoint_metadata(
                best_path, config, PSM_FILE, MGF_DIR, epoch, avg_val_loss, audit
            )
            logger.info("Saved BEST model → %s", best_path)
        else:
            trigger_times += 1
            if trigger_times >= patience:
                logger.info("Early stopping at epoch %d", epoch)
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

    # Final save
    final_path = os.path.join(CHECKPOINT_DIR, "neoepitope_production_final.pth")
    torch.save(model.state_dict(), final_path)
    _save_checkpoint_metadata(
        final_path, config, PSM_FILE, MGF_DIR, epoch, avg_val_loss, audit
    )
    logger.info("FINAL model saved → %s", final_path)


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
        default=os.path.join(PROJECT_ROOT, "results", "checkpoints_curated31_v2"),
    )
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--epochs", type=int, default=40)
    parser.add_argument("--lr", type=float, default=0.0005)
    parser.add_argument("--num-workers", type=int, default=0)
    args = parser.parse_args()

    train_production(
        mgf_dir=args.mgf_dir,
        psm_file=args.psm_file,
        checkpoint_dir=args.checkpoint_dir,
        batch_size=args.batch_size,
        epochs=args.epochs,
        lr=args.lr,
        num_workers=args.num_workers,
    )
