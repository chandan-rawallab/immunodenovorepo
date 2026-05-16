import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split, Subset
import os
import sys
import numpy as np
import glob
import time

# =========================================================================
# DYNAMIC PATH RESOLUTION
# =========================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.append(PROJECT_ROOT)

from dl_module.model import NeoepitopeSeq2Seq
from dl_module.production_dataset import ProductionMgfDataset, VOCAB_SIZE, AMINO_ACIDS

def train_medium_lite():
    # 1. Configuration
    MGF_DIR = "/home/amity/hla_data_mgf"
    RESULTS_DIR = os.path.join(PROJECT_ROOT, "production_results")
    CHECKPOINT_DIR = os.path.join(PROJECT_ROOT, "results", "checkpoints")
    LOG_DIR = os.path.join(PROJECT_ROOT, "logs")
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)
    
    LOG_FILE = os.path.join(LOG_DIR, "train_v2_medium_lite.log")
    
    def log(msg):
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
        formatted = f"[{timestamp}] {msg}"
        print(formatted)
        with open(LOG_FILE, "a") as f:
            f.write(formatted + "\n")

    DEVICE = torch.device("cpu") 
    log(f"Initializing MEDIUM-LITE Production Training on {DEVICE}...")
    log(f"Strategy: 1000 Samples, 0.5 Bin Resolution (Balanced Resource Mode)")

    # 2. Load Dataset with 0.5 bin size (2x more detail than Ultra-Lite)
    # We'll use the first 5 MGF files to increase diversity
    mgf_files = sorted(glob.glob(os.path.join(MGF_DIR, "*.mgf")))[:5]
    
    TEMP_MGF_DIR = os.path.join(PROJECT_ROOT, "temp", "medium_lite_mgf")
    os.makedirs(TEMP_MGF_DIR, exist_ok=True)
    for f in mgf_files:
        dst = os.path.join(TEMP_MGF_DIR, os.path.basename(f))
        if not os.path.exists(dst):
            os.symlink(f, dst)

    log(f"Indexing MGF files: {', '.join([os.path.basename(f) for f in mgf_files])}")
    full_dataset = ProductionMgfDataset(TEMP_MGF_DIR, RESULTS_DIR, bin_size=0.5)
    total_len = len(full_dataset)
    
    if total_len == 0:
        log("ERROR: No labeled spectra found in selected files.")
        return

    # 3. Create Subset (1,000 samples)
    subset_size = min(2000, total_len)
    indices = np.random.choice(total_len, subset_size, replace=False)
    subset_dataset = Subset(full_dataset, indices)
    
    log(f"Subsampled {subset_size} samples from {total_len} total available pairs.")

    # Split into Train and Validation (80/20)
    train_size = int(0.8 * subset_size)
    val_size = subset_size - train_size
    train_ds, val_ds = random_split(subset_dataset, [train_size, val_size])

    train_loader = DataLoader(train_ds, batch_size=16, shuffle=True, num_workers=0)
    val_loader = DataLoader(val_ds, batch_size=16, shuffle=False, num_workers=0)

    # 4. Initialize Model
    # Input vector_size will be 2000 / 0.5 = 4000
    model = NeoepitopeSeq2Seq().to(DEVICE)
    criterion = nn.CrossEntropyLoss(ignore_index=0) 
    optimizer = torch.optim.AdamW(model.parameters(), lr=0.001)

    # 5. Training Loop (5 Epochs)
    epochs = 5
    log(f"Starting Training: {train_size} train, {val_size} val.")
    
    for epoch in range(1, epochs + 1):
        model.train()
        total_loss = 0
        for batch_idx, (spectra, sequences) in enumerate(train_loader):
            spectra, sequences = spectra.to(DEVICE), sequences.to(DEVICE)
            
            optimizer.zero_grad()
            outputs = model(spectra).view(-1, VOCAB_SIZE)
            sequences = sequences.view(-1)
            
            loss = criterion(outputs, sequences)
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
            
            if (batch_idx + 1) % 10 == 0:
                log(f"Epoch {epoch} | Batch {batch_idx+1}/{len(train_loader)} | Loss: {loss.item():.4f}")

        avg_train_loss = total_loss / len(train_loader)
        log(f"Epoch {epoch}/{epochs} Completed | Avg Train Loss: {avg_train_loss:.4f}")

    # 6. Save Model
    save_path = os.path.join(CHECKPOINT_DIR, "neoepitope_medium_lite.pth")
    torch.save(model.state_dict(), save_path)
    log(f"--- MEDIUM-LITE model saved to {save_path} ---")

if __name__ == "__main__":
    train_medium_lite()
