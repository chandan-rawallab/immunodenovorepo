import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split, Subset
import os
import sys
import numpy as np
import glob

# =========================================================================
# DYNAMIC PATH RESOLUTION
# =========================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.append(PROJECT_ROOT)

from dl_module.model import NeoepitopeSeq2Seq
from dl_module.production_dataset import ProductionMgfDataset, VOCAB_SIZE, AMINO_ACIDS

def train_ultra_lite():
    # 1. Configuration
    MGF_DIR = "/home/amity/hla_data_mgf"
    RESULTS_DIR = os.path.join(PROJECT_ROOT, "production_results")
    CHECKPOINT_DIR = os.path.join(PROJECT_ROOT, "results", "checkpoints")
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    
    DEVICE = torch.device("cpu") # Force CPU to avoid any GPU overhead/drivers issues
    print(f"Initializing ULTRA-LITE Production Training on {DEVICE}...")
    print(f"Strategy: High-level Subsampling + Resolution Reduction (100x Less Compute)")

    # 2. Load Dataset with 1.0 bin size (10x less RAM than 0.1)
    # We also manually limit the MGF discovery to the first 2 files to save index RAM
    mgf_files = sorted(glob.glob(os.path.join(MGF_DIR, "*.mgf")))[:2]
    
    # We'll patch ProductionMgfDataset to only use these files
    # Instead of modifying the class, we'll just pass a temp dir with symlinks or similar?
    # Simpler: Create a temp directory with just 2 MGF files.
    TEMP_MGF_DIR = os.path.join(PROJECT_ROOT, "temp", "lite_mgf")
    os.makedirs(TEMP_MGF_DIR, exist_ok=True)
    for f in mgf_files:
        dst = os.path.join(TEMP_MGF_DIR, os.path.basename(f))
        if not os.path.exists(dst):
            os.symlink(f, dst)

    full_dataset = ProductionMgfDataset(TEMP_MGF_DIR, RESULTS_DIR, bin_size=1.0)
    total_len = len(full_dataset)
    
    if total_len == 0:
        print("ERROR: No labeled spectra found in selected files.")
        return

    # 3. Create Subset (100 samples total)
    subset_size = min(100, total_len)
    indices = np.random.choice(total_len, subset_size, replace=False)
    subset_dataset = Subset(full_dataset, indices)
    
    print(f"Subsampled {subset_size} samples from {total_len} total.")

    # Split into Train and Validation (80/20)
    train_size = int(0.8 * subset_size)
    val_size = subset_size - train_size
    train_ds, val_ds = random_split(subset_dataset, [train_size, val_size])

    train_loader = DataLoader(train_ds, batch_size=8, shuffle=True, num_workers=0)
    val_loader = DataLoader(val_ds, batch_size=8, shuffle=False, num_workers=0)

    # 4. Initialize Model
    # Note: Model forward pass is adaptive, so it handles vector_size=2000 automatically
    model = NeoepitopeSeq2Seq().to(DEVICE)
    criterion = nn.CrossEntropyLoss(ignore_index=0) 
    optimizer = torch.optim.AdamW(model.parameters(), lr=0.001)

    # 5. Training Loop (3 Epochs)
    epochs = 3
    print(f"Starting Ultra-Lite Training: {train_size} train, {val_size} val.")
    
    for epoch in range(1, epochs + 1):
        model.train()
        total_loss = 0
        for batch_idx, (spectra, sequences) in enumerate(train_loader):
            # Adapt input if necessary (though NeoepitopeSeq2Seq is already generic in its CNN)
            # Actually, the CNN layer has a fixed in_channels=1.
            # The input shape is (Batch, 1, 2000).
            
            spectra, sequences = spectra.to(DEVICE), sequences.to(DEVICE)
            
            optimizer.zero_grad()
            outputs = model(spectra).view(-1, VOCAB_SIZE)
            sequences = sequences.view(-1)
            
            loss = criterion(outputs, sequences)
            loss.backward()
            optimizer.step()
            total_loss += loss.item()

        avg_train_loss = total_loss / len(train_loader)
        print(f"Epoch {epoch}/{epochs} | Train Loss: {avg_train_loss:.4f}")

    # 6. Save Model
    lite_path = os.path.join(CHECKPOINT_DIR, "neoepitope_ultra_lite.pth")
    torch.save(model.state_dict(), lite_path)
    print(f"--- ULTRA-LITE model saved to {lite_path} ---")

if __name__ == "__main__":
    train_ultra_lite()
