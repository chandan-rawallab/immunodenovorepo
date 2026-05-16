import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split
import os
import sys

# =========================================================================
# DYNAMIC PATH RESOLUTION (FOR ACADEMIC STRUCTURE)
# =========================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
# Add PROJECT_ROOT to sys.path so we can import from dl_module
sys.path.append(PROJECT_ROOT)

from cnnlstm.cnnlstm_model import NeoepitopeSeq2Seq
from cnnlstm.spectral_dataset import SpectralDataset, VOCAB_SIZE, AMINO_ACIDS

def train_production(
    mgf_dir=None,
    psm_file=None,
    checkpoint_dir=None,
    batch_size=32,
    epochs=100,
    lr=0.001,
    num_workers=4,
):
    # 1. Configuration (Relative to Project Structure)
    MGF_DIR = mgf_dir or os.path.join(PROJECT_ROOT, "data", "mgf")
    PSM_FILE = psm_file or os.path.join(PROJECT_ROOT, "results", "immunopeptidome_psms.tsv")
    CHECKPOINT_DIR = checkpoint_dir or os.path.join(PROJECT_ROOT, "results", "checkpoints")
    
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    
    BATCH_SIZE = batch_size
    EPOCHS = epochs
    LR = lr
    DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    print(f"Initializing Production Training on {DEVICE}...")
    print(f"Reading Spectra from: {MGF_DIR}")
    print(f"Reading PSMs from: {PSM_FILE}")

    # 2. Load Dataset
    full_dataset = SpectralDataset(PSM_FILE, MGF_DIR, bin_size=0.1)
    
    if len(full_dataset) == 0:
        print("ERROR: No labeled spectra found. Has 04_extract_psms.py finished yet?")
        print(f"Check your .mgf files in {MGF_DIR} and PSM file {PSM_FILE}")
        return

    # Split into Train and Validation (80/20)
    train_size = int(0.8 * len(full_dataset))
    val_size = len(full_dataset) - train_size
    train_dataset, val_dataset = random_split(full_dataset, [train_size, val_size])

    train_loader = DataLoader(train_dataset, batch_size=BATCH_SIZE, shuffle=True, num_workers=num_workers)
    val_loader = DataLoader(val_dataset, batch_size=BATCH_SIZE, shuffle=False)

    # 3. Initialize Model, Loss, Optimizer, Scheduler
    model = NeoepitopeSeq2Seq().to(DEVICE)
    criterion = nn.CrossEntropyLoss(ignore_index=0) 
    optimizer = torch.optim.AdamW(model.parameters(), lr=LR)
    scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=20, gamma=0.5)

    # 4. Training Loop
    print(f"Starting Training: {train_size} train samples, {val_size} val samples.")
    best_val_loss = float('inf')
    patience = 10
    trigger_times = 0
    
    for epoch in range(1, EPOCHS + 1):
        model.train()
        total_loss = 0
        for batch_idx, (spectra, sequences) in enumerate(train_loader):
            spectra, sequences = spectra.to(DEVICE), sequences.to(DEVICE)
            
            optimizer.zero_grad()
            outputs = model(spectra) 
            
            # Reshape for CrossEntropy
            outputs = outputs.view(-1, VOCAB_SIZE)
            sequences = sequences.view(-1)
            
            loss = criterion(outputs, sequences)
            loss.backward()
            optimizer.step()
            
            total_loss += loss.item()

            if batch_idx % 100 == 0:
                print(f"Epoch [{epoch}/{EPOCHS}] | Batch {batch_idx}/{len(train_loader)} | Loss: {loss.item():.4f}")

        avg_train_loss = total_loss / len(train_loader)
        scheduler.step()
        
        # Validation Phase
        model.eval()
        val_loss = 0
        with torch.no_grad():
            for spectra, sequences in val_loader:
                spectra, sequences = spectra.to(DEVICE), sequences.to(DEVICE)
                outputs = model(spectra).view(-1, VOCAB_SIZE)
                sequences = sequences.view(-1)
                loss = criterion(outputs, sequences)
                val_loss += loss.item()
        
        avg_val_loss = val_loss / len(val_loader)
        
        # Visual Verification: Print one sample prediction
        with torch.no_grad():
            sample_spectra, sample_seq = next(iter(val_loader))
            sample_spectra = sample_spectra[0].unsqueeze(0).to(DEVICE)
            pred = model(sample_spectra).argmax(dim=-1).squeeze().tolist()
            # Convert tokens back to amino acids, handling special tokens
            pred_chars = []
            for t in pred:
                if t == 0: continue # PAD
                if t == 21: pred_chars.append("<S>") # START
                elif t == 22: pred_chars.append("<E>") # END
                elif 0 < t <= len(AMINO_ACIDS):
                    pred_chars.append(AMINO_ACIDS[t-1])
                else:
                    pred_chars.append("?")
            print(f"Sample Prediction (Epoch {epoch}): {''.join(pred_chars)}")

        print(f"==> Epoch {epoch} Complete | Train Loss: {avg_train_loss:.4f} | Val Loss: {avg_val_loss:.4f} | LR: {scheduler.get_last_lr()[0]}")

        # Checkpointing and Early Stopping
        if avg_val_loss < best_val_loss:
            best_val_loss = avg_val_loss
            trigger_times = 0
            best_path = os.path.join(CHECKPOINT_DIR, "neoepitope_production_best.pth")
            torch.save(model.state_dict(), best_path)
            print(f"--- Saved BEST model to {best_path} ---")
        else:
            trigger_times += 1
            if trigger_times >= patience:
                print(f"Early stopping at epoch {epoch}")
                break

        if epoch % 10 == 0:
            ckpt_path = os.path.join(CHECKPOINT_DIR, f"neoepitope_production_epoch_{epoch}.pth")
            torch.save(model.state_dict(), ckpt_path)
            print(f"--- Saved periodic checkpoint to {ckpt_path} ---")

    # Final Save
    final_path = os.path.join(CHECKPOINT_DIR, "neoepitope_production_final.pth")
    torch.save(model.state_dict(), final_path)
    print(f"--- FINAL model saved to {final_path} ---")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Train the Objective 3 CNN/LSTM model on production MGF + MaxQuant labels.")
    parser.add_argument("--mgf-dir", default="/home/amity/hla_data_mgf")
    parser.add_argument("--results-dir", default=os.path.join(PROJECT_ROOT, "production_results"))
    parser.add_argument("--checkpoint-dir", default=os.path.join(PROJECT_ROOT, "results", "checkpoints"))
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--lr", type=float, default=0.001)
    parser.add_argument("--num-workers", type=int, default=4)
    args = parser.parse_args()

    train_production(
        mgf_dir=args.mgf_dir,
        results_dir=args.results_dir,
        checkpoint_dir=args.checkpoint_dir,
        batch_size=args.batch_size,
        epochs=args.epochs,
        lr=args.lr,
        num_workers=args.num_workers,
    )
