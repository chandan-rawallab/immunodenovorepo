import torch
import torch.nn as nn
import numpy as np
import os
import sys
import glob

# =========================================================================
# DYNAMIC PATH RESOLUTION
# =========================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.append(PROJECT_ROOT)

from src.cnnlstm.cnnlstm_model import NeoepitopeSeq2Seq
from src.cnnlstm.spectral_dataset import SpectralDataset, AA_TO_INT, AMINO_ACIDS, VOCAB_SIZE

# Create a reverse dictionary for decoding
INT_TO_AA = {v: k for k, v in AA_TO_INT.items()}

def predict_production():
    # 1. Configuration
    CHECKPOINT_DIR = os.path.join(PROJECT_ROOT, "results", "checkpoints")
    # Using the new standardized checkpoint name from the plan
    MODEL_PATH = os.path.join(PROJECT_ROOT, "results", "model_checkpoint.pth")
    MGF_DIR = os.path.join(PROJECT_ROOT, "data", "mgf")
    PSM_FILE = os.path.join(PROJECT_ROOT, "results", "immunopeptidome_psms.tsv")
    DEVICE = torch.device("cpu")
    BIN_SIZE = 0.1 # Standardized

    if not os.path.exists(MODEL_PATH):
        # Check archive if not in root
        MODEL_PATH = os.path.join(PROJECT_ROOT, "results", "_archive", "checkpoints", "neoepitope_medium_lite.pth")
        if os.path.exists(MODEL_PATH):
            print(f"WARNING: Using archived checkpoint {MODEL_PATH}. This might have bin_size=0.5.")
            BIN_SIZE = 0.5
        else:
            print("ERROR: No checkpoints found.")
            return

    # 2. Load Model
    print(f"Loading Model from {os.path.basename(MODEL_PATH)}...")
    model = NeoepitopeSeq2Seq().to(DEVICE)
    try:
        model.load_state_dict(torch.load(MODEL_PATH, map_location=DEVICE))
    except Exception as e:
        print(f"ERROR: Could not load model. It may have been trained with a different bin_size. {e}")
        return
    model.eval()

    # 3. Load Dataset
    print(f"Initializing dataset (bin_size={BIN_SIZE})...")
    # Use SpectralDataset from cnnlstm package
    dataset = SpectralDataset(PSM_FILE, MGF_DIR, bin_size=BIN_SIZE)
    
    if len(dataset) == 0:
        print("ERROR: No spectra found.")
        return

    # 4. Perform Prediction on 10 random samples
    print(f"Predicting 10 samples...")
    print("-" * 60)
    
    indices = np.random.choice(len(dataset), 10, replace=False)
    
    for i, idx in enumerate(indices):
        spectrum_tensor, actual_labels = dataset[idx]
        
        with torch.no_grad():
            input_tensor = spectrum_tensor.unsqueeze(0).to(DEVICE)
            outputs = model(input_tensor)
            predicted_indices = torch.argmax(outputs, dim=-1).squeeze(0).cpu().numpy()

        # Decode Predicted
        predicted_seq = []
        for p_idx in predicted_indices:
            token = INT_TO_AA.get(p_idx, "")
            if token == "<END>":
                break
            if token not in ["<START>", "<PAD>"]:
                predicted_seq.append(token)
        predicted_str = "".join(predicted_seq)
        
        # Decode Actual
        actual_seq = []
        for a_idx in actual_labels.numpy():
            token = INT_TO_AA.get(a_idx, "")
            if token == "<END>":
                break
            if token not in ["<START>", "<PAD>"]:
                actual_seq.append(token)
        actual_str = "".join(actual_seq)

        print(f"Sample {i+1} [Index {idx}]:")
        print(f"  ACTUAL:    {actual_str}")
        print(f"  PREDICTED: {predicted_str}")
    
    print("-" * 60)

if __name__ == "__main__":
    predict_production()
