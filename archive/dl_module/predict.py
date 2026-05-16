import torch
import torch.nn as nn
import numpy as np
from model import NeoepitopeSeq2Seq
from dataset import MgfDataset, AA_TO_INT, AMINO_ACIDS
import os
import sys

# Ensure weight path is correct Relative to current file
MODEL_PATH = os.path.join(os.path.dirname(__file__), "..", "neoepitope_model_epoch_50.pth")

# Create a reverse dictionary for decoding
INT_TO_AA = {v: k for k, v in AA_TO_INT.items()}

def predict():
    # 1. Configuration
    MGF_PATH = "/home/amity/hla_data_mgf/20160513_TIL1_R2.mgf"
    DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    if not os.path.exists(MODEL_PATH):
        print(f"ERROR: Model file {MODEL_PATH} not found.")
        return

    # 2. Load Model
    print("Loading model architecture and weights...")
    model = NeoepitopeSeq2Seq().to(DEVICE)
    model.load_state_dict(torch.load(MODEL_PATH, map_location=DEVICE))
    model.eval()

    # 3. Load Dataset (To get a real spectrum)
    print("Loading a sample spectrum from MGF...")
    dataset = MgfDataset(MGF_PATH) 
    
    # Let's pick a random spectrum that likely has signal
    sample_idx = 100 # Can be any valid index
    spectrum_tensor, actual_labels = dataset[sample_idx]
    
    # 4. Perform Prediction
    print(f"Predicting sequence for Spectrum #{sample_idx}...")
    with torch.no_grad():
        input_tensor = spectrum_tensor.unsqueeze(0).to(DEVICE) # Add batch dim
        outputs = model(input_tensor) # (1, 30, VocabSize)
        
        # Get the most likely indices
        # outputs is (Batch, SeqLen, VocabSize)
        predicted_indices = torch.argmax(outputs, dim=-1).squeeze(0).cpu().numpy()

    # 5. Decode results
    predicted_seq = []
    for idx in predicted_indices:
        token = INT_TO_AA.get(idx, "")
        if token == "<END>":
            break
        if token not in ["<START>", "<PAD>"]:
            predicted_seq.append(token)
    
    predicted_str = "".join(predicted_seq)
    
    print("\n" + "="*40)
    print(f"RESULTS FOR SPECTRUM INDEX: {sample_idx}")
    print(f"DE NOVO SEQUENCE: {predicted_str}")
    print("="*40)
    print("Note: This is a de novo prediction. For neoantigen prioritization,")
    print("compare this sequence against the human reference proteome to find mutations.")

if __name__ == "__main__":
    predict()
