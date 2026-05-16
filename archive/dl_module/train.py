import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
from dataset import MgfDataset
from model import NeoepitopeSeq2Seq
import os

def train():
    # 1. Configuration
    MGF_PATH = "/home/amity/hla_data_mgf/20160513_TIL1_R2.mgf"
    MSMS_PATH = "/home/amity/hla_data_raw/combined/txt/msms.txt"
    BATCH_SIZE = 16
    EPOCHS = 50
    LEARNING_RATE = 0.001
    DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    
    # 2. Check for data existence
    if not os.path.exists(MGF_PATH):
        print(f"ERROR: No MGF file found at {MGF_PATH}. Waiting for MaxQuant or conversion to finish.")
        # For demonstration/testing, we might use a mock run or skip loading
        return

    # 3. Data Loading
    print("Loading dataset...")
    dataset = MgfDataset(MGF_PATH, msms_path=MSMS_PATH)
    dataloader = DataLoader(dataset, batch_size=BATCH_SIZE, shuffle=True)

    # 4. Model Setup
    model = NeoepitopeSeq2Seq().to(DEVICE)
    optimizer = optim.AdamW(model.parameters(), lr=LEARNING_RATE)
    criterion = nn.CrossEntropyLoss(ignore_index=0) # Ignore <PAD> tokens

    # 5. Training Loop
    print(f"Starting training on {DEVICE}...")
    for epoch in range(EPOCHS):
        model.train()
        total_loss = 0
        for batch_idx, (spectra, targets) in enumerate(dataloader):
            spectra = spectra.to(DEVICE)
            targets = targets.to(DEVICE)
            
            optimizer.zero_grad()
            outputs = model(spectra) # Output shape: (Batch, SeqLocs, VocabSize)
            
            # CrossEntropyLoss expects (Batch, C, d1, d2...) or flat (Batch * SeqLocs, VocabSize)
            # targets shape: (Batch, Seq_Len)
            outputs = outputs.view(-1, outputs.size(-1))
            targets = targets.view(-1)
            
            loss = criterion(outputs, targets)
            loss.backward()
            optimizer.step()
            
            total_loss += loss.item()
            
        avg_loss = total_loss / len(dataloader)
        print(f"Epoch [{epoch+1}/{EPOCHS}], Average Loss: {avg_loss:.4f}")
        
        # Save checkpoint
        torch.save(model.state_dict(), f"neoepitope_model_epoch_{epoch+1}.pth")

if __name__ == "__main__":
    train()
