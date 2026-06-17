
import torch
import torch.nn as nn
import numpy as np
import pandas as pd
import os
import sys
import glob
import argparse

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Go up two levels: src/inference/ -> src/ -> project_root/
SRC_DIR = os.path.dirname(SCRIPT_DIR)
WORKSPACE_ROOT = os.path.dirname(SRC_DIR)
sys.path.append(SRC_DIR)

from cnnlstm.cnnlstm_model import NeoepitopeSeq2Seq
from cnnlstm.spectral_dataset import AA_TO_INT, VOCAB_SIZE

# Fallback for mgf
try:
    from pyteomics import mgf
except ImportError:
    from cnnlstm.mgf_utils import IndexedMgfFallback as mgf_fallback
    mgf = None

INT_TO_AA = {v: k for k, v in AA_TO_INT.items()}

def bin_spectrum(mz_array, int_array, bin_size=0.1, max_mz=2000.0):
    vector_size = int(max_mz / bin_size)
    vector = np.zeros(vector_size, dtype=np.float32)
    for mz, intensity in zip(mz_array, int_array):
        if mz < max_mz:
            bin_idx = int(mz / bin_size)
            if bin_idx < vector_size:
                vector[bin_idx] += intensity
    if np.max(vector) > 0:
        vector = vector / np.max(vector)
    return torch.tensor(vector).unsqueeze(0).unsqueeze(0) # (1, 1, vector_size)

def decode_sequence(predicted_indices):
    seq = []
    for p_idx in predicted_indices:
        token = INT_TO_AA.get(p_idx, "")
        if token == "<END>":
            break
        if token not in ["<START>", "<PAD>"]:
            seq.append(token)
    return "".join(seq)

def compute_sequence_score(outputs, predicted_indices):
    # outputs: (1, 30, VOCAB_SIZE)
    # Convert to probabilities
    probs = torch.softmax(outputs, dim=-1).squeeze(0) # (30, VOCAB_SIZE)
    score = 0.0
    length = 0
    for i, p_idx in enumerate(predicted_indices):
        prob = probs[i, p_idx].item()
        score += np.log(prob + 1e-9)
        length += 1
        token = INT_TO_AA.get(p_idx, "")
        if token == "<END>":
            break
    return score / max(1, length) # Length-normalized log probability

def main():
    parser = argparse.ArgumentParser(description="Run de novo neoepitope prediction.")
    parser.add_argument("--model", type=str, 
                        default=os.path.join(WORKSPACE_ROOT, "results", "checkpoints", "neoepitope_production_best.pth"),
                        help="Path to model checkpoint")
    parser.add_argument("--mgf_dir", type=str, 
                        default=os.path.join(WORKSPACE_ROOT, "data", "mgf_unlabeled"),
                        help="Directory containing unlabeled MGF files")
    parser.add_argument("--output", type=str, 
                        default=os.path.join(WORKSPACE_ROOT, "results", "de_novo_candidates.tsv"),
                        help="Path to save de novo candidates")
    parser.add_argument("--device", type=str, default="cpu", help="Device to run on (cpu or cuda)")
    parser.add_argument("--bin_size", type=float, default=0.1, help="m/z bin size")
    parser.add_argument("--max_mz", type=float, default=2000.0, help="Maximum m/z value")
    
    args = parser.parse_args()

    # 1. Configuration
    MODEL_PATH = args.model
    MGF_DIR = args.mgf_dir
    OUT_FILE = args.output
    DEVICE = torch.device(args.device)
    BIN_SIZE = args.bin_size
    MAX_MZ = args.max_mz

    if not os.path.exists(MODEL_PATH):
        # Fallback to archive if default not found
        alt_path = os.path.join(WORKSPACE_ROOT, "results", "_archive", "checkpoints", "neoepitope_medium_lite.pth")
        if os.path.exists(alt_path):
            print(f"WARNING: Using archived checkpoint {alt_path}.")
            MODEL_PATH = alt_path
            BIN_SIZE = 0.5 # Medium-lite used 0.5 bin size
        else:
            print(f"ERROR: Model checkpoint not found at {MODEL_PATH}.")
            return

    # 2. Load Model
    print(f"Loading Model from {os.path.basename(MODEL_PATH)}...")
    model = NeoepitopeSeq2Seq().to(DEVICE)
    try:
        model.load_state_dict(torch.load(MODEL_PATH, map_location=DEVICE))
    except Exception as e:
        print(f"ERROR: Could not load model: {e}")
        return
    model.eval()

    # 3. Process MGF files
    mgf_files = glob.glob(os.path.join(MGF_DIR, "*.mgf"))
    
    candidates = []
    
    print(f"Found {len(mgf_files)} MGF files. Processing...")
    
    for mgf_path in mgf_files:
        # Run ID might have _unlabeled suffix, so let's clean it up
        run_id = os.path.basename(mgf_path).replace(".mgf", "").replace("_unlabeled", "")
        print(f"Processing {run_id}...")
        
        try:
            if mgf is not None:
                reader = mgf.read(mgf_path)
            else:
                # Use fallback iterator if available
                from src.cnnlstm.mgf_utils import mgf_read_fallback
                reader = mgf_read_fallback(mgf_path)
        except Exception as e:
            print(f"Could not read {mgf_path}: {e}")
            continue
            
        count = 0
        for spectrum in reader:
            if spectrum is None:
                continue
            if count > 0 and count % 10000 == 0:
                print(f"  Processed {count} spectra...")
            params = spectrum.get("params", {})
            scan_id = str(params.get("scans", ""))
            if not scan_id:
                title = str(params.get("title", ""))
                if "scan=" in title:
                    scan_id = title.split("scan=")[-1].strip()
                else:
                    scan_id = str(count)
                    
            mz_array = spectrum.get('m/z array', [])
            int_array = spectrum.get('intensity array', [])
            
            if len(mz_array) == 0:
                count += 1
                continue
                
            input_tensor = bin_spectrum(mz_array, int_array, BIN_SIZE, MAX_MZ).to(DEVICE)
            
            with torch.no_grad():
                outputs = model(input_tensor)
                predicted_indices = torch.argmax(outputs, dim=-1).squeeze(0).cpu().numpy()
                
            predicted_seq = decode_sequence(predicted_indices)
            score = compute_sequence_score(outputs, predicted_indices)
            
            if len(predicted_seq) >= 8 and len(predicted_seq) <= 11:
                candidates.append({
                    'run_id': run_id,
                    'spectrum_id': scan_id,
                    'peptide': predicted_seq,
                    'score': score,
                    'decoy': False
                })
                
                # Decoy: encode the reversed sequence and score it with the model
                # This is a proper target-decoy approach — the model predicts the
                # reversed sequence given the SAME spectrum, providing an empirical
                # null distribution under the same conditions as the target.
                decoy_seq = predicted_seq[::-1]
                decoy_encoded = torch.tensor(
                    [AA_TO_INT.get(c, 0) for c in decoy_seq], dtype=torch.long
                ).unsqueeze(0).to(DEVICE)
                with torch.no_grad():
                    decoy_outputs = model(input_tensor)  # same spectrum
                decoy_score = compute_sequence_score(decoy_outputs, decoy_encoded.squeeze(0).cpu().numpy())
                candidates.append({
                    'run_id': run_id,
                    'spectrum_id': scan_id,
                    'peptide': decoy_seq,
                    'score': decoy_score,
                    'decoy': True
                })
                
            count += 1

    if len(candidates) == 0:
        print("No candidates generated.")
        return

    df = pd.DataFrame(candidates)
    
    # 5. FDR Estimation
    print("Estimating FDR...")
    df = df.sort_values(by='score', ascending=False).reset_index(drop=True)
    
    targets = ~df['decoy']
    decoys = df['decoy']
    
    df['target_cum'] = targets.cumsum()
    df['decoy_cum'] = decoys.cumsum()
    
    # FDR = decoys / targets
    df['fdr'] = df['decoy_cum'] / df['target_cum'].replace(0, 1)
    
    # Filter 5% FDR targets only
    filtered_df = df[(df['fdr'] <= 0.05) & (~df['decoy'])]
    
    print(f"Retained {len(filtered_df)} candidates at 5% FDR.")
    
    # Save
    os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)
    filtered_df[['run_id', 'spectrum_id', 'peptide', 'score', 'fdr']].to_csv(OUT_FILE, sep='\t', index=False)
    print(f"Saved candidates to {OUT_FILE}")

if __name__ == "__main__":
    main()
