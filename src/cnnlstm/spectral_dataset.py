import torch
from torch.utils.data import Dataset
import numpy as np
import pandas as pd
import os
import glob

# Try to import pyteomics
try:
    from pyteomics import mgf
except ImportError:
    print("WARNING: pyteomics library not found. Falling back to local mgf_utils.")
    mgf = None

# Local imports from the same package
try:
    from .mgf_utils import IndexedMgfFallback
except ImportError:
    # Fallback for direct script execution
    from mgf_utils import IndexedMgfFallback

# Standard Amino Acid dictionary for Neoepitopes
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_INT = {aa: i + 1 for i, aa in enumerate(AMINO_ACIDS)}
AA_TO_INT['<PAD>'] = 0
AA_TO_INT['<START>'] = 21
AA_TO_INT['<END>'] = 22
VOCAB_SIZE = len(AA_TO_INT)

class SpectralDataset(Dataset):
    """
    A dataset that pairs spectral vectors from MGF files with peptide labels from PSM extractions.
    """
    def __init__(self, psm_file, mgf_dir, bin_size=0.1, max_mz=2000.0, max_seq_len=30):
        super().__init__()
        self.bin_size = bin_size
        self.max_mz = max_mz
        self.vector_size = int(max_mz / bin_size)
        self.max_seq_len = max_seq_len
        
        # Load extracted PSMs
        print(f"Loading PSMs from {psm_file}...")
        self.psms = pd.read_csv(psm_file, sep='\t')
        
        # Cache for MGF readers
        self.mgf_dir = mgf_dir
        self.readers = {} 
        
        # Validate that runs exist in MGF dir
        unique_runs = self.psms['run_id'].unique()
        print(f"Dataset has {len(self.psms)} PSMs across {len(unique_runs)} runs.")
        
        # Filter PSMs for which MGF file exists
        valid_indices = []
        for idx, row in self.psms.iterrows():
            run_id = row['run_id']
            mgf_path = os.path.join(mgf_dir, f"{run_id}.mgf")
            if os.path.exists(mgf_path):
                valid_indices.append(idx)
        
        if len(valid_indices) < len(self.psms):
            print(f"WARNING: {len(self.psms) - len(valid_indices)} PSMs skipped because MGF file was missing.")
            self.psms = self.psms.iloc[valid_indices].reset_index(drop=True)
            
        print(f"Final SpectralDataset initialized with {len(self.psms)} valid pairs.")
        
    def __len__(self):
        return len(self.psms)
        
    def _bin_spectrum(self, mz_array, int_array):
        vector = np.zeros(self.vector_size, dtype=np.float32)
        for mz, intensity in zip(mz_array, int_array):
            if mz < self.max_mz:
                bin_idx = int(mz / self.bin_size)
                if bin_idx < self.vector_size:
                    vector[bin_idx] += intensity
        
        if np.max(vector) > 0:
            vector = vector / np.max(vector)
        return torch.tensor(vector).unsqueeze(0) # (1, vector_size)
        
    def _encode_sequence(self, sequence):
        tokens = [AA_TO_INT.get(c, 0) for c in list(sequence)]
        tokens = [AA_TO_INT['<START>']] + tokens + [AA_TO_INT['<END>']]
        pad_len = self.max_seq_len - len(tokens)
        if pad_len > 0:
            tokens.extend([AA_TO_INT['<PAD>']] * pad_len)
        return torch.tensor(tokens[:self.max_seq_len], dtype=torch.long)

    def __getitem__(self, idx):
        row = self.psms.iloc[idx]
        run_id = row['run_id']
        scan_id = str(row['spectrum_id'])
        peptide = str(row['peptide'])
        
        mgf_path = os.path.join(self.mgf_dir, f"{run_id}.mgf")
        
        # Use a cached reader per file
        if mgf_path not in self.readers:
            if IndexedMgfFallback is not None:
                self.readers[mgf_path] = IndexedMgfFallback(mgf_path)
            elif mgf is not None:
                self.readers[mgf_path] = mgf.IndexedMGF(mgf_path)
            else:
                raise RuntimeError("No MGF reader available.")
        
        reader = self.readers[mgf_path]
        try:
            # Try different ID formats
            spectrum = None
            for key in [f"SCANS={scan_id}", scan_id, f"TITLE=scan={scan_id}"]:
                try:
                    spectrum = reader.get_by_id(key)
                    if spectrum: break
                except:
                    continue
            
            if spectrum is None:
                # Return zero vector if spectrum not found to avoid crashing training
                return torch.zeros((1, self.vector_size)), self._encode_sequence("")
                
            mz_array = spectrum.get('m/z array', [])
            int_array = spectrum.get('intensity array', [])
            
            x = self._bin_spectrum(mz_array, int_array)
            y = self._encode_sequence(peptide)
            return x, y
            
        except Exception as e:
            # Fail gracefully
            return torch.zeros((1, self.vector_size)), self._encode_sequence("")
