import torch
from torch.utils.data import Dataset
import numpy as np
import pandas as pd
import os

try:
    from pyteomics import mgf
except ImportError:
    print("WARNING: pyteomics library not found. Please run: pip install pyteomics numpy torch pandas")

# Standard Amino Acid dictionary for Neoepitopes
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_INT = {aa: i + 1 for i, aa in enumerate(AMINO_ACIDS)}
AA_TO_INT['<PAD>'] = 0
AA_TO_INT['<START>'] = 21
AA_TO_INT['<END>'] = 22
VOCAB_SIZE = len(AA_TO_INT)

class MgfDataset(Dataset):
    """
    A PyTorch Dataset that lazily parses large MGF spectra files and maps them to MaxQuant identifications.
    """
    def __init__(self, mgf_path, msms_path=None, bin_size=0.1, max_mz=2000.0, max_seq_len=30):
        super().__init__()
        self.mgf_path = mgf_path
        self.bin_size = bin_size
        self.max_mz = max_mz
        self.vector_size = int(max_mz / bin_size)
        self.max_seq_len = max_seq_len
        
        # Load Labels from MaxQuant msms.txt
        self.labels = {}
        if msms_path and os.path.exists(msms_path):
            print(f"Loading identifications from {msms_path}...")
            # msms.txt is tab-separated. We need Scan number (39) and Sequence (3)
            df = pd.read_csv(msms_path, sep='\t', low_memory=False)
            # Filter for this specific raw file if necessary, but here we assume single run
            for _, row in df.iterrows():
                scan = str(int(row['Scan number']))
                seq = str(row['Sequence'])
                self.labels[scan] = seq
            print(f"Mapped {len(self.labels)} peptide sequences.")

        print(f"Indexing MGF spectra lazily from {mgf_path}...")
        self.reader = mgf.read(mgf_path)
        self.spectra = []
        
        # To maintain efficiency, we only keep spectra that have a corresponding label
        # or we keep all for inference (but we need training data now)
        for spec in self.reader:
            params = spec.get('params', {})
            scan = params.get('scans')
            if msms_path:
                if scan in self.labels:
                    self.spectra.append(spec)
            else:
                self.spectra.append(spec)
                
        print(f"Dataset ready with {len(self.spectra)} labeled training pairs.")
        
    def __len__(self):
        return len(self.spectra)
        
    def _bin_spectrum(self, mz_array, int_array):
        vector = np.zeros(self.vector_size, dtype=np.float32)
        for mz, intensity in zip(mz_array, int_array):
            if mz < self.max_mz:
                bin_idx = int(mz / self.bin_size)
                vector[bin_idx] += intensity
        
        if np.max(vector) > 0:
            vector = vector / np.max(vector)
        return torch.tensor(vector).unsqueeze(0) # Add Channel Dim: (1, vector_size)
        
    def _encode_sequence(self, sequence):
        tokens = [AA_TO_INT.get(c, 0) for c in list(sequence)]
        tokens = [AA_TO_INT['<START>']] + tokens + [AA_TO_INT['<END>']]
        pad_len = self.max_seq_len - len(tokens)
        if pad_len > 0:
            tokens.extend([AA_TO_INT['<PAD>']] * pad_len)
        return torch.tensor(tokens[:self.max_seq_len], dtype=torch.long)

    def __getitem__(self, idx):
        spectrum = self.spectra[idx]
        mz_array = spectrum.get('m/z array', [])
        int_array = spectrum.get('intensity array', [])
        
        # NN Input
        x = self._bin_spectrum(mz_array, int_array)
        
        # NN Target
        scan = spectrum.get('params', {}).get('scans')
        target_seq = self.labels.get(scan, "PAD")
        y = self._encode_sequence(target_seq)
        
        return x, y
