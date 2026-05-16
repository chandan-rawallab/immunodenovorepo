import torch
from torch.utils.data import Dataset
import numpy as np
import pandas as pd
import os
import glob

try:
    from pyteomics import mgf
except ImportError:
    print("WARNING: pyteomics library not found. Please run: pip install pyteomics numpy torch pandas")
    mgf = None

try:
    from objective3.mgf_utils import IndexedMgfFallback
except ImportError:
    try:
        import sys
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
        from objective3.mgf_utils import IndexedMgfFallback
    except ImportError:
        IndexedMgfFallback = None

# Standard Amino Acid dictionary for Neoepitopes
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_INT = {aa: i + 1 for i, aa in enumerate(AMINO_ACIDS)}
AA_TO_INT['<PAD>'] = 0
AA_TO_INT['<START>'] = 21
AA_TO_INT['<END>'] = 22
VOCAB_SIZE = len(AA_TO_INT)

class ProductionMgfDataset(Dataset):
    """
    A Production-scale PyTorch Dataset that manages multiple MGF files and MaxQuant results.
    """
    def __init__(self, mgf_dir, results_dir, bin_size=0.1, max_mz=2000.0, max_seq_len=30):
        super().__init__()
        self.bin_size = bin_size
        self.max_mz = max_mz
        self.vector_size = int(max_mz / bin_size)
        self.max_seq_len = max_seq_len
        
        self.samples = [] # List of (mgf_path, scan_id, sequence)
        self.readers = {} # Cache for mgf readers
        
        # 1. Discover all MGF files
        mgf_files = glob.glob(os.path.join(mgf_dir, "*.mgf"))
        print(f"Found {len(mgf_files)} MGF files in {mgf_dir}")
        
        # 2. Match with corresponding production results
        for m_path in mgf_files:
            basename = os.path.basename(m_path).replace(".mgf", "")
            # Our scale_search.py saves results as msms_{basename}.raw.txt
            # Wait, the scale_search.py script used: f"msms_{raw_filename}.txt" where raw_filename included .raw
            # So it looks like: msms_20160513_TIL1_R2.raw.txt
            msms_pattern = os.path.join(results_dir, f"msms_{basename}*.txt")
            msms_files = glob.glob(msms_pattern)
            
            if not msms_files:
                print(f"Skipping {basename}: No identification results found in {results_dir}")
                continue
                
            msms_path = msms_files[0]
            try:
                # Use only required columns to save memory
                df = pd.read_csv(msms_path, sep='\t', low_memory=False, usecols=['Scan number', 'Sequence'])
                for _, row in df.iterrows():
                    # Ensure scan number is correctly formatted for lookup
                    scan = str(int(row['Scan number']))
                    seq = str(row['Sequence'])
                    # Store reference to this specific pair
                    self.samples.append((m_path, scan, seq))
                del df # Free memory immediately
            except Exception as e:
                print(f"Error reading {msms_path}: {e}")

        print(f"Production dataset initialized with {len(self.samples)} labeled spectra pairs.")
        
    def __len__(self):
        return len(self.samples)
        
    def _bin_spectrum(self, mz_array, int_array):
        vector = np.zeros(self.vector_size, dtype=np.float32)
        for mz, intensity in zip(mz_array, int_array):
            if mz < self.max_mz:
                bin_idx = int(mz / self.bin_size)
                # Ensure bin_idx is within bounds
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
        mgf_path, scan_id, sequence = self.samples[idx]
        
        # Use a cached reader per file
        if mgf_path not in self.readers:
            if IndexedMgfFallback is not None:
                self.readers[mgf_path] = IndexedMgfFallback(mgf_path)
            elif mgf is not None:
                self.readers[mgf_path] = mgf.IndexedMGF(mgf_path)
            else:
                raise RuntimeError("No MGF reader available.")
        
        reader = self.readers[mgf_path]
        # Get spectrum by scan number (IndexedMgfFallback uses 'scans' field as key)
        try:
            spectrum = reader.get_by_id(f"SCANS={scan_id}")
        except (KeyError, AssertionError):
            try:
                spectrum = reader.get_by_id(scan_id)
            except (KeyError, AssertionError):
                return torch.zeros((1, self.vector_size)), self._encode_sequence("")

        mz_array = spectrum.get('m/z array', [])
        int_array = spectrum.get('intensity array', [])
        
        x = self._bin_spectrum(mz_array, int_array)
        y = self._encode_sequence(sequence)
        
        return x, y
