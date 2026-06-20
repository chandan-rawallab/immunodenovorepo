"""SpectralDataset: pairs MGF spectra with PSM peptide labels.

Audit fixes applied (2026-06-20):
  - Missing spectra are tracked and reported (count + fraction) instead of
    silently returning zero-padded tensors without any record.
  - Broad bare ``except`` replaced with structured logging that includes the
    exception type and scan ID so root-cause debugging is possible.
  - Unified spectrum lookup helper replaces ad-hoc multi-key probing loop.
  - Log/sqrt intensity transformation option added (default: log1p).
  - Top-N peak filtering applied before vectorisation (default: top 200 peaks).
  - Unknown amino acids (not in AA_TO_INT) now map to explicit UNK token (index 23)
    instead of collapsing silently into PAD (index 0).
  - Peptide truncation is logged with a warning when it occurs.
"""

from __future__ import annotations

import logging
import os

import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset

logger = logging.getLogger(__name__)

# Try to import pyteomics
try:
    from pyteomics import mgf
except ImportError:
    logger.warning("pyteomics not available. Falling back to local mgf_utils.")
    mgf = None

# Local imports from the same package
try:
    from .mgf_utils import IndexedMgfFallback
except ImportError:
    from mgf_utils import IndexedMgfFallback

# ---------------------------------------------------------------------------
# Vocabulary
# ---------------------------------------------------------------------------
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_INT: dict[str, int] = {aa: i + 1 for i, aa in enumerate(AMINO_ACIDS)}
AA_TO_INT["<PAD>"] = 0
AA_TO_INT["<START>"] = 21
AA_TO_INT["<END>"] = 22
AA_TO_INT["<UNK>"] = 23   # Explicit unknown token — audit fix
VOCAB_SIZE = len(AA_TO_INT)

# Keys tried in order when looking up a spectrum by scan ID
_SCAN_KEY_FORMATS = ("{scan}", "SCANS={scan}", "TITLE=scan={scan}")


class SpectralDataset(Dataset):
    """Dataset that pairs spectral vectors from MGF files with PSM peptide labels.

    Parameters
    ----------
    psm_file:
        Tab-separated file with columns: run_id, spectrum_id, peptide.
    mgf_dir:
        Directory containing per-run ``.mgf`` files.
    bin_size:
        m/z bin width (Da).
    max_mz:
        Maximum m/z value to consider.
    max_seq_len:
        Padded sequence length (tokens).
    top_n_peaks:
        Retain only the N most intense peaks before binning.
        ``None`` disables peak filtering.
    intensity_transform:
        ``"log1p"`` (default), ``"sqrt"``, or ``"none"``.
    """

    def __init__(
        self,
        psm_file: str,
        mgf_dir: str,
        bin_size: float = 0.1,
        max_mz: float = 2000.0,
        max_seq_len: int = 30,
        top_n_peaks: int | None = 200,
        intensity_transform: str = "log1p",
    ):
        super().__init__()
        self.bin_size = bin_size
        self.max_mz = max_mz
        self.vector_size = int(max_mz / bin_size)
        self.max_seq_len = max_seq_len
        self.top_n_peaks = top_n_peaks
        self.intensity_transform = intensity_transform

        # --- Load PSMs ---
        logger.info("Loading PSMs from %s", psm_file)
        self.psms = pd.read_csv(psm_file, sep="\t")

        self.mgf_dir = mgf_dir
        self.readers: dict[str, IndexedMgfFallback] = {}

        # --- Filter PSMs for which the MGF file exists ---
        unique_runs = self.psms["run_id"].unique()
        logger.info(
            "Dataset has %d PSMs across %d runs.", len(self.psms), len(unique_runs)
        )

        valid_indices = [
            idx
            for idx, row in self.psms.iterrows()
            if os.path.exists(os.path.join(mgf_dir, f"{row['run_id']}.mgf"))
        ]

        skipped = len(self.psms) - len(valid_indices)
        if skipped > 0:
            logger.warning(
                "%d PSMs skipped because the corresponding MGF file was not found.", skipped
            )
            self.psms = self.psms.iloc[valid_indices].reset_index(drop=True)

        # Audit counters
        self._missing_spectrum_count = 0
        self._truncation_count = 0

        logger.info(
            "SpectralDataset initialised with %d valid PSM pairs.", len(self.psms)
        )

    # ------------------------------------------------------------------
    # Properties for external audit reporting
    # ------------------------------------------------------------------
    @property
    def missing_spectrum_count(self) -> int:
        return self._missing_spectrum_count

    @property
    def truncation_count(self) -> int:
        return self._truncation_count

    def __len__(self) -> int:
        return len(self.psms)

    # ------------------------------------------------------------------
    # Spectrum preprocessing
    # ------------------------------------------------------------------
    def _apply_top_n_filter(
        self, mz_array: list, int_array: list
    ) -> tuple[list, list]:
        """Keep only the top-N most intense peaks."""
        if self.top_n_peaks is None or len(mz_array) <= self.top_n_peaks:
            return mz_array, int_array
        pairs = sorted(zip(int_array, mz_array), reverse=True)[: self.top_n_peaks]
        intensities, mzs = zip(*pairs)
        return list(mzs), list(intensities)

    def _transform_intensity(self, values: np.ndarray) -> np.ndarray:
        if self.intensity_transform == "log1p":
            return np.log1p(values)
        if self.intensity_transform == "sqrt":
            return np.sqrt(np.maximum(values, 0.0))
        return values  # "none"

    def _bin_spectrum(self, mz_array: list, int_array: list) -> torch.Tensor:
        mz_array, int_array = self._apply_top_n_filter(mz_array, int_array)

        vector = np.zeros(self.vector_size, dtype=np.float32)
        for mz, intensity in zip(mz_array, int_array):
            if mz < self.max_mz:
                bin_idx = int(mz / self.bin_size)
                if bin_idx < self.vector_size:
                    vector[bin_idx] += intensity

        vector = self._transform_intensity(vector)

        max_val = vector.max()
        if max_val > 0:
            vector /= max_val

        return torch.tensor(vector).unsqueeze(0)  # (1, vector_size)

    # ------------------------------------------------------------------
    # Sequence encoding
    # ------------------------------------------------------------------
    def _encode_sequence(self, sequence: str) -> torch.Tensor:
        raw_len = len(sequence) + 2  # +START +END
        if raw_len > self.max_seq_len:
            self._truncation_count += 1
            logger.warning(
                "Peptide '%s' (len=%d) will be truncated to fit max_seq_len=%d.",
                sequence,
                len(sequence),
                self.max_seq_len - 2,
            )

        unk_idx = AA_TO_INT["<UNK>"]
        tokens = [AA_TO_INT.get(c, unk_idx) for c in sequence]

        if any(t == unk_idx for t in tokens):
            unknown_chars = [c for c in sequence if c not in AA_TO_INT]
            logger.warning(
                "Peptide '%s' contains unknown amino acids %s; "
                "they will be encoded as <UNK> (index %d).",
                sequence,
                unknown_chars,
                unk_idx,
            )

        tokens = [AA_TO_INT["<START>"]] + tokens + [AA_TO_INT["<END>"]]
        pad_len = self.max_seq_len - len(tokens)
        if pad_len > 0:
            tokens.extend([AA_TO_INT["<PAD>"]] * pad_len)
        return torch.tensor(tokens[: self.max_seq_len], dtype=torch.long)

    # ------------------------------------------------------------------
    # Spectrum lookup
    # ------------------------------------------------------------------
    def _get_reader(self, mgf_path: str) -> IndexedMgfFallback:
        if mgf_path not in self.readers:
            self.readers[mgf_path] = IndexedMgfFallback(mgf_path)
        return self.readers[mgf_path]

    def _lookup_spectrum(self, reader: IndexedMgfFallback, scan_id: str) -> dict | None:
        """Try multiple key formats; return None if spectrum is absent."""
        for fmt in _SCAN_KEY_FORMATS:
            key = fmt.format(scan=scan_id)
            try:
                return reader.get_by_id(key)
            except KeyError:
                continue
            except Exception as exc:  # noqa: BLE001
                logger.debug(
                    "Unexpected error probing key '%s': %s: %s",
                    key,
                    type(exc).__name__,
                    exc,
                )
        return None

    # ------------------------------------------------------------------
    # __getitem__
    # ------------------------------------------------------------------
    def __getitem__(self, idx: int) -> tuple[torch.Tensor, torch.Tensor]:
        row = self.psms.iloc[idx]
        run_id = row["run_id"]
        scan_id = str(row["spectrum_id"])
        peptide = str(row["peptide"])

        mgf_path = os.path.join(self.mgf_dir, f"{run_id}.mgf")

        try:
            reader = self._get_reader(mgf_path)
            spectrum = self._lookup_spectrum(reader, scan_id)
        except Exception as exc:  # noqa: BLE001
            logger.error(
                "Failed to open/index MGF '%s': %s: %s",
                mgf_path,
                type(exc).__name__,
                exc,
            )
            spectrum = None

        if spectrum is None:
            self._missing_spectrum_count += 1
            logger.debug(
                "Spectrum not found: run=%s scan=%s (cumulative missing: %d)",
                run_id,
                scan_id,
                self._missing_spectrum_count,
            )
            return torch.zeros((1, self.vector_size)), self._encode_sequence("")

        mz_array = spectrum.get("m/z array", [])
        int_array = spectrum.get("intensity array", [])
        x = self._bin_spectrum(mz_array, int_array)
        y = self._encode_sequence(peptide)
        return x, y

    def audit_summary(self) -> dict:
        """Return a summary dict of data-quality counters for logging/reporting."""
        total = len(self.psms)
        return {
            "total_psms": total,
            "missing_spectrum_count": self._missing_spectrum_count,
            "missing_spectrum_fraction": self._missing_spectrum_count / total if total else 0.0,
            "truncation_count": self._truncation_count,
            "intensity_transform": self.intensity_transform,
            "top_n_peaks": self.top_n_peaks,
        }
