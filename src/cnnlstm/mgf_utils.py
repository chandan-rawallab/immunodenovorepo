"""Lightweight MGF parsing utilities.

The local Python environment may not have pyteomics installed, so the Objective 3
bridge keeps a tiny reader for SCANS/TITLE metadata and peak pairs.

Audit fixes applied (2026-06-20):
  - Duplicate spectrum detection and reporting instead of silent overwrite.
  - Structured KeyError with helpful message for missing spectrum IDs.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class MgfSpectrum:
    spectrum_id: str
    mz_array: list[float]
    intensity_array: list[float]
    params: dict[str, str]


def iter_mgf(path: Path):
    params: dict[str, str] = {}
    mz_array: list[float] = []
    intensity_array: list[float] = []
    in_block = False
    ordinal = 0

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            upper = line.upper()

            if upper == "BEGIN IONS":
                params = {}
                mz_array = []
                intensity_array = []
                in_block = True
                ordinal += 1
                continue

            if upper == "END IONS" and in_block:
                spectrum_id = (
                    params.get("SCANS")
                    or _scan_from_title(params.get("TITLE", ""))
                    or str(ordinal)
                )
                yield MgfSpectrum(
                    spectrum_id=spectrum_id,
                    mz_array=mz_array,
                    intensity_array=intensity_array,
                    params=params,
                )
                in_block = False
                continue

            if not in_block:
                continue

            if "=" in line:
                key, value = line.split("=", 1)
                params[key.strip().upper()] = value.strip()
                continue

            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                mz_array.append(float(parts[0]))
                intensity_array.append(float(parts[1]))
            except ValueError:
                continue


def _scan_from_title(title: str) -> str:
    for token in title.split():
        if token.startswith("scan="):
            return token.split("=", 1)[1]
    return ""


class IndexedMgfFallback:
    """Small get_by_id-compatible fallback for SpectralDataset.

    Audit fix: duplicate spectrum IDs are now detected and logged instead of
    silently overwriting the earlier entry. The total duplicate count is reported
    when the index is first built.
    """

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self._index: dict[str, MgfSpectrum] | None = None
        self.duplicate_count: int = 0

    def _ensure_index(self) -> None:
        if self._index is not None:
            return
        index: dict[str, MgfSpectrum] = {}
        duplicates: list[str] = []

        for spectrum in iter_mgf(self.path):
            sid = spectrum.spectrum_id
            # Register under bare ID and "SCANS=<id>" alias
            for key in (sid, f"SCANS={sid}"):
                if key in index:
                    duplicates.append(key)
                else:
                    index[key] = spectrum

        self._index = index
        self.duplicate_count = len(duplicates)

        if duplicates:
            logger.warning(
                "MGF file '%s': %d duplicate spectrum IDs detected "
                "(first 10: %s). Later entries were discarded.",
                self.path.name,
                len(duplicates),
                duplicates[:10],
            )
        else:
            logger.debug(
                "MGF index built for '%s': %d unique spectra.",
                self.path.name,
                len(index) // 2,  # each spectrum stored under 2 keys
            )

    def get_by_id(self, spectrum_id: str):
        self._ensure_index()
        assert self._index is not None

        spectrum = self._index.get(spectrum_id)
        if spectrum is None:
            available_sample = list(self._index.keys())[:5]
            raise KeyError(
                f"Spectrum ID '{spectrum_id}' not found in '{self.path.name}'. "
                f"Index size: {len(self._index) // 2} spectra. "
                f"Sample IDs present: {available_sample}"
            )

        return {
            "m/z array": spectrum.mz_array,
            "intensity array": spectrum.intensity_array,
            "params": {key.lower(): value for key, value in spectrum.params.items()},
        }


def mgf_read_fallback(path: str | Path):
    """Fallback generator yielding dicts compatible with pyteomics.mgf.read."""
    for spectrum in iter_mgf(Path(path)):
        yield {
            "m/z array": spectrum.mz_array,
            "intensity array": spectrum.intensity_array,
            "params": {key.lower(): value for key, value in spectrum.params.items()},
        }
