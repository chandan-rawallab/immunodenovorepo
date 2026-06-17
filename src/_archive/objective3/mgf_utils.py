"""Lightweight MGF parsing utilities.

The local Python environment may not have pyteomics installed, so the Objective 3
bridge keeps a tiny reader for SCANS/TITLE metadata and peak pairs.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


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
                spectrum_id = params.get("SCANS") or _scan_from_title(params.get("TITLE", "")) or str(ordinal)
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
    """Small get_by_id-compatible fallback for ProductionMgfDataset."""

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self._index: dict[str, MgfSpectrum] | None = None

    def _ensure_index(self) -> None:
        if self._index is not None:
            return
        index: dict[str, MgfSpectrum] = {}
        for spectrum in iter_mgf(self.path):
            index[spectrum.spectrum_id] = spectrum
            index[f"SCANS={spectrum.spectrum_id}"] = spectrum
        self._index = index

    def get_by_id(self, spectrum_id: str):
        self._ensure_index()
        assert self._index is not None
        spectrum = self._index[spectrum_id]
        return {
            "m/z array": spectrum.mz_array,
            "intensity array": spectrum.intensity_array,
            "params": {key.lower(): value for key, value in spectrum.params.items()},
        }

