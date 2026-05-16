### [PROJECT STATUS] Production Ready (v2.1)
- **Directory Structure**: Fully reorganized into Standard Academic Layout (`src/`, `bin/`, `data/`, `results/`).
- **Portability**: All scripts (`scale_search.py`, `train_production.py`, `conversion_pipeline.py`) transitioned to Relative Pathing.
- **Acquisition**: Upgraded `dataretrieval.py` to Universal Mode (PRIDE + MassIVE support).
- **Environment**: .NET 8.0 and Native Linux ThermoRawFileParser (v2.0.0-dev) validated.

---

### [CHANGELOG] 2026-04-27
#### [UPGRADE] Universal Data Retrieval
- Added FTP scraping logic for UCSD MassIVE repository.
- Unified EBI/UCSD downloading into a single CLI tool (`dataretrieval.py`).
- Implemented `wget`--recursive--no-parent mirroring for MSV IDs.

#### [RESTRUCTURE] Academic Pipeline Architecture
- Migrated all code logic from root to `src/`.
- Moved project configurations to `src/configs/`.
- Isolated deep learning architecture into `src/dl_module/`.
- Standardized data paths to `data/raw/` and `data/mgf/`.

---

## Phase 1: Data Acquisition
The pipeline supports automated retrieval of immunopeptidomics datasets from two major repositories: PRIDE (via EBI API) and MassIVE (via FTP scraping).

*   **Command:** `python3 src/utils/dataretrieval.py -i <Accession_ID>`
*   **Storage Path:** `data/raw/<ID>/`

## Phase 2: Spectra Processing (RAW to MGF)
While MaxQuant can ingest standard `.raw` files natively, our neural networks require standardized `.mgf` format for custom manipulation.

*   **Engine**: ThermoRawFileParser v2.0.0-dev (Native Linux binary).
*   **Command:** `python3 src/utils/conversion_pipeline.py`
*   **Output Path:** `data/mgf/` (Recursive discovery of all downloaded RAW files).

## Phase 3: Sequential Production Search
To prevent memory overflow (OOM) on large datasets (even on 64-core servers), we implemented a sequential orchestrator that search files one-by-one and aggregates results.

*   **Command:** `python3 src/scale_search.py`
*   **Results:** `results/identifications/msms_<filename>.txt`

## Phase 4: Production Deep Learning Training
High-memory training utilizing the standard identification labels matched with the standardized spectra.

*   **Command:** `python3 src/train_production.py`
*   **Checkpoints:** `results/checkpoints/`
