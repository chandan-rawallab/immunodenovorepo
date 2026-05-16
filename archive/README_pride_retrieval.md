# PRIDE Data Retrieval Module

## Overview
This module securely retrieves raw mass spectrometry files required for the Neoepitope pipeline from the PRIDE archive (e.g., project `PXD005231`).

It has been custom-engineered to bypass Python's memory limitations by offloading the actual binary transfers to native OS mechanisms (`wget` via subprocess) while keeping the flexibility of a single Python entry point.

## Execution
Run the orchestrator script directly in your terminal. You must specify the PRIDE ID and the absolute path to your desired output sandbox.

```bash
cd /home/amity/Documents/experiments
python3 dataretrieval.py -i PXD005231 -o /home/amity/hla_data_raw
```

- **Logging**: The script automatically generates a `retrieval.log` file inside your PXD folder to track every step and transfer.

## ⚠️ Critical Architecture Warning
**Always designate an output directory OUTSIDE the main IDE workspace.**
Do not download massive analytical binaries directly into your code workspace (e.g., `Documents/experiments`). If you do, your code editor's background indexing engine will attempt to parse dozens of gigabytes of proprietary binary data, resulting in catastrophic Out-of-Memory (OOM) crashes and system freezing. 

Keep pipeline code in the workspace; keep raw big-data in isolated folders.

## Resumption & Stability
The execution process is idempotent and crash-proof:
- If your internet drops or the system restarts, simply run the exact same python command again. 
- The module will query the remote servers, safely skip any fully downloaded files, and seamlessly resume any partially downloaded binaries right where it left off.

## Phase 2: Mass Spec Conversion (mzML & MGF)
After downloading the `.raw` files, they must be converted into open formats for PEAKS X and Deep Learning (CNN/LSTM) analysis.

### Execution
Use the `conversion_pipeline.py` script to generate both formats for all files in a single pass. The script uses `ThermoRawFileParser` natively on Linux.

```bash
cd /home/amity/Documents/experiments
python3 conversion_pipeline.py \
    -i /home/amity/hla_data_raw/PXD005231 \
    --mzml /home/amity/hla_data_mzml \
    --mgf /home/amity/hla_data_mgf
```

- **Outputs**: Loss-less gzipped `.mzML` for database searching and flat-text `.mgf` for deep learning ingestion.
- **Resource Management**: The script automatically applies `renice` and `ionice` priorities so that the heavy CPU/IO parsing operations happen in the background without freezing your desktop.

## Finding Dataset IDs for Objective 3
**Objective 3** requires testing/training data from three specific studies. To test the pipeline on these other datasets, you must locate their `PXD` identifiers:
1. **Bassani-Sternberg et al., 2016:** This is the default dataset you just downloaded! Its PRIDE ID is `PXD005231`.
2. **Carreno et al., 2015 & Laumont et al., 2018:** To find the IDs for these remaining studies:
   - Go to the **Data Availability** section at the very end of those respective published papers.
   - Alternatively, search the [PRIDE Archive Search Portal](https://www.ebi.ac.uk/pride/) for the authors' names ("Laumont" or "Carreno") or the paper titles.
   - Once you locate the `PXDxxxxxx` code, simply run this script again passing that new ID!
