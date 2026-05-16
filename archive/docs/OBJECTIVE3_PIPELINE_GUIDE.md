# Objective 3: Hybrid Neoepitope Discovery Pipeline Guide

This document provides step-by-step instructions for the complete manual setup and execution of the Immunopeptidomics pipeline, including the Deep Learning (DL) de novo sequencing modules.

---

## 1. Prerequisites & System Requirements
- **OS**: Ubuntu 22.04+ (or compatible Linux)
- **Hardware**: 32GB+ RAM, NVIDIA GPU (Recommended for DL steps)
- **Software**: .NET 8.0 SDK/Runtime, Conda/Miniconda.

### Step 0: Manual Tool Extraction (OneDrive Downloads)
The large binaries in the `bin/` folder are zipped for transport. You MUST extract them before running any scripts:
```bash
cd /home/amity/Documents/experiments/bin/
unzip MaxQuant_v2.8.0.0.zip
unzip ThermoRawFileParser_v2.0.0.zip
unzip dotnet-sdk-8.0.zip
```

---

## 2. Environment Setup
Create the dedicated environment for the Deep Learning modules:
```bash
conda create -n obj3-dl python=3.10
conda activate obj3-dl
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
pip install pandas numpy pyteomics tqdm scikit-learn
```

---

## 3. Data Acquisition (PRIDE/MassIVE)
Download RAW data directly into the sandbox:
```bash
# PRIDE Example
python3 src/utils/dataretrieval.py -i PXD005231
# MassIVE Example
python3 src/utils/dataretrieval.py -i MSV000080620
```

---

## 4. Database Search (MaxQuant)
Run the sequential search to establish the baseline of "known" peptides:
```bash
python3 src/scale_search.py
```
*Outputs: `production_results/msms.txt`.*

---

## 5. Deep Learning (DL) Pipeline
This section handles the discovery of non-canonical (mutated) peptides.

### Step A: Spectral Preprocessing
Convert RAW files to MGF format for neural network consumption:
```bash
python3 src/utils/conversion_pipeline.py
```

### Step B: De Novo Prediction (Inference)
The `NeoepitopeSeq2Seq` model predicts amino acid sequences directly from MS spectra.
```bash
# Use the CLI entrypoint
python3 -m objective3.cli predict-denovo \
  --mgf-dir /home/amity/hla_data_mgf \
  --model-path src/models/neoepitope_weights.pth
```

### Step C: Database Subtraction
We filter the DL predictions against the MaxQuant results to ensure we only keep novel candidates.
```bash
python3 -m objective3.cli subtract-hits \
  --mq-results production_results/ \
  --dl-results results/objective3/denovo_out.tsv
```

### Step D: Multi-Evidence Ranking
Final neoantigens are ranked using a composite score:
1.  **Model Confidence**: Probability from the Seq2Seq model.
2.  **Binding Affinity**: (External tool integration like MHCflurry).
3.  **Spectral Support**: Intensity and fragment match scores.

```bash
python3 -m objective3.cli rank-neoantigens \
  --input results/objective3/subtracted_candidates.tsv \
  --output results/ranked_neoantigens.tsv
```

---

## 6. Troubleshooting
- **MKL Errors**: If you see `iJIT_NotifyEvent` errors, run:
  `export MKL_SERVICE_FORCE_INTEL=1`
- **Memory Issues**: MaxQuant is memory-intensive. The `scale_search.py` script is set to sequential mode by default to prevent system crashes.
- **Missing FASTAs**: Ensure the `uniprot_human_reviewed.fasta` is in the `src/configs/` directory.

### Process Management (Tracking & Killing)
When running long pipeline scripts like `scale_search.py` in the background (e.g., using `nohup ... &`), you cannot see the output directly in the terminal.

**1. How to Track Progress:**
To watch a background process live, view the end of its log file:
```bash
# Watch the log update in real-time (Press Ctrl+C to stop watching)
tail -f scale_search_output.log
```

**2. How to Check if it is Running:**
To see if the Python orchestrator or the MaxQuant engine are active:
```bash
ps aux | grep -E "python3 src/scale_search.py|dotnet" | grep -v grep
```

**3. How to Force Kill (If Hung):**
If you accidentally start multiple instances (which causes deadlocks) or if a process is permanently frozen:
```bash
# Force kill the python orchestrator
pkill -9 -f "python3 src/scale_search.py"

# Force kill all background MaxQuant worker threads
pkill -9 -f "dotnet"
```
*(Note: Always check for and delete corrupted temporary folders like `data/raw/combined` if you forcefully kill a process.)*
