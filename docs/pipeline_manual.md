# Objective 3: Production-Scale Neoepitope Discovery Manual
**Project:** Immunogenicity Prediction & De Novo Sequencing Pipeline
**Goal:** High-throughput identification of mutated peptides for personalized vaccines.

---

## 📅 Chronological Execution Log (Date-wise)

### 2026-05-11 to 2026-05-12: Data Acquisition (PRIDE PXD005231)
*   **Activity:** Retrieval of ~100GB of patient-derived mass spectrometry raw data.
*   **Manual Steps:** 
    1.  Go to [PRIDE Archive](https://www.ebi.ac.uk/pride/archive/projects/PXD005231).
    2.  Use the "FTP Download" option or Aspera Client.
    3.  Select all `.raw` files from the "Files" tab.
    4.  Manually download and organize into folders by patient ID.
*   **Background Action:** The agent used `dataretrieval.py` with `wget` subprocesses to handle multi-threaded downloads and automatic resumption.
*   **Screenshot Placeholder:** [Image: PRIDE Project Page showing PXD005231 and the list of associated RAW files].

### 2026-05-13: High-Throughput MaxQuant Screening
*   **Activity:** Automated identification of peptides from 40+ patient-derived `.raw` files to create the "Ground Truth" dataset.
*   **Manual Steps:**
    1.  Copy `.raw` files into the `data/raw/` directory.
    2.  Open MaxQuant GUI.
    3.  Load files and configure `mqpar.xml` (Enzyme: Trypsin, Mods: Oxidation, Fixed: Carbamidomethyl).
    4.  Click "Start" and wait for "combined/" folder generation.
*   **Background Action:** The agent utilized `scale_search.py` to index and process multiple runs in parallel using the `MaxQuantCmd.exe`.
*   **Screenshot Placeholder:** [Image: MaxQuant "msms.txt" file showing Peptide Sequences and PEP scores].

### 2026-05-14: Spectral Signal Processing (RAW → MGF)
*   **Activity:** Extraction of peak lists (M/Z, Intensity) from binary raw data for Deep Learning input.
*   **Manual Steps:**
    1.  Open MSConvert (ProteoWizard).
    2.  Add `.raw` files to the input list.
    3.  Add "Peak Picking" filter (Vendor, Level 1-2).
    4.  Set output format to "Mascot generic (mgf)".
    5.  Start conversion.
*   **Background Action:** Automated conversion using `ThermoRawFileParser` via the `conversion_pipeline.py` script.
*   **Screenshot Placeholder:** [Image: A spectrum visualization in a viewer like SeeMS showing M/Z peaks].

### 2026-05-15: Resource-Optimized Model Training (Ultra-Lite)
*   **Activity:** Training the `NeoepitopeSeq2Seq` model under strict 16GB RAM constraints.
*   **Manual Steps:**
    1.  Open a Python IDE (e.g., PyCharm or VSCode).
    2.  Write a DataLoader that reads MGF files and maps them to MaxQuant sequences.
    3.  Define a Transformer or CNN-LSTM architecture in PyTorch.
    4.  Manually adjust batch size and resolution if the system starts swapping (thrashing).
*   **Background Action:** 
    - **Optimization:** Initialized **"Ultra-Lite" Training** (`train_lite.py`).
    - **Technique:** Reduced bin resolution from 0.1 to 1.0 (10x RAM saving) and used a 100-sample pilot subset to ensure the pipeline functions without OOM crashes.
*   **Screenshot Placeholder:** [Image: Console output showing "Epoch 1/5... Loss: 2.45" and memory usage monitor].

### 2026-05-15 (Current): Medium-Lite Scaling (1,000 Samples)
*   **Activity:** Scaling the de novo model to 1,000 samples and 0.5 bin resolution for improved sequence prediction accuracy.
*   **Manual Steps:**
    1.  Increase the number of MGF files indexed in the DataLoader (use at least 5 patient samples).
    2.  Set `bin_size=0.5` to capture double the spectral detail compared to Ultra-Lite.
    3.  Run the training loop for 5 epochs using a larger batch size (e.g., 16).
*   **Background Action:** The agent launched `train_v2.py` in the background using `nohup` to ensure continuity.
*   **Screenshot Placeholder:** [Image: Log output showing "Subsampled 1000 samples" and successful training start].

---

## 🛠️ Reimplementation Guide (Step-by-Step)

### Step 1: Environment Preparation
```bash
# Activate the specialized bioinformatics environment
conda activate objective3-casanovo
pip install torch numpy pandas pyteomics
```

### Step 2: Data Preparation (Manual equivalent)
1.  Ensure your MGF files are in `/home/amity/hla_data_mgf`.
2.  Ensure MaxQuant `msms.txt` results are in `/home/amity/Documents/experiments/production_results/`.

### Step 3: Training the Model (Medium-Lite Mode)
To replicate the current optimized training state:
1.  Navigate to the project root.
2.  Run the medium-lite training script:
    ```bash
    python3 src/train_v2.py
    ```
    *Note: This script uses 1,000 samples and 0.5 binning, providing a balance between detail and RAM usage on 16GB machines.*

### Step 4: Verification (Inference)
To test if the model can predict a sequence from a spectrum:
```bash
python3 src/predict_lite.py
```

---

## 📜 Standing Rules for Objective 3
1.  **Memory Guard:** Always check RAM (`free -h`) before training. If `< 2GB` free, use `Lite` scripts.
2.  **Traceability:** MGF files MUST keep the original RAW filename prefix (e.g., `Patient1_SampleA.mgf`).
3.  **No Placeholders:** All code must be functional; if a dataset is too large, use stratified sampling rather than empty placeholders.
4.  **Audit Logs:** All background processes must output to `/logs/` with a unique timestamp.

---

## 📌 Maintenance & Next Steps
*This manual is a living document. The current focus is validating the **Medium-Lite model** (1,000 samples) and implementing the **Proteome Comparison** phase to filter non-human (neoepitope) hits.*
