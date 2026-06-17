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

## 🔬 Manual MaxQuant Setup & Execution on Linux (Without `scale_search.py`)

> **When to use this:** You have freshly installed MaxQuant on a Linux machine and want to run a database search against your `.raw` files **without** our automated `scale_search.py` orchestrator. This covers everything from zero.

### Prerequisites

#### A. Install .NET 8 Runtime (Required — MaxQuant is a C# application)

MaxQuant v2.8+ runs on .NET 8. It will **not** work with older .NET versions or Mono.

```bash
# Download the official .NET installer
wget -q https://dot.net/v1/dotnet-install.sh -O dotnet-install.sh
chmod +x dotnet-install.sh

# Install .NET 8 runtime into a local directory (no sudo needed)
./dotnet-install.sh --channel 8.0 --install-dir ./bin/dotnet_runtime

# Verify it works
./bin/dotnet_runtime/dotnet --info
```

You should see output containing `Microsoft.NETCore.App 8.x.x`. If you see an error about `libicu` or `libssl`, install them:
```bash
# Ubuntu/Debian
sudo apt install -y libicu-dev libssl-dev

# Fedora/RHEL
sudo dnf install -y libicu openssl-devel
```

#### B. Install MaxQuant v2.8.0.0

MaxQuant does **not** have a native Linux installer. You must:

1.  Download from [maxquant.org](https://www.maxquant.org/download_asset/maxquant/latest) (requires free registration).
2.  Unzip it into a known directory:
    ```bash
    mkdir -p bin/MaxQuant_v2.8.0.0
    unzip MaxQuant_2.8.0.0.zip -d bin/MaxQuant_v2.8.0.0
    ```
3.  Verify the DLL exists:
    ```bash
    ls bin/MaxQuant_v2.8.0.0/bin/MaxQuantCmd.dll
    # Should print the path without errors
    ```

#### C. Download the Reference Proteome (FASTA)

MaxQuant needs a protein database to search against. We use the UniProt Human Reviewed (Swiss-Prot) proteome:

```bash
mkdir -p data/reference
wget -c "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(reviewed:true)+AND+(model_organism:9606)" \
     -O data/reference/uniprot_human_reviewed.fasta
```

---

### Creating `mqpar.xml` from Scratch on Linux

Since MaxQuant has **no GUI on Linux**, you must create the parameter file (`mqpar.xml`) by hand. There are two approaches:

#### Option 1: Copy and Edit the Production Template (Recommended)

We already have a working template at `configs/mqpar_production.xml`. Copy it and modify the key fields:

```bash
cp configs/mqpar_production.xml configs/my_search.xml
```

Then edit `configs/my_search.xml` with any text editor (`nano`, `vim`, `code`, etc.).

#### Option 2: Generate on Windows, Transfer to Linux

If you have temporary access to a Windows machine:
1.  Open MaxQuant GUI on Windows.
2.  Load one `.raw` file and configure your settings.
3.  Go to **File → Save Parameters** to export `mqpar.xml`.
4.  Transfer the file to Linux via `scp` or USB.
5.  **You MUST then fix all paths** (see below).

#### Option 3: Write From Scratch

Use the minimal template below as a starting point. This is the absolute minimum required XML structure:

```xml
<?xml version="1.0" encoding="utf-8"?>
<MaxQuantParams xmlns:xsd="http://www.w3.org/2001/XMLSchema"
                xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">

   <!-- ═══════════════════════════════════════════════════ -->
   <!-- SECTION 1: Reference Database (FASTA)              -->
   <!-- ═══════════════════════════════════════════════════ -->
   <fastaFiles>
      <FastaFileInfo>
         <!-- ⚠️  ABSOLUTE PATH REQUIRED — relative paths WILL fail -->
         <fastaFilePath>/home/amity/Documents/experiments/data/reference/uniprot_human_reviewed.fasta</fastaFilePath>
         <identifierParseRule>>[^|]*\|(.*?)\|</identifierParseRule>
         <descriptionParseRule>>(.*)</descriptionParseRule>
         <taxonomyParseRule>OX=(\d+)</taxonomyParseRule>
         <variationParseRule></variationParseRule>
         <modificationParseRule></modificationParseRule>
         <taxonomyId></taxonomyId>
      </FastaFileInfo>
   </fastaFiles>

   <!-- ═══════════════════════════════════════════════════ -->
   <!-- SECTION 2: Input Raw Files                         -->
   <!-- ═══════════════════════════════════════════════════ -->
   <!-- ⚠️  ABSOLUTE PATHS REQUIRED for EVERY raw file    -->
   <filePaths>
      <string>/home/amity/hla_data_raw/20160513_TIL1_R2.raw</string>
      <!-- Add more <string>...</string> entries for additional files -->
   </filePaths>

   <!-- Experiment names: one per raw file, in the same order -->
   <experiments>
      <string>TIL1_R2</string>
   </experiments>

   <!-- Fractions: one per raw file (use 32767 for "no fractionation") -->
   <fractions>
      <short>32767</short>
   </fractions>

   <!-- Parameter group index: one per raw file (usually all 0) -->
   <paramGroupIndices>
      <int>0</int>
   </paramGroupIndices>

   <!-- ═══════════════════════════════════════════════════ -->
   <!-- SECTION 3: Global Search Settings                  -->
   <!-- ═══════════════════════════════════════════════════ -->
   <numThreads>2</numThreads>        <!-- Keep low (2-4) on 16GB RAM -->
   <decoyMode>revert</decoyMode>
   <includeContaminants>True</includeContaminants>
   <minPeptideLength>7</minPeptideLength>
   <peptideFdr>0.01</peptideFdr>
   <proteinFdr>0.01</proteinFdr>
   <matchBetweenRuns>False</matchBetweenRuns>
   <maxQuantVersion>2.8.0.0</maxQuantVersion>

   <!-- ═══════════════════════════════════════════════════ -->
   <!-- SECTION 4: Parameter Group (Enzyme & Modifications)-->
   <!-- ═══════════════════════════════════════════════════ -->
   <parameterGroups>
      <parameterGroup>
         <enzymeMode>2</enzymeMode>  <!-- 2 = semi-specific (for immunopeptidomics) -->
         <maxMissedCleavages>2</maxMissedCleavages>

         <enzymes>
            <string>Trypsin/P</string>
         </enzymes>

         <!-- No fixed mods for immunopeptidomics (no alkylation step) -->
         <fixedModifications>
         </fixedModifications>

         <!-- Variable modifications -->
         <variableModifications>
            <string>Oxidation (M)</string>
            <string>Acetyl (Protein N-term)</string>
         </variableModifications>

         <!-- Mass tolerance -->
         <firstSearchTol>20</firstSearchTol>     <!-- ppm -->
         <mainSearchTol>4.5</mainSearchTol>       <!-- ppm -->
         <searchTolInPpm>True</searchTolInPpm>
         <maxCharge>7</maxCharge>
         <multiplicity>1</multiplicity>           <!-- 1 = label-free -->
         <lcmsRunType>Standard</lcmsRunType>
      </parameterGroup>
   </parameterGroups>

</MaxQuantParams>
```

> **Note:** The minimal template above omits many optional sections (DIA, cross-linking, TIMS, etc.) that MaxQuant fills with defaults. For a complete production-ready file, use `configs/mqpar_production.xml` as your base.

---

### Critical XML Tags You Must Edit

Every time you create or reuse a `mqpar.xml`, check these tags:

| Tag | What It Does | ⚠️ Common Mistake |
|:---|:---|:---|
| `<fastaFilePath>` | Path to UniProt FASTA | Must be **absolute** (`/home/...`). Relative paths silently fail. |
| `<filePaths><string>` | Path(s) to `.raw` files | Must be **absolute**. One `<string>` per file. |
| `<experiments><string>` | Experiment label per raw file | Must have **exactly** as many entries as `<filePaths>`. |
| `<fractions><short>` | Fraction number per raw file | Use `32767` for single-shot (no fractionation). |
| `<paramGroupIndices><int>` | Links each file to a parameter group | Usually `0` for all files. |
| `<numThreads>` | CPU threads for the search | Set to `2` on 16GB RAM. Use `4–8` on 32GB+. |
| `<enzymeMode>` | `0`=specific, `2`=semi-specific | Use `2` for immunopeptidomics (MHC-I peptides are not tryptic). |
| `<fixedModifications>` | Always-present chemical mods | Leave **empty** for immunopeptidomics (no iodoacetamide step). |

---

### Running MaxQuant (Single File, No `scale_search.py`)

Once your `mqpar.xml` is ready:

```bash
# Step 1: Set the .NET runtime path
export PATH=$PWD/bin/dotnet_runtime:$PATH

# Step 2: Run MaxQuant via the command-line DLL
./bin/dotnet_runtime/dotnet ./bin/MaxQuant_v2.8.0.0/bin/MaxQuantCmd.dll \
    configs/my_search.xml
```

**What happens:**
- MaxQuant creates a `combined/` folder **next to your raw files** (e.g., `/home/amity/hla_data_raw/combined/`).
- Inside `combined/txt/` you'll find `msms.txt` — this is the PSM output our pipeline needs.
- Runtime: ~30–90 minutes per `.raw` file depending on size and threads.

**To monitor progress:**
```bash
# MaxQuant writes progress to stdout. If running in background:
nohup ./bin/dotnet_runtime/dotnet ./bin/MaxQuant_v2.8.0.0/bin/MaxQuantCmd.dll \
    configs/my_search.xml > logs/maxquant_run.log 2>&1 &

# Watch the log
tail -f logs/maxquant_run.log
```

---

### Running Multiple Files Sequentially (Without `scale_search.py`)

To avoid OOM on 16GB RAM, process files **one at a time** with a simple bash loop:

```bash
export PATH=$PWD/bin/dotnet_runtime:$PATH

for RAW_FILE in /home/amity/hla_data_raw/*.raw; do
    BASENAME=$(basename "$RAW_FILE" .raw)
    echo "━━━ Processing: $BASENAME ━━━"

    # 1. Create a per-file mqpar.xml from the template
    sed "s|<string>/home/amity/hla_data_raw/.*\.raw</string>|<string>${RAW_FILE}</string>|" \
        configs/mqpar_production.xml > "/tmp/mqpar_${BASENAME}.xml"

    # 2. Update the experiment name
    sed -i "s|<string>.*</string>\(.*experiments\)|<string>${BASENAME}</string>\1|" \
        "/tmp/mqpar_${BASENAME}.xml"

    # 3. Run MaxQuant
    ./bin/dotnet_runtime/dotnet ./bin/MaxQuant_v2.8.0.0/bin/MaxQuantCmd.dll \
        "/tmp/mqpar_${BASENAME}.xml" 2>&1 | tee "logs/maxquant_${BASENAME}.log"

    # 4. Copy the msms.txt result into our production_results directory
    COMBINED_DIR="$(dirname "$RAW_FILE")/combined/txt"
    if [ -f "$COMBINED_DIR/msms.txt" ]; then
        cp "$COMBINED_DIR/msms.txt" "production_results/msms_${BASENAME}.txt"
        echo "✓ Saved: production_results/msms_${BASENAME}.txt"
        # Clean up MaxQuant temp files to free disk space
        rm -rf "$(dirname "$RAW_FILE")/combined"
    else
        echo "✗ ERROR: No msms.txt found for $BASENAME"
    fi
done
```

---

### After MaxQuant: Feeding Results into the Pipeline

Once you have `msms_*.txt` files in `production_results/` (or `data/psms/`), the pipeline picks them up automatically:

```bash
# Extract PSMs from MaxQuant results
PYTHONPATH=src python3 src/data_prep/04_extract_psms.py \
    --input-dir data/psms \
    --output-file results/immunopeptidome_psms.tsv \
    --manifest configs/sample_manifest.tsv
```

This produces `results/immunopeptidome_psms.tsv` — the "ground truth" normal peptide list that the rest of the pipeline (Steps 05–09) depends on.

---

### Troubleshooting

| Problem | Cause | Fix |
|:---|:---|:---|
| `NullReferenceException` on startup | Relative paths in `mqpar.xml` | Change ALL paths to absolute (`/home/amity/...`) |
| `Cannot find FASTA file` | Wrong `<fastaFilePath>` | Verify with `ls -la <the_path>` |
| OOM / Killed by kernel | Too many threads or files | Set `<numThreads>2</numThreads>`, process one file at a time |
| `msms.txt` is empty | Wrong enzyme mode for immunopeptidomics | Set `<enzymeMode>2</enzymeMode>` (semi-specific) |
| `dotnet: command not found` | .NET not in PATH | Run `export PATH=$PWD/bin/dotnet_runtime:$PATH` first |
| Search takes >4 hours | Normal for large files with semi-specific digest | Use `<numThreads>4</numThreads>` if RAM allows |

---

## 📜 Standing Rules for Objective 3
1.  **Memory Guard:** Always check RAM (`free -h`) before training. If `< 2GB` free, use `Lite` scripts.
2.  **Traceability:** MGF files MUST keep the original RAW filename prefix (e.g., `Patient1_SampleA.mgf`).
3.  **No Placeholders:** All code must be functional; if a dataset is too large, use stratified sampling rather than empty placeholders.
4.  **Audit Logs:** All background processes must output to `/logs/` with a unique timestamp.

---

## 📌 Maintenance & Next Steps
*This manual is a living document. The current focus is validating the **Medium-Lite model** (1,000 samples) and implementing the **Proteome Comparison** phase to filter non-human (neoepitope) hits.*
