# Pipeline File Structure & Reference Card

A concise reference for every source file currently active in the Objective 3 neoantigen discovery pipeline, organized by module. Each entry covers the exact inputs consumed, outputs produced, and what the file does.

---

## Project Tree

```
experiments/
│
├── configs/
│   ├── sample_manifest.tsv               # Master sample registry (31 active runs, 9 patients)
│   ├── excluded_runs.tsv                 # Runs removed during curation
│   └── mqpar.xml / mqpar_production.xml  # MaxQuant search parameter configurations
│
├── data/
│   ├── mgf/                              # Labeled MS2 spectra (database-identified runs)
│   ├── mgf_unlabeled/                    # Dark-matter spectra (1,437,700 unidentified spectra)
│   ├── psms/                             # Raw MaxQuant msms.txt output (one file per run)
│   ├── reference/
│   │   └── uniprot_human_reviewed.only_human.fasta  # Filtered human UniProt proteome
│   └── expression/
│       └── <patient>_tpm.tsv             # Per-patient UniProt → TPM tables (CCLE surrogate)
│
├── results/
│   ├── immunopeptidome_psms.tsv          # Consolidated high-confidence WT PSMs (4,901 rows)
│   ├── de_novo_candidates.tsv            # Raw model predictions (1,375,140 rows)
│   ├── filtered_neoantigens.tsv          # Post-filter missense candidates (~591 missense / ~2,697 total)
│   ├── ranked_neoantigens.tsv            # HLA-ranked, expression-annotated final candidates
│   ├── checkpoints_curated31_v2/
│   │   ├── neoepitope_production_best.pth  # Active CNN-LSTM model weights
│   │   └── test_set_psms.tsv               # Held-out test PSMs (saved at split time)
│   ├── model_accuracy_curated31_v2.json  # CNN-LSTM evaluation metrics (25.1% exact accuracy)
│   ├── casanovo_baseline.log             # Casanovo unlabeled sequencing run log
│   ├── provenance_audit_current.md       # Data lineage audit report
│   └── preflight_report.md              # Cohort readiness check report
│
└── src/
    ├── data_prep/    # Phase 1 & 2: Data ingestion and clinical annotation
    ├── cnnlstm/      # Model architecture definition and dataset loaders
    ├── training/     # Model training orchestration
    ├── inference/    # De novo prediction engine
    ├── evaluation/   # Accuracy benchmarking and clinical validation
    ├── postprocess/  # Biological filtering and candidate ranking
    └── validation/   # Preflight integrity checks and provenance auditing
```

---

## Module: `src/data_prep/` — Data Ingestion & Annotation

---

### `00_acquire_data.py`
- **Takes**: PRIDE accession ID (`PXD005231`), local target directories for raw and PSM files.
- **Outputs**: `data/raw/*.raw` (Thermo binary), `data/psms/*_msms.txt` (MaxQuant results).
- **Work**: Queries the PRIDE EBI REST API to enumerate all files for the accession. Downloads raw instrument files and MaxQuant result tables in parallel threads with MD5 checksum validation. Skips already-present files to save bandwidth. This is the cold-start entry point for the entire pipeline.

---

### `01_download_pride.py`
- **Takes**: A PRIDE FTP/HTTPS URL and target directory path.
- **Outputs**: Downloaded files written to the target directory.
- **Work**: Helper module invoked by `00_acquire_data.py`. Handles the low-level FTP/HTTPS transfer logic with retry support. Not designed for standalone invocation.

---

### `00b_build_manifest.py`
- **Takes**: `data/raw/` and `data/psms/` directory contents.
- **Outputs**: `configs/sample_manifest.tsv` (initial draft with `TBD` placeholders for HLA alleles).
- **Work**: Scans directory filenames to extract patient subject codes and instrument run identifiers. Assembles an initial structured registry mapping each run to a patient. Initializes placeholder columns for HLA alleles, expression file paths, and cohort roles.

---

### `00f_rebuild_curated_manifest.py`
- **Takes**: `configs/sample_manifest.tsv.bak` (original draft manifest), known bad run ID list.
- **Outputs**: Final `configs/sample_manifest.tsv` (31 curated active runs), `configs/excluded_runs.tsv`.
- **Work**: Filters out unverified, mislabeled, or empty runs from the initial manifest draft. Normalizes patient ID inconsistencies (e.g., `CM647` → `CM467`). Acts as the final QC gate before any heavy computation to ensure the cohort is clean.

---

### `00c_autotype_hla.py`
- **Takes**: `results/immunopeptidome_psms.tsv`, `configs/sample_manifest.tsv`.
- **Outputs**: Updated `configs/sample_manifest.tsv` with the `hla_alleles` column populated.
- **Work**: For patients missing clinical NGS HLA typing, runs MHCflurry binding predictions against a library of common HLA-A, B, and C alleles using that patient's own WT peptides from the immunopeptidome. Selects alleles with the highest binding density (percentile rank ≤ 2%). Curated publication-sourced alleles take priority over computed ones.

---

### `00d_link_expression.py`
- **Takes**: `configs/sample_manifest.tsv`, human reference FASTA.
- **Outputs**: Per-patient `data/expression/<patient>_tpm.tsv`, updated manifest `rna_expr_path` column.
- **Work**: Connects RNA expression profiles to each patient. If real RNA-seq is unavailable, generates a mock log-normal TPM profile (10% genes set to high expression) as a diagnostic placeholder. Updates the manifest to point each sample at its expression file path.

---

### `00e_map_ccle_to_uniprot.py`
- **Takes**: `CCLE_expression.csv` (DepMap), `Model.csv` (CCLE sample metadata), `HUMAN_9606_idmapping_selected.tab.gz` (UniProt ID mapping), `configs/sample_manifest.tsv`.
- **Outputs**: Per-patient `data/expression/<patient>_tpm.tsv` with UniProt accession IDs as the gene key.
- **Work**: Replaces mock expression profiles with real CCLE surrogate profiles. Maps CCLE Entrez gene symbols to UniProt accessions via the official UniProt ID mapping table. Selects tissue-appropriate CCLE cell lines (LCL lines for B-cell patients; skin/melanoma lines for TIL patients) and averages TPM values across matched lines. Strips whitespace from manifest values to prevent `KeyError` mismatches on patient ID lookup.

> ⚠️ **Limitation**: All current expression values are CCLE cell-line surrogates, not patient-matched RNA. The `rna_source` column in the manifest is set to `ccle_lcl_surrogate` for all 31 runs.

---

### `00e_filter_human_fasta.py`
- **Takes**: `data/reference/uniprot_human_reviewed.fasta` (may contain BSA/trypsin contaminants).
- **Outputs**: `data/reference/uniprot_human_reviewed.only_human.fasta`.
- **Work**: Reads the UniProt FASTA and retains only entries with `OS=Homo sapiens` or `OX=9606` in their header. Removes common lab contaminants (bovine serum albumin, porcine trypsin) to prevent false protein alignments during downstream mutation detection.

---

### `02_convert_raw_to_mgf.py`
- **Takes**: `data/raw/*.raw` (Thermo Fisher proprietary binary format).
- **Outputs**: `data/mgf/<run_id>.mgf` (Mascot Generic Format text peak lists).
- **Work**: Shells out to `ThermoRawFileParser` (via .NET 8) to convert proprietary binary files into open-format text. Filters to MS2-level scans only. Each output MGF file contains `BEGIN IONS` blocks with `SCANS=`, `PEPMASS=`, `CHARGE=`, and per-peak `m/z intensity` pairs.

---

### `04_extract_psms.py`
- **Takes**: `data/psms/*_msms.txt` (MaxQuant/Andromeda output), `configs/sample_manifest.tsv`.
- **Outputs**: `results/immunopeptidome_psms.tsv` (4,901 rows: the patient immunopeptidome).
- **Work**: Reads MaxQuant PSM files for all 31 curated runs. Applies three confidence filters: Andromeda PEP ≤ 0.01 (1% local FDR), peptide length 8–11 amino acids (HLA-I canonical range), and no decoy hits. Concatenates all runs into a single consolidated immunopeptidome table that forms the training ground-truth and the subtraction baseline.

---

### `05_extract_unlabeled_spectra.py`
- **Takes**: `data/mgf/*.mgf` (all spectra), `results/immunopeptidome_psms.tsv` (identified scan IDs), manifest.
- **Outputs**: `data/mgf_unlabeled/<run_id>.mgf` — 1,437,700 spectra in total across all 31 runs.
- **Work**: Performs spectral subtraction. For each run, builds the set of scan IDs successfully identified by MaxQuant. Reads the corresponding MGF and writes only the **unmatched** scans to the unlabeled output directory. The result is the "dark matter" immunopeptidome — spectra the database search could not explain — which feeds the de novo discovery step.

---

## Module: `src/cnnlstm/` — Model Definition & Utilities

---

### `cnnlstm_model.py`
- **Takes**: Binned spectrum tensor `[batch, 20000]`, target peptide token sequence.
- **Outputs**: Per-position log-probability distributions over the 23-token vocabulary `[batch, max_len, 23]`.
- **Work**: Defines the `NeoepitopeSeq2Seq` architecture. A 1D CNN encoder extracts local fragmentation ladder features. Adaptive average pooling compresses features to 30 time steps. A 2-layer Bidirectional LSTM decoder produces the sequence. A 2-layer FC classification head maps to vocabulary logits. Current config: `cnn_channels=128`, `lstm_hidden=256`.

---

### `psm_dataset.py`
- **Takes**: `results/immunopeptidome_psms.tsv`, `data/mgf/` directory.
- **Outputs**: `torch.Dataset` yielding `(spectrum_tensor [20000], peptide_token_tensor)` pairs.
- **Work**: Implements `MgfDataset`. Loads PSM records, matches each to its spectrum by scan ID, bins raw m/z peaks into a fixed 20,000-bin vector (0–2000 Da, 0.1 Da resolution, intensity normalized to base peak), and encodes the peptide string as integer token indices. Used exclusively for **labeled** training and evaluation.

---

### `spectral_dataset.py`
- **Takes**: `data/mgf_unlabeled/` directory.
- **Outputs**: `torch.Dataset` yielding `(spectrum_tensor [20000],)` — no labels.
- **Work**: Implements `SpectralDataset` for inference on unlabeled spectra. Mirrors the exact same binning and normalization logic as `psm_dataset.py` but reads from unlabeled MGFs without requiring ground-truth peptide labels. Falls back to the local `mgf_utils.IndexedMgfFallback` reader if `pyteomics` is unavailable.

---

### `mgf_utils.py`
- **Takes**: Path to any MGF file.
- **Outputs**: `MgfSpectrum` dataclass objects (spectrum ID, m/z array, intensity array, params dict).
- **Work**: Provides lightweight MGF parsing independent of the `pyteomics` library. `IndexedMgfFallback` scans an MGF file once to build a byte-offset index keyed by scan ID, enabling fast random-access reads without loading the full file into memory. Critical for environments where `pyteomics` is not installed.

---

### `mass_filter.py`
- **Takes**: `results/de_novo_candidates.tsv`, `data/mgf/` directory, ppm tolerance (default: 20 ppm).
- **Outputs**: Filtered candidates TSV retaining only mass-consistent predictions.
- **Work**: For each predicted peptide, calculates its theoretical monoisotopic mass from standard amino acid residue masses and compares it against the precursor m/z and charge state from the original spectrum. Rejects predictions where the mass error exceeds the ppm tolerance. Can be used as a standalone script or imported as `filter_dataframe()`.

---

## Module: `src/training/` — Model Training

---

### `05_train_denovo_model.py`
- **Takes**: `data/mgf/` (labeled spectra), `results/immunopeptidome_psms.tsv` (ground truth), `--checkpoint-dir`, `--epochs`, `--lr`.
- **Outputs**: `results/checkpoints_curated31_v2/neoepitope_production_best.pth` (best weights), per-epoch checkpoint files, `test_set_psms.tsv` (held-out split saved alongside).
- **Work**: Loads the PSM table and instantiates `SpectralDataset`. Performs a stratified 70/15/15 split by peptide length at the **PSM level** (⚠️ known issue: 62.4% of test peptide sequences also appear in training, inflating the reported 25.1% accuracy). Trains `NeoepitopeSeq2Seq` with Cross-Entropy loss and AdamW for 40 epochs with teacher forcing. Saves the checkpoint with the lowest validation loss.

---

### `test_train_split.py`
- **Takes**: Nothing.
- **Outputs**: Prints "test passed" to stdout.
- **Work**: A single-line import check that verifies `sklearn.model_selection.train_test_split` is available in the active Python environment. Not a functional pipeline component.

---

## Module: `src/inference/` — De Novo Prediction

---

### `06_predict_denovo.py`
- **Takes**: `--model` (trained `.pth` weights), `--mgf_dir` (`data/mgf_unlabeled/`), `--output`.
- **Outputs**: `results/de_novo_candidates.tsv` — columns: `run_id`, `spectrum_id`, `peptide`, `score`, `fdr` (1,375,140 rows).
- **Work**: The core discovery engine. Loads the trained CNN-LSTM checkpoint. Iterates over all unlabeled MGF files. For each spectrum, runs greedy argmax decoding through the LSTM decoder and records the predicted sequence and log-probability score. Applies a target-decoy FDR sweep using reversed spectrum vectors as decoys, filtering the final output at 5% FDR. Converts 1,437,700 input spectra into 1,375,140 raw candidate predictions.

---

### `convert_casanovo_output.py`
- **Takes**: `results/casanovo_all.mztab` (Casanovo mzTab output), `configs/sample_manifest.tsv`.
- **Outputs**: `results/de_novo_candidates_casanovo.tsv` (identical schema to `de_novo_candidates.tsv`).
- **Work**: Parses the PSM section of a Casanovo mzTab file. Maps `spectra_ref` identifiers back to `run_id` and `spectrum_id` using the manifest. Renames Casanovo confidence columns to the pipeline's `score` column. Enables the full filter and ranking chain to process Casanovo predictions identically to CNN-LSTM predictions for a fair head-to-head comparison.

---

## Module: `src/evaluation/` — Accuracy & Clinical Benchmarking

---

### `10_evaluate_denovo_model.py`
- **Takes**: `--checkpoint` (`.pth` weights), `--test-psms` (`test_set_psms.tsv`), `--mgf-dir` (labeled MGFs).
- **Outputs**: `results/model_accuracy_curated31_v2.json` — exact match accuracy, token accuracy, edit distance distribution, per-length breakdown.
- **Work**: Loads the trained model in inference mode. For each PSM in the held-out test split, looks up the raw MGF spectrum by scan ID, runs the greedy decoder, and compares the predicted peptide to the ground-truth sequence. Computes exact accuracy (with and without I/L mass-equivalence collapse), per-position token accuracy, and the fraction of predictions within edit distance ≤ 1.

---

### `09_evaluate_neoantigens.py`
- **Takes**: `results/ranked_neoantigens.tsv`, `data/reference/.../Dataset1.txt` (Bassani-Sternberg 2016 validated peptides), manifest.
- **Outputs**: `results/evaluation_report.md`.
- **Work**: Intersects the pipeline's top-ranked neoantigen candidates against peptides experimentally validated in the Bassani-Sternberg 2016 publication (the same PXD005231 dataset). Calculates Precision and Recall at Top-10, Top-25, and Top-50 cutoffs on a per-patient basis. Produces a markdown report of clinical validation performance.

---

## Module: `src/postprocess/` — Biological Filtering & Ranking

---

### `07_filter_neoantigens.py`
- **Takes**: `results/de_novo_candidates.tsv`, `results/immunopeptidome_psms.tsv`, `data/reference/uniprot_human_reviewed.only_human.fasta`, manifest.
- **Outputs**: `results/filtered_neoantigens.tsv` — adds columns: `wildtype_peptide`, `mutation_pos`, `wt_aa`, `mut_aa`, `mutation_type`, `source_protein`.
- **Work**: 5-stage biological filter cascade:
  1. **Score filter**: Drops predictions below `--score_cutoff` (default: -0.5).
  2. **Length filter**: Keeps only 8–11 AA sequences matching the HLA-I canonical pattern.
  3. **Database subtraction**: Removes peptides already found in the patient's WT immunopeptidome for that sample.
  4. **PSM support**: Requires ≥ 2 independent spectra supporting the same peptide.
  5. **Levenshtein-1 mutation check**: Pre-indexes all 8–11mer subsequences from the reference human proteome. For each candidate, searches exhaustively for a single-amino-acid substitution producing a known WT peptide. Records mutation position, WT amino acid, mutant amino acid, and source UniProt protein IDs. Excludes mutations at flanking positions (pos 1 or last) as likely cleavage artifacts.

---

### `08_rank_candidates.py`
- **Takes**: `results/filtered_neoantigens.tsv`, `configs/sample_manifest.tsv`, `data/expression/<patient>_tpm.tsv`.
- **Outputs**: `results/ranked_neoantigens.tsv` — adds columns: `binding_rank`, `expression_tpm`, `evidence_class`, `rna_source`.
- **Work**: For each filtered candidate, calls MHCflurry to predict HLA binding percentile rank against the patient's alleles from the manifest. Looks up expression by splitting `source_protein` (semicolon-separated UniProt IDs), querying the patient TPM file, and returning `max(TPM)` across all matching proteins. Groups candidates: **Class A** (missense + binding rank ≤ 2% + TPM ≥ 1.0), **Class B** (strong binder + expressed, non-missense), **Class C** (weak or unexpressed). Tags each row with `rna_source` to distinguish surrogate from patient-matched expression.

---

## Module: `src/validation/` — Integrity & Provenance

---

### `preflight_validate.py`
- **Takes**: Manifest, reference FASTA path, expected number of active runs.
- **Outputs**: `results/preflight_report.md`.
- **Work**: Runs a suite of pre-pipeline integrity checks: verifies every `run_id` has a matching PSM file and MGF file, checks for duplicate run IDs, validates HLA allele format strings, confirms expression file paths resolve on disk, and checks that the active run count matches the expected value. Flags all discrepancies in a human-readable markdown report.

---

### `provenance_audit.py`
- **Takes**: Manifest, local results directory, PRIDE REST API (live query).
- **Outputs**: `results/provenance_audit_current.md`, `results/provenance_audit_current.json`.
- **Work**: Audits data lineage end-to-end. Queries the PRIDE EBI API to confirm the expected raw file count for PXD005231. Checks that local PSM and MGF file counts match. Hashes model checkpoint weights to confirm no accidental overwrite. Records file sizes, modification timestamps, and row counts for all major pipeline outputs.

---

## Data Flow Summary

```
PRIDE (PXD005231)
       ↓
data/raw/*.raw + data/psms/*.txt
       ↓
04_extract_psms.py ──────────────────────────────────────────────┐
       ↓                                                          │
results/immunopeptidome_psms.tsv (4,901 WT PSMs)                 │
       ↓                         ↓                               │
02_convert_raw_to_mgf.py   05_extract_unlabeled_spectra.py       │
       ↓                         ↓                               │
data/mgf/*.mgf         data/mgf_unlabeled/*.mgf                  │
(labeled)               (1,437,700 dark-matter spectra)          │
       ↓                         ↓                               │
05_train_denovo_model.py  06_predict_denovo.py                   │
       ↓                         ↓                               │
checkpoints_curated31_v2/  de_novo_candidates.tsv (1,375,140)    │
neoepitope_production_best.pth   ↓                               │
       ↓               07_filter_neoantigens.py ←────────────────┘
10_evaluate_denovo_model.py      ↓                   (uses FASTA +
       ↓               filtered_neoantigens.tsv        WT PSMs)
model_accuracy_curated31_v2.json (~591 missense)
                                 ↓
                       08_rank_candidates.py ← expression/*_tpm.tsv
                                 ↓
                       ranked_neoantigens.tsv (~987 total)
                                 ↓
                       09_evaluate_neoantigens.py
                                 ↓
                       evaluation_report.md
```
