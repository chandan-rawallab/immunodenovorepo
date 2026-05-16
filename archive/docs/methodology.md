# Methodology: Hybrid Immunopeptidomics Pipeline for Neoepitope Discovery

## 1. Overview

We present a hybrid computational pipeline integrating database-driven peptide identification with deep learning-based de novo sequencing to discover candidate neoepitopes from liquid chromatography-tandem mass spectrometry (LC-MS/MS) immunopeptidomics data. The pipeline combines MaxQuant database search against the canonical human proteome with a convolutional neural network-long short-term memory (CNN-LSTM) sequence-to-sequence model to identify peptides absent from reference databases, followed by multi-evidence ranking for neoantigen prioritisation.

## 2. Dataset

Raw LC-MS/MS data were obtained from the PRIDE archive (accession PXD005231) [1]. The dataset comprises HLA class I immunopeptidomics measurements from tumour samples, acquired on Thermo Fisher Orbitrap instruments. A total of 20 raw files were included in the analysis.

## 3. Data Preprocessing

Proprietary Thermo Fisher raw files were converted to the open Mascot Generic Format (MGF) using ThermoRawFileParser v2.0.0-dev [2], a native Linux implementation that does not require Mono or Wine. Conversion was performed in batch: each `.raw` file produced a corresponding `.mgf` file containing the extracted MS/MS peak lists with mass-to-charge (m/z) values and associated intensities.

## 4. Database Search (Immunopeptidome Construction)

### 4.1 Search Engine

Database searching was performed using MaxQuant v2.8.0.0 [3], a widely used proteomics search engine. MaxQuant was executed via its command-line interface (`MaxQuantCmd.dll`) using .NET 8.0 runtime on Linux.

### 4.2 Reference Database

Spectra were searched against the UniProt/Swiss-Prot human reference proteome (reviewed entries only, downloaded April 2026). A reversed-sequence decoy database was used for false discovery rate (FDR) estimation.

### 4.3 Search Parameters

Key MaxQuant parameters were configured as follows:

- **Enzyme**: Trypsin/P (cleavage C-terminal to Lysine and Arginine, even before Proline) with a maximum of 2 missed cleavages
- **Variable modifications**: Oxidation of methionine; Acetylation of protein N-termini
- **Fixed modifications**: None
- **Precursor mass tolerance**: 20 ppm for first search, 4.5 ppm for main search
- **Fragment mass tolerance**: 20 ppm (FTMS), 0.5 Da (ITMS)
- **Minimum peptide length**: 7 amino acids
- **Maximum peptide mass**: 4,600 Da
- **PSM FDR**: 1% (peptide, protein, and site level)
- **Quantification mode**: Label-free quantification (LFQ) disabled for immunopeptidomics
- **Number of threads**: 1 (sequential processing to prevent out-of-memory errors on large datasets)

### 4.4 Sequential Search Strategy

Raw files were processed sequentially using a custom orchestrator script (`scale_search.py`). Each file was searched independently with a dynamically generated parameter file. For each raw file, the MaxQuant output (`msms.txt`) was extracted and stored as `msms_{filename}.txt` in a dedicated production results directory. This approach avoided memory exhaustion from concurrent processing of multiple large raw files.

## 5. Deep Learning Model Architecture

We developed a custom sequence-to-sequence model (`NeoepitopeSeq2Seq`) that translates binned mass spectra into amino acid sequences de novo. The model was implemented in PyTorch and comprises three main components:

### 5.1 Spectral Encoding (CNN)

The input spectrum is binned into a fixed-length vector of 20,000 bins spanning m/z 0-2,000 Da with a bin width of 0.1 Da and normalised to unit maximum intensity. This vector is processed by a 1D convolutional neural network:

- **Layer 1**: Conv1d(1, 64, kernel=5, stride=2, padding=2), BatchNorm1d, ReLU, MaxPool1d(k=2)
- **Layer 2**: Conv1d(64, 128, kernel=5, stride=2, padding=2), BatchNorm1d, ReLU, MaxPool1d(k=2)

The CNN reduces the spectral dimensionality from 20,000 to 1,250 feature maps, which are then adaptively pooled to exactly 30 time steps via adaptive average pooling.

### 5.2 Sequence Decoding (BiLSTM)

The temporally-structured CNN features are processed by a bidirectional LSTM:

- **Input size**: 128
- **Hidden size**: 128 per direction
- **Number of layers**: 2
- **Directionality**: Bidirectional

The bidirectional LSTM captures both forward and backward dependencies in the spectral sequence, analogous to how b- and y-ion series provide complementary information in tandem mass spectrometry.

### 5.3 Output Layer

The bidirectional LSTM outputs are concatenated and passed through a dropout layer (p=0.3) followed by a linear projection to the vocabulary size (23 tokens: 20 standard amino acids plus special tokens for padding, start-of-sequence, and end-of-sequence).

### 5.4 Vocabulary

The model operates over a vocabulary of 23 tokens: the 20 canonical amino acids (A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y) plus `<PAD>` (index 0), `<START>` (index 21), and `<END>` (index 22). Sequences are encoded with `<START>` and `<END>` sentinels and padded or truncated to a maximum length of 30 tokens.

## 6. Training Procedure

### 6.1 Dataset Construction

Training data consisted of paired MGF spectra and MaxQuant-derived peptide sequences. Only spectra with confident database identifications were used. To filter peptide-spectrum matches (PSMs), we required that each spectrum have a corresponding MaxQuant identification. The dataset was constructed across all available .mgf files by matching each spectrum scan number to its corresponding entry in the MaxQuant msms.txt output.

### 6.2 Data Split

The labelled dataset was randomly split into training (80%) and validation (20%) sets.

### 6.3 Hyperparameters

- **Optimiser**: AdamW
- **Learning rate**: 0.001
- **Batch size**: 32
- **Number of epochs**: 100
- **Loss function**: Cross-entropy loss (ignoring padding index)
- **Hardware**: CPU (CUDA-compatible GPU not available on the analysis server)

### 6.4 Checkpointing

Model weights were saved every 5 epochs during training to enable evaluation at multiple training stages and to provide fallback checkpoints in case of interruption.

## 7. De Novo Prediction

After training, the model performs inference on all MGF files (regardless of whether they have database search results). For each spectrum, the model generates a predicted amino acid sequence by greedy decoding: at each time step, the token with the highest softmax probability is selected, and decoding continues until the `<END>` token is produced or the maximum sequence length is reached. The confidence score for each prediction is computed as the mean per-residue probability.

Predictions are filtered to retain only canonical HLA class I peptide lengths (8-11 amino acids) with sequences matching the pattern of standard amino acids.

## 8. Database Subtraction and Candidate Ranking

### 8.1 Subtraction Logic

Candidate neoepitopes are identified by subtracting all database-identified peptides from the de novo predictions. Specifically, a de novo-predicted peptide is retained as a candidate neoepitope if and only if:

1. Its sequence does not appear in the MaxQuant identifications for the same sample at a PSM-level FDR threshold of 1%
2. The de novo prediction confidence score exceeds a minimum cutoff (default: 0.70)

This subtraction ensures that only peptides absent from the reference proteome are considered as potential neoepitopes.

### 8.2 Multi-Evidence Ranking

Retained candidates are ranked using a composite scoring system incorporating:

1. **HLA binding affinity**: Predicted binding affinity to patient-specific HLA alleles (using MHCflurry or pre-computed binding predictions)
2. **Gene expression**: RNA-seq expression evidence at the gene level (TPM)
3. **Variant evidence**: Known somatic mutation information from whole-exome or targeted sequencing
4. **PSM support**: Number of spectra supporting each de novo prediction

Candidates are classified into evidence tiers:
- **Class A**: Supported by mutation evidence, binding affinity (rank ≤ 2%), and expression (TPM ≥ 1)
- **Class B**: Supported by binding and expression but lacking direct mutation evidence
- **Class C**: Supported by fewer than two orthogonal data types

Ranking within each evidence class is performed by binding rank (ascending), expression level (descending), and de novo confidence score (descending).

## 9. Implementation Details

### 9.1 Software Environment

The pipeline was implemented and executed on Ubuntu Linux (22.04). Key software dependencies included:

- MaxQuant v2.8.0.0 (via .NET 8.0 runtime)
- ThermoRawFileParser v2.0.0-dev
- Python 3.10+ with PyTorch, NumPy, pandas, pyteomics
- Custom Python modules in `src/dl_module/` (model, dataset, training) and `src/objective3/` (pipeline orchestration, CLI)

### 9.2 Code Organisation

The source code is organised as follows:

- `src/dl_module/`: Deep learning model definition (`model.py`), dataset classes (`dataset.py`, `production_dataset.py`), training script (`train.py`), and prediction script (`predict.py`)
- `src/objective3/`: Pipeline orchestration (`pipeline.py`), command-line interface (`cli.py`), MGF parsing utilities (`mgf_utils.py`), and I/O helpers (`io_utils.py`)
- `src/configs/`: MaxQuant parameter templates (`mqpar_production.xml`)
- `src/scale_search.py`: Sequential MaxQuant search orchestrator

## 10. Limitations

Several limitations should be acknowledged. First, the deep learning model was trained exclusively on CPU due to the absence of a CUDA-compatible GPU on the analysis server, limiting training speed and the feasibility of hyperparameter optimisation. Second, the MaxQuant database search covered only 11 of 20 raw files at the time of analysis, meaning the immunopeptidome reference may be incomplete. Third, de novo predictions are subject to higher false discovery rates than database searches, necessitating stringent confidence thresholds and orthogonal validation.

## References

[1] Bassani-Sternberg, M., et al. (2016). Mass spectrometry of HLA-I peptidomes reveals strong effects of protein abundance and turnover on antigen presentation. *Molecular & Cellular Proteomics*.

[2] Hulstaert, N., et al. (2020). ThermoRawFileParser: modular, scalable and cross-platform RAW file conversion. *Journal of Proteome Research*, 19(1), 537-542.

[3] Cox, J. & Mann, M. (2008). MaxQuant enables high peptide identification rates, individualized p.p.b.-range mass accuracies and proteome-wide protein quantification. *Nature Biotechnology*, 26(12), 1367-1372.

[4] Tran, N.H., et al. (2017). De novo peptide sequencing by deep learning. *Proceedings of the National Academy of Sciences*, 114(31), 8247-8252.

[5] Tran, N.H., et al. (2019). Deep learning enables de novo peptide sequencing from data-independent-acquisition mass spectrometry. *Nature Methods*, 16(1), 63-66.
