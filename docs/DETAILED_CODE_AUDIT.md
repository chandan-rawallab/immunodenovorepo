# Detailed Code Audit

## Scope
This document is intended to become the permanent engineering and scientific audit record for the Objective 3 ImmunoDeNovo pipeline.

## Confirmed Audit Findings (Completed Reviews)

### cnnlstm_model.py
- CNN + BiLSTM architecture is reasonable for spectrum-to-sequence learning.
- BatchNorm and Dropout are present.
- No attention mechanism currently implemented.
- Model itself is not mass-aware; relies on downstream filtering.

Recommendations:
- Evaluate attention-based decoder.
- Integrate precursor-mass awareness or enforce mass filtering during inference.

Risk: Medium.

---

### psm_dataset.py
- Dataset is NOT truly lazy despite comments claiming lazy parsing.
- Entire MGF file is iterated during initialization and matching spectra are stored in memory.
- Duplicate peptide observations may contribute to train/test leakage if split incorrectly.

Recommendations:
- Replace in-memory spectrum storage with indexed access.
- Audit peptide duplication across train/test partitions.

Risk: High.

---

### 05_train_denovo_model.py
- PAD masking, gradient clipping, early stopping and reproducible splits are implemented.
- Dedicated held-out test set is generated.
- Potential train/test leakage exists because splitting is spectrum-based rather than peptide-based.
- Checkpoint provenance tracking could be improved.

Recommendations:
- Use peptide-grouped train/test splitting.
- Save training_config.json and manifest metadata with checkpoints.
- Verify final held-out evaluation path.

Risk: High (scientific validation).

---

### 00d_link_expression.py
- Mock RNA generation only occurs when --debug-expression is explicitly enabled.
- Deterministic RNG seeding improves reproducibility.

Risk: Medium.

---

### 00f_rebuild_curated_manifest.py
- Rebuilds cohort from curated publication manifest.
- Tracks excluded runs separately.

Risk: Medium.

---

### 05_extract_unlabeled_spectra.py
- Uses streaming MGF read/write workflow.
- Avoids loading complete MGF collections into memory.

Risk: Low-Medium.

---

### mass_filter.py
- Implements precursor-mass validation.
- Production integration path remains to be verified.

Risk: Low.
