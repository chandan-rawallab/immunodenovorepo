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

### spectral_dataset.py
- Uses indexed MGF access and reader caching, which is positive for scalability.
- Missing spectra are silently replaced with zero-valued spectra and empty peptide labels.
- Broad exception handling suppresses root-cause debugging.
- Spectrum lookup relies on several fragile ID formats.
- Peptide truncation is not logged.
- Dataset initialization validates MGF existence using a row-wise loop that can be vectorized.

Recommendations:
- Track and report missing-spectrum counts.
- Replace silent exception suppression with structured logging.
- Build a unified spectrum index instead of probing multiple key formats.
- Log peptide truncation events.
- Vectorize MGF run validation.
- Evaluate DataLoader multiprocessing behaviour with cached readers.

Risk: High.

---

### mgf_utils.py
- Lightweight fallback implementation removes mandatory pyteomics dependency.
- MGF parser is streaming and memory efficient.
- IndexedMgfFallback builds a full in-memory index on first access, which can become expensive for very large MGF files.
- Duplicate spectrum identifiers silently overwrite previous entries in the index.
- get_by_id raises KeyError without diagnostic context.
- TITLE parsing only recognizes lowercase 'scan=' tokens and may miss alternative vendor formats.
- No validation exists for malformed MGF blocks.

Recommendations:
- Add duplicate-spectrum detection and reporting.
- Improve error messages for missing spectrum IDs.
- Support broader TITLE parsing patterns.
- Add parser statistics (spectra parsed, malformed blocks, duplicate IDs).
- Consider optional on-disk indexing for large production datasets.

Risk: Medium.

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