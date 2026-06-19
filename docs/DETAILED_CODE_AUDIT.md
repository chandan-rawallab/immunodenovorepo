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
- Raw peak intensities are used without transformation.
- No peak filtering is applied before vectorization.
- Unknown amino acids collapse into the PAD token.

Recommendations:
- Track and report missing-spectrum counts.
- Replace silent exception suppression with structured logging.
- Build a unified spectrum index instead of probing multiple key formats.
- Evaluate log/square-root intensity transformation.
- Evaluate top-N peak filtering.
- Add explicit UNK token handling.
- Add validation for peptide truncation.

Risk: High.

---

### mgf_utils.py
- Lightweight fallback implementation removes mandatory pyteomics dependency.
- MGF parser is streaming and memory efficient.
- IndexedMgfFallback builds a full in-memory index on first access.
- Duplicate spectrum identifiers silently overwrite previous entries.

Recommendations:
- Add duplicate-spectrum detection and reporting.
- Improve error messages for missing spectrum IDs.

Risk: Medium.

---

### 05_train_denovo_model.py
- Potential train/test leakage exists because splitting is spectrum-based rather than peptide-based.
- Dedicated held-out test set is generated.

Recommendations:
- Use peptide-grouped train/test splitting.
- Save training configuration and manifest provenance with checkpoints.

Risk: High.

---

### 06_predict_denovo.py
- Uses greedy decoding only.
- No beam-search implementation.
- FDR methodology should be independently validated.

Recommendations:
- Add checkpoint metadata validation.
- Evaluate beam-search decoding.
- Review target-decoy methodology and q-value estimation.

Risk: High.

---

### 07_filter_neoantigens.py
- Implements biologically sensible filtering.
- Repeated FASTA parsing may become expensive.
- Mutation detection is limited primarily to substitution-style inference.

Recommendations:
- Cache protein lookup information.
- Add audit metrics for candidate counts at every filter stage.

Risk: Medium-High.

---

### 08_rank_candidates.py
- Implements biologically meaningful ranking using HLA binding predictions plus RNA-expression evidence.
- RNA evidence gating is stronger than many prototype pipelines because debug/mock RNA sources are explicitly prevented from supporting biological evidence classes.
- Duplicate candidate collapsing prevents peptide inflation in downstream reports.
- MHCflurry execution is performed sample-by-sample and may become slow for larger cohorts.
- Expression lookup construction repeatedly iterates over RNA tables and builds a Python dictionary rather than using indexed joins.
- Extensive use of DataFrame apply(axis=1) may become a bottleneck at scale.
- Sample matching relies primarily on sample_id == patient_id and falls back to run_id matching.
- Missing MHCflurry predictions silently downgrade candidates.

Recommendations:
- Add ranking-stage audit metrics.
- Replace row-wise expression assignment with vectorized joins.
- Validate manifest sample matching assumptions.
- Add provenance fields for MHCflurry version and predictor configuration.

Risk: Medium-High.

---

### 10_evaluate_denovo_model.py
- Evaluation is performed against a dedicated held-out PSM test split.
- Metrics include exact peptide accuracy, token accuracy, edit distance, length accuracy, and amino-acid positional accuracy.
- Missing spectra are tracked and reported rather than silently ignored.
- Entire MGF files are streamed sequentially for scan lookup, which can become expensive on large cohorts.
- Evaluation recreates spectral vectorization logic instead of reusing a shared preprocessing module.
- Default model instantiation uses NeoepitopeSeq2Seq() without validating architecture metadata against the checkpoint.
- No confidence calibration, top-k accuracy, or beam-search evaluation metrics are reported.
- Duplicate peptide observations may inflate performance if train/test leakage exists upstream.

Recommendations:
- Save architecture metadata alongside checkpoints and validate before loading.
- Reuse a shared spectrum preprocessing implementation to prevent training/evaluation drift.
- Consider indexed MGF access for faster evaluation.
- Add top-k peptide recovery metrics and confidence calibration statistics.
- Explicitly audit peptide overlap between training and held-out sets.

Risk: Medium-High.