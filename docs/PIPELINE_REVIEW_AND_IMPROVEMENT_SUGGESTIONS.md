# Pipeline Review and Improvement Suggestions

## Purpose
This document records findings from inspection of the Objective 3 pipeline and identifies areas requiring review, validation, or future improvement.

## High Priority Findings

### 1. Mass Filter Not Integrated Into Main Pipeline
Status: Review Required

The repository contains `src/cnnlstm/mass_filter.py`, which performs precursor mass validation of de novo predictions.

Current observed flow:

06_predict_denovo.py -> de_novo_candidates.tsv -> 07_filter_neoantigens.py

Recommended flow:

06_predict_denovo.py -> de_novo_candidates.tsv -> mass_filter.py -> de_novo_candidates_massfiltered.tsv -> 07_filter_neoantigens.py

Reason:
- Removes physically inconsistent peptide predictions.
- Likely reduces false positives.
- Provides additional confidence layer.

Action:
- Integrate mass filtering into run_pipeline.sh or inference workflow.
- Compare candidate counts before and after filtering.

---

### 2. Silent Missing Spectrum Handling
Status: Review Required

Observed in spectral_dataset.py.

When a spectrum cannot be located, the dataset currently returns a zero-filled spectrum tensor.

Risk:
- Missing scans become valid training examples.
- Can silently reduce model quality.
- Dataset integrity issues become difficult to detect.

Recommendation:
- Track missing spectrum counts.
- Log affected scans.
- Fail if missing rate exceeds threshold.

---

### 3. De Novo FDR Strategy Requires Validation
Status: Scientific Review Required

Current inference uses reversed peptide sequences as decoys.

Potential concerns:
- Reversed predictions are not independently generated decoys.
- Target-decoy estimates may not be publication-grade.

Recommendation:
- Validate against established de novo FDR approaches.
- Compare with Casanovo or DeepNovo methodologies.

---

### 4. RNA Evidence Currently Debug Only
Status: Expected Limitation

Manifest rebuild assigns RNA source as mock_debug.

Impact:
- Expression evidence cannot support biological validation.
- Final candidates should not be described as expression-supported neoantigens.

Recommendation:
- Clearly label outputs as lacking patient-matched RNA evidence.
- Prioritize acquisition of real patient RNA datasets.

---

## Architecture Validation Tasks

### CNN-LSTM Architecture Audit
Files:
- src/cnnlstm/cnnlstm_model.py
- src/training/05_train_denovo_model.py
- src/inference/06_predict_denovo.py

Questions:
- Is decoder architecture identical between training and inference?
- Are token definitions identical?
- Are sequence length assumptions consistent?
- Is teacher forcing behavior documented?

---

### Dataset Integrity Audit
Files:
- src/cnnlstm/psm_dataset.py
- src/cnnlstm/spectral_dataset.py

Review:
- Scan mapping accuracy.
- Duplicate spectrum handling.
- Missing scan handling.
- MGF indexing performance.

---

## Data Quality Checks

### Manifest Consistency
Review:
- patient_id consistency.
- validation_id consistency.
- HLA provenance.
- cohort assignments.

### Expression Data
Review:
- Mock TPM generation labeling.
- Expression lookup correctness.
- UniProt identifier normalization.

---

## Performance Opportunities

### MGF Parsing
Potential improvements:
- Spectrum indexing cache.
- Persistent lookup tables.
- Faster scan retrieval.

### Training
Potential improvements:
- Learning rate scheduling review.
- Early stopping validation.
- Class imbalance analysis.
- Train/test leakage audit.

---

## Documentation To Create

Recommended future docs:

1. docs/code_reference/training.md
2. docs/code_reference/inference.md
3. docs/code_reference/postprocessing.md
4. docs/code_reference/data_prep.md
5. docs/scientific_assumptions.md
6. docs/model_limitations.md
7. docs/publication_readiness_checklist.md

---

## Publication Readiness Checklist

Before claiming biological neoantigen discovery:

- Verify patient-specific mutation source.
- Verify patient-matched RNA evidence.
- Verify HLA assignment provenance.
- Validate de novo FDR methodology.
- Validate mass consistency filtering.
- Confirm reference proteome filtering.
- Confirm ranking methodology.

---

## Overall Assessment

Strengths:
- End-to-end Objective 3 pipeline exists.
- Manifest system is substantially improved.
- HLA provenance tracking exists.
- RNA evidence gating exists.
- Mass filter implementation exists.

Remaining Review Areas:
- CNN-LSTM architecture audit.
- FDR methodology audit.
- Mass filter integration.
- Dataset integrity validation.
- Publication-grade evidence review.
