# Detailed Code Audit

## Scope
This document is intended to become the permanent engineering and scientific audit record for the Objective 3 ImmunoDeNovo pipeline.

Repository Areas:
- src/data_prep
- src/cnnlstm
- src/training
- src/inference
- src/postprocess
- src/evaluation
- src/validation

---

# Executive Summary

Current strengths:
- End-to-end pipeline exists.
- Curated manifest workflow implemented.
- HLA provenance tracking present.
- Expression linkage framework implemented.
- Training and inference separation is clear.
- Validation and provenance auditing modules exist.

Primary review priorities:
1. Model architecture validation.
2. Dataset integrity validation.
3. De novo FDR methodology validation.
4. Mass-filter integration.
5. Publication-readiness review.

---

# Module Audit Template

For each script evaluate:

## Metadata
- Purpose
- Inputs
- Outputs
- Dependencies

## Engineering Review
- Error handling
- Logging
- Scalability
- Memory efficiency
- Reproducibility

## Scientific Review
- Assumptions
- Validation evidence
- Potential bias
- Publication risks

## Risk Rating
- Low
- Medium
- High
- Critical

---

# Data Preparation Audit

## 00b_build_manifest.py

Purpose:
Build study manifest and patient mapping.

Strengths:
- HLA normalization logic.
- SDRF support.
- Manual override support.
- Multiple fallback strategies.

Review Items:
- Verify filename heuristics.
- Verify patient extraction rules.
- Verify SDRF column detection.

Risk:
Low.

---

## 00d_link_expression.py

Purpose:
Link RNA evidence.

Strengths:
- Explicit debug mode.
- Reproducible mock generation.
- Manifest integration.

Review Items:
- Ensure mock data never used in scientific claims.
- Verify RNA provenance labels.

Risk:
Medium.

---

## 00f_rebuild_curated_manifest.py

Purpose:
Rebuild curated publication cohort.

Strengths:
- Controlled cohort definition.
- Explicit exclusion tracking.

Review Items:
- Validate all patient/HLA mappings against publication.
- Verify excluded runs remain documented.

Risk:
Medium.

---

## 02_convert_raw_to_mgf.py

Review Items:
- Conversion failure recovery.
- MGF completeness verification.
- Spectrum count validation.

Risk:
Medium.

---

## 05_extract_unlabeled_spectra.py

Strengths:
- Streaming processing.
- Manifest filtering.

Review Items:
- Scan identifier consistency.
- Duplicate handling.

Risk:
Medium.

---

# CNN-LSTM Audit

## psm_dataset.py

Strengths:
- Lazy loading approach.
- Spectrum binning.

Review Items:
- Dataset memory footprint.
- Sequence truncation effects.
- Scan mapping correctness.

Risk:
High.

---

## spectral_dataset.py

Review Items:
- Missing spectrum behavior.
- Zero-vector fallback validation.
- Scan lookup robustness.

Risk:
High.

---

## cnnlstm_model.py

Review Items:
- Architecture documentation.
- Decoder correctness.
- Token vocabulary consistency.
- Output sequence handling.

Risk:
High.

---

## mass_filter.py

Strengths:
- Physically motivated validation.
- Independent filtering layer.

Review Items:
- Verify spectrum title matching.
- Benchmark candidate reduction.
- Integrate into production workflow.

Risk:
Low.

Status:
Not yet confirmed as part of main pipeline.

---

# Training Audit

## 05_train_denovo_model.py

Review Items:
- Train/test leakage.
- Split strategy.
- Early stopping.
- Checkpoint selection.
- Metric reporting.

Risk:
Critical.

---

# Inference Audit

## 06_predict_denovo.py

Review Items:
- Beam search implementation.
- Confidence score calibration.
- Decoy generation.
- Candidate export consistency.

Risk:
Critical.

---

# Postprocessing Audit

## 07_filter_neoantigens.py

Review Items:
- Reference proteome completeness.
- Exact peptide matching logic.
- False positive controls.

Risk:
Critical.

---

## 08_rank_candidates.py

Review Items:
- Ranking criteria.
- HLA evidence weighting.
- Expression evidence weighting.

Risk:
High.

---

# Evaluation Audit

## 09_evaluate_neoantigens.py

Review Items:
- Biological interpretation limits.
- Ground truth availability.

Risk:
High.

---

## 10_evaluate_denovo_model.py

Review Items:
- Accuracy metrics.
- Sequence-level accuracy.
- Spectrum-level accuracy.

Risk:
High.

---

# Validation Audit

## preflight_validate.py

Purpose:
Pipeline readiness checks.

Recommendation:
Expand to verify:
- Manifest consistency.
- MGF/PSM parity.
- Expression provenance.
- Model checkpoint availability.

---

## provenance_audit.py

Purpose:
Scientific traceability.

Recommendation:
Require provenance report generation before final ranking.

---

# Publication Readiness Gates

Gate 1:
Manifest validated.

Gate 2:
HLA provenance validated.

Gate 3:
MGF/PSM mapping validated.

Gate 4:
Mass filter integrated.

Gate 5:
FDR methodology validated.

Gate 6:
Reference proteome filtering validated.

Gate 7:
Expression evidence provenance validated.

Gate 8:
Final scientific review completed.

---

# Next Actions

1. Perform line-by-line audit of cnnlstm_model.py.
2. Perform line-by-line audit of 05_train_denovo_model.py.
3. Perform line-by-line audit of 06_predict_denovo.py.
4. Verify mass_filter integration path.
5. Verify ranked_neoantigens.tsv provenance chain.
6. Produce publication-readiness report.
