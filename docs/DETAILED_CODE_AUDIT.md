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

# Confirmed Audit Findings (Completed Reviews)

## psm_dataset.py

Confirmed:
- Dataset is NOT truly lazy despite comments claiming lazy parsing.
- Entire MGF file is iterated during initialization and matching spectra are stored in self.spectra.
- Memory consumption scales with number of labeled spectra.
- Scan→sequence mapping logic appears correct.
- Sequence truncation occurs through max_seq_len clipping and requires validation against peptide length distribution.

Recommendation:
- Replace in-memory spectrum storage with indexed access or streaming retrieval.

Risk:
High.

---

## 00d_link_expression.py

Confirmed:
- Mock RNA generation only occurs when --debug-expression is explicitly enabled.
- Deterministic RNG seeding improves reproducibility.
- RNA provenance labels are recorded.
- Real expression overrides are supported.

Recommendation:
- Prevent mock_debug RNA from entering publication-facing outputs.

Risk:
Medium.

---

## 00f_rebuild_curated_manifest.py

Confirmed:
- Rebuilds cohort from curated publication manifest.
- Tracks excluded runs in a separate audit table.
- Enforces expected cohort size checks.
- Improves reproducibility versus heuristic manifest construction.

Recommendation:
- Validate all curated HLA assignments against publication supplementary material.

Risk:
Medium.

---

## 05_extract_unlabeled_spectra.py

Confirmed:
- Uses streaming MGF read/write workflow.
- Avoids loading entire MGF collections into memory.
- Manifest filtering implemented.
- Stale output cleanup implemented.

Recommendation:
- Verify scan identifier extraction against all MGF title formats.

Risk:
Low-Medium.

---

## mass_filter.py

Confirmed:
- Implements physically meaningful precursor-mass validation.
- Computes theoretical peptide mass.
- Computes precursor neutral mass.
- Filters candidates using ppm tolerance.
- Keeps candidates conservatively when precursor metadata is unavailable.

Outstanding Review:
- Verify production integration path from 06_predict_denovo.py.
- Benchmark candidate reduction and recovery statistics.

Risk:
Low.

---

# Remaining Audit Sections

(Existing audit sections retained below and will be progressively replaced with confirmed findings as modules are inspected.)