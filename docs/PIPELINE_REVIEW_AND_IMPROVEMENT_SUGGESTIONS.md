# Pipeline Review and Improvement Suggestions

## Purpose
This document records findings from inspection of the Objective 3 pipeline and identifies areas requiring review, validation, or future improvement.

## Confirmed Findings From Repository Inspection

### Confirmed Finding A: psm_dataset.py is not truly lazy-loaded
Status: Confirmed

Observation:
- The class description claims lazy parsing.
- During initialization, all matching spectra are appended into self.spectra.
- Large cohorts may therefore consume substantial memory.

Recommendation:
- Store scan offsets or identifiers rather than complete spectrum objects.
- Load spectra on demand in __getitem__.

Priority: Medium

---

### Confirmed Finding B: Expression linker correctly gates mock RNA generation
Status: Confirmed Good Practice

Observation:
- Mock TPM generation only occurs when --debug-expression is explicitly supplied.
- Deterministic seeding improves reproducibility.
- RNA provenance is written into manifest metadata.

Recommendation:
- Preserve current behavior.
- Prevent mock RNA evidence from appearing in publication-grade outputs.

Priority: Low

---

### Confirmed Finding C: Curated manifest rebuild improves reproducibility
Status: Confirmed Good Practice

Observation:
- Active cohort is explicitly reconstructed.
- Excluded runs are tracked.
- HLA provenance source is preserved.

Recommendation:
- Cross-check all patient-to-HLA mappings against publication supplementary material.

Priority: Medium

---

### Confirmed Finding D: Mass filter exists but production integration remains unverified
Status: Needs Validation

Observation:
- src/cnnlstm/mass_filter.py exists.
- Physical precursor-mass validation is implemented.
- Pipeline invocation path has not yet been confirmed.

Recommendation:
- Verify execution path from run_pipeline.sh.
- Measure candidate reduction before and after filtering.

Priority: High

---

### Confirmed Finding E: Unlabeled spectrum extraction uses efficient streaming
Status: Confirmed Good Practice

Observation:
- Generator-based processing avoids loading complete MGF files.
- Manifest filtering reduces accidental processing of inactive runs.

Recommendation:
- Add duplicate-spectrum reporting.
- Add scan-ID validation report.

Priority: Low

---

## High Priority Findings

### 1. Mass Filter Not Integrated Into Main Pipeline
Status: Review Required

The repository contains `src/cnnlstm/mass_filter.py`, which performs precursor mass validation of de novo predictions.

### 2. Silent Missing Spectrum Handling
Status: Review Required

Observed in spectral_dataset.py.

### 3. De Novo FDR Strategy Requires Validation
Status: Scientific Review Required

### 4. RNA Evidence Currently Debug Only
Status: Expected Limitation

## Overall Assessment

Current repository quality is substantially stronger than a prototype pipeline.

Most remaining risks are scientific-validation risks rather than software-engineering failures.

Highest-priority remaining audits:
1. cnnlstm_model.py
2. 05_train_denovo_model.py
3. 06_predict_denovo.py
4. 07_filter_neoantigens.py
5. 08_rank_candidates.py
