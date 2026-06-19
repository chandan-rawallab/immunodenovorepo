# Pipeline Review and Improvement Suggestions

## Newly Confirmed Findings

### Critical: Silent spectrum lookup failures in spectral_dataset.py
Observation:
- Missing spectra are silently converted into zero-valued training examples.
- Broad exception handling hides indexing and parsing problems.

Recommendation:
- Add explicit failure accounting.
- Log missing run_id/scan_id combinations.
- Fail training when failure rate exceeds threshold.

Priority: P1

---

### High: Verify train/test peptide leakage
Observation:
- Current training workflow appears to split spectra rather than unique peptide sequences.
- Repeated peptide observations may therefore occur in both train and test sets.

Recommendation:
- Perform grouped splitting by peptide sequence.
- Audit overlap between training peptides and test_set_psms.tsv.

Priority: P1

---

### High: Improve spectrum preprocessing
Observation:
- No intensity transformation is currently applied.
- No top-N peak filtering is currently applied.

Recommendation:
- Evaluate log1p or sqrt intensity scaling.
- Benchmark top-200 and top-300 peak filtering.

Priority: P1-P2

---

### High: Verify mass-filter integration
Observation:
- mass_filter.py is implemented and scientifically valuable.
- Integration into the production inference path has not yet been fully verified.

Recommendation:
- Trace execution from run_pipeline.sh through inference outputs.
- Measure candidate reduction statistics.

Priority: P1

---

### High: Optimize neoantigen filtering performance
Observation:
- 07_filter_neoantigens.py reparses FASTA data during source-protein discovery.
- Multiple row-wise pandas operations may become bottlenecks at scale.

Recommendation:
- Cache protein mappings.
- Reduce repeated FASTA scans.

Priority: P1-P2

---

### High: Strengthen ranking-stage provenance and auditing
Observation:
- 08_rank_candidates.py performs critical final evidence classification.
- Missing HLA predictions, missing RNA data, and duplicate-candidate collapsing are not comprehensively summarized.

Recommendation:
- Emit ranking statistics and evidence-stage QC reports.
- Record MHCflurry version and prediction provenance.

Priority: P1

---

### Medium: Improve checkpoint provenance
Recommendation:
- Save training configuration.
- Save manifest hash/version.
- Save git commit identifier.

Priority: P2

---

### Medium: psm_dataset memory scaling
Observation:
- Dataset currently stores matched spectra in memory.

Recommendation:
- Implement indexed retrieval or true lazy loading.

Priority: P2

---

### Medium: cnnlstm_model architecture improvements
Recommendation:
- Evaluate attention mechanism.
- Evaluate mass-aware decoding.

Priority: P2

---

## Next Audit Targets
1. src/evaluation/10_evaluate_denovo_model.py
2. src/validation/preflight_validate.py
3. src/validation/provenance_audit.py
4. run_pipeline.sh