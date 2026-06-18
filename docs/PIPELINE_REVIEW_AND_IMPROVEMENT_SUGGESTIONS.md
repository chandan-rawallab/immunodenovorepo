# Pipeline Review and Improvement Suggestions

## Newly Confirmed Findings

### Critical: Verify train/test peptide leakage
Observation:
- Current training workflow appears to split spectra rather than unique peptide sequences.
- Repeated peptide observations may therefore occur in both train and test sets.

Recommendation:
- Perform grouped splitting by peptide sequence.
- Audit overlap between training peptides and test_set_psms.tsv.

Priority: P1

---

### High: Verify mass-filter integration
Observation:
- mass_filter.py is implemented and scientifically valuable.
- Integration into the production inference path has not yet been confirmed.

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
- Benchmark runtime on larger cohorts.

Priority: P1-P2

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
1. src/postprocess/08_rank_candidates.py
2. src/evaluation/10_evaluate_denovo_model.py
3. src/validation/preflight_validate.py
4. src/validation/provenance_audit.py
5. run_pipeline.sh