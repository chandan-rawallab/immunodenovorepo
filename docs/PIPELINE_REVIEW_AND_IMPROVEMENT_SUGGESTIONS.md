# Pipeline Review and Improvement Suggestions

## Newly Confirmed Findings

### High: Strengthen ranking-stage provenance and auditing
Observation:
- 08_rank_candidates.py performs critical final evidence classification.
- Missing HLA predictions, missing RNA data, and duplicate-candidate collapsing are not comprehensively summarized in an audit output.

Recommendation:
- Emit ranking statistics and evidence-stage QC reports.
- Record MHCflurry version and prediction provenance.
- Record duplicate-collapse counts and evidence transitions.

Priority: P1

---

### High: Optimize expression assignment
Observation:
- Expression scoring relies on row-wise DataFrame apply operations.
- Runtime may grow substantially with larger candidate sets.

Recommendation:
- Replace repeated row-wise evaluation with indexed merges/vectorized lookups.

Priority: P1-P2

---

### Medium: Validate sample mapping assumptions
Observation:
- Ranking logic assumes sample_id maps directly to patient_id and falls back to run_id matching.

Recommendation:
- Explicitly validate manifest relationships across future cohorts.

Priority: P2

---

## Next Audit Targets
1. src/evaluation/10_evaluate_denovo_model.py
2. src/validation/preflight_validate.py
3. src/validation/provenance_audit.py
4. run_pipeline.sh