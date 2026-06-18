# Detailed Code Audit

## Scope
This document is intended to become the permanent engineering and scientific audit record for the Objective 3 ImmunoDeNovo pipeline.

## Confirmed Audit Findings (Completed Reviews)

### 08_rank_candidates.py
- Implements biologically meaningful ranking using HLA binding predictions plus RNA-expression evidence.
- RNA evidence gating is stronger than many prototype pipelines because debug/mock RNA sources are explicitly prevented from supporting biological evidence classes.
- Duplicate candidate collapsing prevents peptide inflation in downstream reports.
- MHCflurry execution is performed sample-by-sample and may become slow for larger cohorts.
- Expression lookup construction repeatedly iterates over RNA tables and builds a Python dictionary rather than using indexed joins.
- Extensive use of DataFrame apply(axis=1) may become a bottleneck at scale.
- Sample matching relies primarily on sample_id == patient_id and falls back to run_id matching, which should be validated against manifest assumptions.
- Missing MHCflurry predictions silently downgrade candidates rather than producing a structured audit report.
- Evidence class assignment currently reduces all candidates into A/B/C classes and may not capture uncertainty from missing HLA typing or missing RNA.

Recommendations:
- Add ranking-stage audit metrics (binding failures, missing RNA, missing HLA, duplicate collapse counts).
- Replace row-wise expression assignment with vectorized joins where practical.
- Validate patient_id/sample_id/run_id mapping assumptions across all cohorts.
- Add provenance fields for MHCflurry version and predictor configuration.
- Consider explicit evidence subclasses for missing biological evidence rather than collapsing into Class C.
- Benchmark runtime on larger candidate collections.

Risk: Medium-High.

---

### Existing audited files
(See previous sections retained in repository history.)