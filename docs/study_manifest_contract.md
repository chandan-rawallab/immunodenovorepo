# Study Manifest Contract

Objective 3 should treat each dataset as a study adapter that emits the same manifest contract. The pipeline may ingest PRIDE, local RAW/MGF/PSM folders, curated publication tables, or future patient cohorts, but downstream steps must consume the same fields.

## Required Columns

| Column | Purpose |
|:--|:--|
| `study_id` | Stable study or cohort identifier, such as `PXD005231`. |
| `run_id` | Raw/MS run identifier used to match RAW, MGF, PSM, and prediction files. |
| `patient_id` | Biological sample or patient identifier used for ranking. |
| `validation_id` | Identifier used for external validation tables when different from `patient_id`. |
| `filename` | Source filename for the run. |
| `cohort` | Cohort/accession label for grouping. |
| `sample_role` | Role such as `hla_peptidome`, `apheresis_control`, or `pooled_til`. |
| `hla_alleles` | Comma-separated HLA alleles in standard form. |
| `hla_source` | Source of HLA evidence, for example `curated_publication`, `clinical_ngs`, `sdrf`, `manual`, `inferred`, or `missing`. |
| `rna_expr_path` | Path to patient expression TSV, if available. |
| `rna_source` | RNA evidence label. See allowed labels below. |
| `raw_source` | Provenance label for raw spectra, such as `pride_raw` or `local_raw`. |
| `psm_source` | Provenance label for database-search PSMs, such as `local_maxquant_psm`. |
| `include_in_pipeline` | Boolean-like flag controlling whether the row is active. |
| `notes` | Human-readable provenance or curation notes. |

## RNA Evidence Labels

Only patient-matched labels may support Class A/B expression evidence in ranking:

- `real_patient_matched`
- `patient_matched`
- `provided_patient_matched`
- `clinical_rna`
- `matched_rna`

These labels are allowed only when explicit debug or override mode is used, and they cannot support biological expression evidence:

- `mock_debug`
- `surrogate`
- `inferred`
- `provided_override`
- `missing`

## Adapter Rule

A new study-specific adapter should stop after producing this manifest plus source files. It should not change filtering, ranking, or evaluation code. The shared pipeline is responsible for validating provenance, rejecting stale or incomplete state, and limiting claims based on available evidence.


The expression linker must fail closed by default: if no patient-matched RNA exists, leave `rna_expr_path` empty and `rna_source = missing` unless debug mode is explicitly enabled.
