# Linked RNA Cohort Plan

Goal: find at least one cohort where the same patient has MS/MS raw data, sample-to-patient mapping, HLA alleles, and patient-matched RNA-seq or expression data.

## Acceptance Criteria

A cohort is usable only if the manifest can carry, for each active run:

- `study_id`
- `run_id`
- `patient_id`
- `validation_id`
- `hla_alleles`
- `hla_source`
- `rna_expr_path`
- `rna_source = real_patient_matched`
- `raw_source`
- `psm_source`
- `include_in_pipeline = True`

## Search Strategy

1. Start from the known MS cohort in PRIDE, then look for patient-matched RNA in the publication supplements, GEO/ArrayExpress links, or associated project metadata.
2. Prefer studies with explicit patient identifiers in both the MS and RNA tables.
3. Reject cohorts where RNA is only a cell line surrogate, pooled control, public reference tissue, or simulated debug profile.
4. Confirm HLA alleles either from the publication, matched clinical typing, or a clearly cited metadata source.
5. Keep a note of any partial cohort where only some patients have RNA; those rows can remain in the manifest, but only the matched rows may support expression evidence.

## Practical Evidence Checks

- RAW/PSM run IDs must map cleanly to the manifest.
- RNA paths must exist locally or be reproducibly downloadable.
- `rna_source` must be a patient-matched label before ranking can use TPM as biological evidence.
- Preflight must pass with no unknown samples and no stale outputs.
- Provenance audit must document where each RNA file came from.

## Output

Produce one short cohort note per candidate study:

- accession or study name
- where the MS data came from
- where the RNA came from
- how patient IDs were matched
- what HLA evidence exists
- whether the cohort is fully usable or only partially usable
