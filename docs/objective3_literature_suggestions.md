# Objective 3 Literature-Informed Suggestions

Date: 2026-06-09

Context: Objective 3 currently has a curated 31-run PXD005231 workflow, a custom CNN/LSTM-style de novo checkpoint at about 25.1% exact held-out peptide accuracy, and a Casanovo baseline run in progress on `data/mgf_unlabeled/*.mgf`.

## Proposal Scope Anchor

The current repo plan states the proposal objective as:

> Application of a Deep Learning System on de novo peptide sequencing (MS) data to identify neo-antigens via patient-specific CNN-LSTM training.

The Objective 3 presentation describes the same core idea: conventional database matching misses tumor-specific neoantigens not present in reference databases, so the project builds an autonomous CNN-LSTM de novo sequencing pipeline from MS/MS spectra, followed by neoantigen filtering and ranking.

Because the proposal is about three years old, the literature below should be used to modernize validation and framing, not to replace the promised objective. In practical terms:

- In scope: patient/cohort-specific CNN-LSTM training, held-out accuracy evaluation, de novo prediction on unlabeled spectra, FDR/filtering, HLA/expression ranking, and comparison to published/known neoantigens.
- In-scope as controls: Casanovo or InstaNovo as external benchmarks that tell us whether the custom CNN-LSTM is competitive enough and where it fails.
- Supporting validation only: Prosit/MS2PIP/DeepLC-style spectral or retention-time checks, model agreement tiers, run-level QC, and stricter FDR/calibration.
- Out of scope for the main claim: replacing the project with a transformer-only pipeline, claiming a full NeoDisc-like clinical proteogenomics system, or making biological expression claims from mock RNA.

## Papers Scanned

- Casanovo: "Sequence-to-sequence translation from mass spectra to peptides with a transformer model", Nature Communications 2024. https://www.nature.com/articles/s41467-024-49731-x
- InstaNovo: "InstaNovo enables diffusion-powered de novo peptide sequencing in large-scale proteomics experiments", Nature Machine Intelligence 2025. https://www.nature.com/articles/s42256-025-01019-5
- MHCquant2: "MHCquant2 refines immunopeptidomics tumor antigen discovery", Genome Biology 2025. https://link.springer.com/article/10.1186/s13059-025-03763-8
- NeoDisc: "A comprehensive proteogenomic pipeline for neoantigen discovery to advance personalized cancer immunotherapy", Nature Biotechnology 2024/2025. https://www.nature.com/articles/s41587-024-02420-y
- Non-canonical antigen proteogenomics: "Integrated proteogenomic deep sequencing and analytics accurately identify non-canonical peptides in tumor immunopeptidomes", Nature Communications 2020. https://www.nature.com/articles/s41467-020-14968-9
- MHCflurry 2.0: "Improved Pan-Allele Prediction of MHC Class I-Presented Peptides by Incorporating Antigen Processing", Cell Systems 2020. https://www.sciencedirect.com/science/article/pii/S2405471220302398
- Prosit: "Prosit: proteome-wide prediction of peptide tandem mass spectra by deep learning", Nature Methods 2019. https://www.nature.com/articles/s41592-019-0426-7
- MS2Rescore: "Data-Driven Rescoring Dramatically Boosts Immunopeptide Identification Rates", Molecular and Cellular Proteomics 2022. https://pmc.ncbi.nlm.nih.gov/articles/PMC9411678/

## Main Takeaways

1. The strongest newer de novo sequencing work is transformer-based and heavily pretrained.

   Casanovo treats MS/MS peaks as a sequence-to-sequence translation problem, and the 2024 paper reports training on tens of millions of labeled spectra plus a non-enzymatic fine-tuned variant. InstaNovo moves in the same direction, adding transformer decoding and diffusion-style refinement. For this proposal, these tools should be used as external baselines and sanity checks, while the main deliverable remains the patient-specific CNN-LSTM pipeline.

2. Immunopeptidomics needs non-tryptic, HLA-length-specific settings.

   HLA-I peptides are not ordinary tryptic peptides. For Objective 3, score interpretation should focus on 8-11/12 aa peptides, precursor mass consistency, charge, peak count, and HLA ligand plausibility. Casanovo's log shows many spectra being skipped for insufficient peaks, so run-level QC is important before biological interpretation.

3. Score calibration is as important as raw exact accuracy.

   A 25.1% exact match rate can still be useful if the high-confidence tail is much cleaner than the full prediction set. The next useful metric is not only global exact accuracy, but accuracy by score decile, edit-distance threshold, peptide length, charge, and run.

4. Downstream antigen discovery should be framed as proteogenomics-aware, not peptide prediction alone.

   NeoDisc and related immunopeptidomics work integrate genomics, transcriptomics, MS evidence, HLA binding, tumor specificity, and immunogenicity ranking. Objective 3 does not need to become NeoDisc, but its outputs should be described cautiously: de novo-derived candidate neoantigens, prioritized by HLA/expression evidence where available. Our current mock RNA is acceptable for debugging only; it should not be used as biological expression evidence.

5. False discovery risk grows fast in non-canonical and open search spaces.

   The non-canonical antigen literature emphasizes group-specific FDR, multiple search engines, retention-time/hydrophobicity checks, and synthetic peptide validation. For this repo, "found by de novo" should be a candidate-generating event, not final evidence.

6. Modern immunopeptidomics pipelines use predicted RT/MS2 features for rescoring.

   MHCquant2, MS2Rescore, Prosit, DeepLC, and MS2PIP show that fragment-intensity and retention-time prediction can raise sensitivity while keeping FDR controlled. Objective 3 can borrow this idea as a candidate validation layer even if the primary discovery remains de novo sequencing.

## Proposal-Aligned Interpretation

The central story should remain:

1. Curate MS/MS immunopeptidomics data and known PSMs.
2. Train a patient/cohort-specific CNN-LSTM de novo sequencing model.
3. Measure held-out peptide reconstruction accuracy.
4. Apply the trained model to unlabeled spectra after spectral subtraction.
5. Filter and rank candidate neoantigens using peptide quality, reference/proteome checks, HLA presentation, and expression metadata.
6. Compare with external baselines and literature/published neoantigen evidence.

The newer literature changes the standards around this story. It means the report should include stronger benchmarking, calibration, and limitations. It does not mean the proposal objective should be rewritten around transformers.

## Suggested Next Actions

### A. Finish the Casanovo baseline cleanly

When `results/casanovo_baseline.mztab` appears:

```bash
.venv/bin/python src/inference/convert_casanovo_output.py \
  --input results/casanovo_baseline.mztab \
  --output results/de_novo_candidates_casanovo.tsv \
  --manifest configs/sample_manifest.tsv
```

Then count converted rows, peptide length distribution, per-run rows, empty sample IDs, and Casanovo score distribution.

### B. Evaluate Casanovo on the same held-out test set

The current Casanovo run is on unlabeled dark spectra, so it cannot tell us exact peptide accuracy. A fair comparison needs Casanovo predictions on the known held-out spectra in `results/checkpoints_curated31_v2/test_set_psms.tsv`, then the same metrics as `10_evaluate_denovo_model.py`.

Recommended new output:

- `results/casanovo_accuracy_curated31_v2.json`
- metrics matching the custom model: exact accuracy, length accuracy, token accuracy, edit distance, edit <= 1 rate
- extra normalization flags: whether I/L are collapsed, whether modifications are stripped

### C. Build a high-confidence candidate tier by model agreement

After custom-model inference and Casanovo conversion exist, create candidate tiers:

- Tier 1: custom model and Casanovo agree exactly or edit distance <= 1, with mass error acceptable and HLA binding plausible.
- Tier 2: one model predicts a high-score peptide that passes mass, length, and HLA checks.
- Tier 3: low-score or biologically weak predictions kept only for debugging.

This reduces reliance on one weak global accuracy number.

### D. Add calibration plots before discovery claims

For the custom model and Casanovo held-out benchmark, generate:

- exact accuracy by score decile
- edit distance <= 1 by score decile
- accuracy by peptide length
- accuracy by run/sample
- fraction of predictions failing precursor mass consistency

The practical decision gate should be: "What score threshold gives an acceptable precision-like estimate?", not only "What is overall exact accuracy?"

### E. Make run-level spectrum QC explicit

Casanovo is skipping many spectra due to insufficient peaks. Add a QC table:

- input spectra per unlabeled MGF
- spectra accepted by Casanovo
- spectra skipped
- skipped fraction
- output predictions per run

Runs with extreme skipped fractions should be flagged before ranking neoantigens.

### F. Tighten immunopeptidomics-specific filters

Before ranking, require:

- HLA-I length focus: 8-11 aa, optionally keep 12 aa as exploratory
- precursor mass agreement
- no empty sample ID
- no exact match to known normal PSMs
- no exact mapping to healthy/reference-only sequences unless labeled as tumor-associated/non-mutated
- HLA binding/presentation evidence from MHCflurry/NetMHCpan

Keep mock RNA labels in output columns, but do not use mock TPM as biological evidence.

### G. Use predicted spectral/RT evidence as a second-pass validator

For the top candidate set, add one validation score from:

- Prosit or MS2PIP predicted fragment-intensity similarity
- DeepLC or equivalent retention-time plausibility
- mirror-plot-ready spectrum annotations for top hits

This mirrors the MHCquant2/MS2Rescore direction without replacing the current pipeline.

### H. Keep pretrained models as baselines, not the main deliverable

If the goal is to stay faithful to the proposal, the priority order should be:

1. Benchmark the custom CNN-LSTM on the curated held-out test set.
2. Benchmark Casanovo on the same held-out test set.
3. Optionally try InstaNovo as a second external baseline if dependencies are manageable.
4. Use pretrained results to contextualize the CNN-LSTM performance, not to erase the CNN-LSTM objective.
5. Only discuss fine-tuning a pretrained model as future work or as an extra experiment, unless the proposal/report explicitly allows extending the method.

## Suggested Decision Gate

Do not advance candidates to biological interpretation unless at least one of these is true:

- the model-specific held-out score threshold gives high-confidence exact/edit-distance performance;
- two independent de novo models agree on the peptide or near-neighbor sequence;
- the peptide passes MS plausibility, HLA presentation, and reference/tumor specificity checks;
- there is external validation evidence such as targeted MS, synthetic peptide comparison, or matched real RNA/genomics.

## Immediate Next Step After Casanovo Finishes

1. Convert `results/casanovo_baseline.mztab`.
2. Generate basic Casanovo candidate QC.
3. Run a same-test-set Casanovo benchmark.
4. Compare custom v2 vs Casanovo.
5. Use that comparison to decide whether the custom CNN-LSTM candidate set is strong enough to rank, or whether the report should emphasize model limitations and further training.
