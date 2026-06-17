# Final Manual Rerun Commands

This runbook assumes the curated 31-run `PXD005231` manifest is active and the 9 unverified TIL2/TIL4/DC runs are excluded in `configs/excluded_runs.tsv`.

## 1. Rebuild Metadata and Validate Inputs

```bash
python3 src/data_prep/00f_rebuild_curated_manifest.py \
    --base configs/sample_manifest.tsv.bak \
    --psm-dir data/psms \
    --output configs/sample_manifest.tsv \
    --excluded-output configs/excluded_runs.tsv

python3 src/data_prep/00e_filter_human_fasta.py \
    --input data/reference/uniprot_human_reviewed.fasta \
    --output data/reference/uniprot_human_reviewed.only_human.fasta

python3 src/data_prep/00d_link_expression.py \
    --manifest configs/sample_manifest.tsv \
    --fasta data/reference/uniprot_human_reviewed.only_human.fasta \
    --output data/expression_matrix.tsv

python3 src/validation/preflight_validate.py \
    --reference-fasta data/reference/uniprot_human_reviewed.only_human.fasta \
    --expected-active-runs 31 \
    --output results/preflight_before_rerun.md

python3 src/validation/provenance_audit.py \
    --output-md results/provenance_audit_before_rerun.md \
    --output-json results/provenance_audit_before_rerun.json
```

`mock_debug` RNA is allowed for pipeline debugging only, and only when `--debug-expression` is used. Do not present TPM-supported Class A/B calls as biological evidence until patient-matched real RNA is supplied.

## 2. Baseline and Spectral Subtraction

```bash
python3 src/data_prep/04_extract_psms.py \
    --input-dir data/psms \
    --output-file results/immunopeptidome_psms.tsv \
    --manifest configs/sample_manifest.tsv

python3 src/data_prep/05_extract_unlabeled_spectra.py \
    --mgf-dir data/mgf \
    --psms results/immunopeptidome_psms.tsv \
    --output-dir data/mgf_unlabeled \
    --manifest configs/sample_manifest.tsv \
    --clean-output
```

## 3. Model Accuracy Check

Current baseline full held-out accuracy from the old checkpoint was:

```text
Exact peptide accuracy:        10.89%
Length accuracy:               78.52%
Token accuracy excluding PAD:  57.67%
Position amino-acid accuracy:  49.89%
Mean edit distance:            4.705
Edit distance <= 1 rate:       20.25%
```

Because exact accuracy is weak, retraining is recommended after the curated manifest cleanup.

Quick smoke check:

```bash
.venv/bin/python -u src/evaluation/10_evaluate_denovo_model.py \
    --checkpoint results/checkpoints/neoepitope_production_best.pth \
    --test-psms results/checkpoints/test_set_psms.tsv \
    --mgf-dir data/mgf \
    --output results/model_accuracy_smoke.json \
    --max-runs 1 \
    --device cpu
```

Full accuracy check:

```bash
.venv/bin/python -u src/evaluation/10_evaluate_denovo_model.py \
    --checkpoint results/checkpoints/neoepitope_production_best.pth \
    --test-psms results/checkpoints/test_set_psms.tsv \
    --mgf-dir data/mgf \
    --output results/model_accuracy_report.json \
    --device cpu
```

## 4. Retrain If You Want To Improve Accuracy

Recommended retraining command:

```bash
.venv/bin/python src/training/05_train_denovo_model.py \
    --mgf-dir data/mgf \
    --checkpoint-dir results/checkpoints_curated31 \
    --epochs 80 \
    --batch-size 16 \
    --lr 0.0005 \
    --num-workers 0
```

After training:

```bash
.venv/bin/python -u src/evaluation/10_evaluate_denovo_model.py \
    --checkpoint results/checkpoints_curated31/neoepitope_production_best.pth \
    --test-psms results/checkpoints_curated31/test_set_psms.tsv \
    --mgf-dir data/mgf \
    --output results/model_accuracy_curated31.json \
    --device cpu
```

Proceed with de novo discovery only if exact peptide accuracy and edit-distance metrics improve meaningfully over the old baseline.

## 5. De Novo Discovery and Ranking

Use the old checkpoint only for debugging, or the curated retrained checkpoint if retraining improved metrics.

```bash
python3 src/inference/06_predict_denovo.py \
    --model results/checkpoints_curated31/neoepitope_production_best.pth \
    --mgf_dir data/mgf_unlabeled \
    --output results/de_novo_candidates.tsv \
    --device cpu

python3 src/postprocess/07_filter_neoantigens.py \
    --input results/de_novo_candidates.tsv \
    --psms results/immunopeptidome_psms.tsv \
    --fasta data/reference/uniprot_human_reviewed.only_human.fasta \
    --manifest configs/sample_manifest.tsv \
    --output results/filtered_neoantigens.tsv \
    --score_cutoff -0.5 \
    --min_psm_support 2

python3 src/postprocess/08_rank_candidates.py \
    --input results/filtered_neoantigens.tsv \
    --manifest configs/sample_manifest.tsv \
    --output results/ranked_neoantigens.tsv \
    --binding_rank_cutoff 2.0 \
    --tpm_cutoff 1.0
```

## 6. Final Checks and Evaluation

```bash
python3 src/validation/preflight_validate.py \
    --reference-fasta data/reference/uniprot_human_reviewed.only_human.fasta \
    --expected-active-runs 31 \
    --output results/preflight_after_rerun.md

python3 src/validation/provenance_audit.py \
    --output-md results/provenance_audit_after_rerun.md \
    --output-json results/provenance_audit_after_rerun.json

python3 src/evaluation/09_evaluate_neoantigens.py \
    --input results/ranked_neoantigens.tsv \
    --validated data/reference/s2_dataset_extracted/Dataset1/Dataset1.txt \
    --manifest configs/sample_manifest.tsv \
    --output results/evaluation_report.md \
    --top_n 10 25 50
```

## 7. Optional Quarantine Cleanup

Dry-run:

```bash
python3 src/maintenance/quarantine_workspace.py
```

Apply:

```bash
python3 src/maintenance/quarantine_workspace.py --apply
```

This moves clutter into `archive/quarantine/` and writes `archive/quarantine/MANIFEST.tsv`. It does not delete files.
