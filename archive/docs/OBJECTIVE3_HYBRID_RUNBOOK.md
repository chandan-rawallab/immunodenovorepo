# Objective 3 Hybrid Runbook

This workspace now uses `experiments` as the primary Objective 3 implementation:

- MaxQuant/Swiss-Prot immunopeptidome construction remains the upstream label source.
- The custom CNN/LSTM model remains the intended de novo sequencing engine.
- The downstream contract, subtraction, ranking, and reporting logic follows the stronger shape from the sibling `../obj3` workspace.

## Current Status

Check MGF/MaxQuant coverage:

```bash
PYTHONPATH=src python3 -m objective3.cli status \
  --mgf-dir /home/amity/hla_data_mgf \
  --results-dir production_results
```

The latest checked state was:

- `20` MGF files in `/home/amity/hla_data_mgf`
- `11` completed MaxQuant `msms_*.txt` files in `production_results`
- `11` matched MGF/MSMS sample names
- `4,331` raw PSM rows

## Normalize MaxQuant PSMs

```bash
PYTHONPATH=src python3 -m objective3.cli normalize-maxquant-psms \
  --results-dir production_results \
  --output results/objective3/immunopeptidome_psms.tsv
```

This writes the shared Objective 3 PSM contract:

```text
sample_id  spectrum_id  peptide  q_value  source_file
```

MaxQuant `PEP` is used as the normalized `q_value` field when a direct q-value column is not present.

## Train The Custom Model

Use a Python environment with working PyTorch:

```bash
PYTHONPATH=src python3 src/train_production.py \
  --mgf-dir /home/amity/hla_data_mgf \
  --results-dir production_results \
  --checkpoint-dir results/checkpoints \
  --epochs 100 \
  --batch-size 32
```

The default system `python3` currently lacks `torch`. The `objective3-casanovo` conda environment has PyTorch installed, but tensor execution currently fails with an MKL symbol error on this machine, so that environment needs repair before full training or prediction.

## Predict De Novo Candidates

After a working PyTorch environment is available:

```bash
PYTHONPATH=src python3 -m objective3.cli predict-denovo \
  --mgf-dir /home/amity/hla_data_mgf \
  --checkpoint results/checkpoints/neoepitope_production_epoch_100.pth \
  --output results/objective3/de_novo_predictions.tsv
```

For a small smoke run:

```bash
PYTHONPATH=src python3 -m objective3.cli predict-denovo \
  --mgf-dir /home/amity/hla_data_mgf \
  --checkpoint neoepitope_model_epoch_50.pth \
  --output results/objective3/de_novo_predictions_smoke.tsv \
  --max-spectra-per-file 1
```

## Postprocess And Rank

```bash
PYTHONPATH=src python3 -m objective3.cli postprocess \
  --psms results/objective3/immunopeptidome_psms.tsv \
  --denovo results/objective3/de_novo_predictions.tsv \
  --outdir results/objective3/final
```

Optional evidence files:

```bash
PYTHONPATH=src python3 -m objective3.cli postprocess \
  --psms results/objective3/immunopeptidome_psms.tsv \
  --denovo results/objective3/de_novo_predictions.tsv \
  --sample-metadata path/to/sample_metadata.tsv \
  --binding-input path/to/binding_predictions.tsv \
  --outdir results/objective3/final
```

The final outputs are:

- `de_novo_candidates.tsv`
- `ranked_neoantigens.tsv`
- `run_report.md`
- `run_report.json`

## Smoke Outputs Already Generated

The current repository contains a non-model smoke de novo input at:

```text
results/objective3/de_novo_predictions_smoke.tsv
```

and verified downstream outputs at:

```text
results/objective3/final_smoke/
```

This confirms that database-supported peptides are subtracted and de novo-only candidates are ranked and reported.
