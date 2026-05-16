# Objective 3 — Master Plan (Revised)
## Aligned to: `archive/workplan/obj3.txt`
## Primary workspace: `/home/amity/Documents/experiments/`

> **Objective:** Application of a Deep Learning System on de novo peptide sequencing (MS) data to identify neo-antigens via patient-specific CNN-LSTM training.

---

## Notes on Key Decisions

**On PSM extraction:** The MaxQuant `msms.txt` files in `production_results/` already exist — a converter is the fastest path since re-running MaxQuant on 61 RAW files would take days of compute. The converter treats MaxQuant as the database-search tool (which it is — PEAKS X equivalent). Future fresh dataset runs can use MaxQuant directly from scratch; the converter is only needed for the current pre-existing results.

**On workspace:** All work is done inside `experiments/`. The `obj3/` workspace is treated as a reference implementation for the ranking logic only — its Python package (`obj3/objective3/`) is borrowed, not its directory structure.

---

## PHASE 0 — Academic Directory Restructuring
*Do this first — clean structure prevents confusion across all subsequent phases*

### Current Problems

The `experiments/` directory has accumulated organic clutter:
- Training scripts `train_lite.py`, `train_v2.py`, `train_production.py` scattered at `src/` root with no version logic
- Results split across `results/`, `results/checkpoints/`, `results/objective3/`, `results/identifications/` with no date or run metadata
- Logs buried in `logs/` with inconsistent naming (`training_output_v2.log`, `train_v2_medium_lite.log`, `train_v2.pid`)
- Temporary working files mixed with actual outputs (`temp/lite_mgf/`, `temp/medium_lite_mgf/`)
- Two separate `objective3` packages: `src/objective3/` (experiments) and `obj3/objective3/` (reference)
- `production_results/` contains MaxQuant msms files named with raw file names — not linked to any patient/sample metadata
- `OBJECTIVE3_MANUAL.md` sitting at the project root with no versioning or section linking

### Naming Convention

All files are named after **the scientific process they execute**, not the implementation detail. Three rules:
1. **No dev-ops jargon** — `production`, `lite`, `scale` describe infrastructure, not biology
2. **No version numbers in filenames** — `v2`, `_final` belong in git commit messages
3. **One vocabulary** — the pipeline step name, matching the 4 proposal activities

| Old name (current) | New name | Reason |
|:---|:---|:---|
| `dataretrieval.py` | `download_pride.py` | Names the action and source |
| `conversion_pipeline.py` | `convert_raw_to_mgf.py` | States input → output |
| `scale_search.py` | `search_database.py` | Names the scientific step |
| `msms_to_psm_tsv.py` *(new)* | `extract_psms.py` | Concise, names the product |
| `dl_module/` *(folder)* | `cnnlstm/` | Names the algorithm, not the abstraction |
| `model.py` | `cnnlstm_model.py` | Unambiguous across the project |
| `production_dataset.py` | `spectral_dataset.py` | Names the data type, not deployment tier |
| `dataset.py` | `psm_dataset.py` | Distinguishes from spectral dataset |
| `train_v2.py` / `train_lite.py` / `train_production.py` | `train_denovo_model.py` | One canonical entry point |
| `train_patient_specific.py` *(new)* | `train_denovo_model.py` | Same — merged into one |
| `predict_lite.py` | `predict_denovo.py` | Names the scientific output |
| `proteome_filter.py` | `filter_neoantigens.py` | Names the biological goal |
| `postprocess.py` | `rank_candidates.py` | Names what it actually does |
| `evaluate.py` *(new)* | `evaluate_neoantigens.py` | Symmetric with filter step |
| `mqpar_production.xml` / `mqpar.xml` | `maxquant_config.xml` | Tool + purpose, no tier jargon |
| `casanovo_cpu.yaml` | `casanovo_config.yaml` | Hardware detail belongs inside the file |
| `sample_metadata.tsv` | `sample_manifest.tsv` | Standard bioinformatics term |
| `OBJECTIVE3_MANUAL.md` | `docs/pipeline_manual.md` | Descriptive, not named after a project code |

### Target Academic Structure

```
experiments/
├── README.md                            ← Project overview, dataset provenance, how to run
├── CHANGELOG.md                         ← Timestamped record of all major actions
│
├── data/
│   ├── raw/                             ← Symlink → /home/amity/hla_data_raw/
│   ├── mgf/                             ← Symlink → /home/amity/hla_data_mgf/
│   ├── psms/
│   │   ├── sample_manifest.tsv          ← patient_id, hla_alleles, cohort, run_ids
│   │   └── psms_<patient_id>.tsv        ← Extracted PSMs per patient (from extract_psms.py)
│   └── reference/
│       ├── uniprot_human_reviewed.fasta
│       └── reference_manifest.json      ← FASTA checksum + version
│
├── src/
│   ├── 01_download_pride.py             ← Activity 1: fetch PXD005231 from PRIDE
│   ├── 02_convert_raw_to_mgf.py         ← Activity 1: RAW → MGF via MSConvert
│   ├── 03_search_database.py            ← Activity 1: MaxQuant parallel launcher
│   ├── 04_extract_psms.py               ← Activity 1: msms.txt → immunopeptidome_psms.tsv
│   ├── 05_train_denovo_model.py         ← Activity 2: CNN-LSTM trained on patient PSMs
│   ├── 06_predict_denovo.py             ← Activity 2/3: apply model to unlabelled spectra
│   ├── 07_filter_neoantigens.py         ← Activity 3: FDR, length, missense mutation filters
│   ├── 08_rank_candidates.py            ← Activity 3: binding affinity + expression ranking
│   ├── 09_evaluate_neoantigens.py       ← Activity 4: compare vs. published neoantigens
│   │
│   ├── cnnlstm/                         ← CNN-LSTM model package
│   │   ├── __init__.py
│   │   ├── cnnlstm_model.py             ← NeoepitopeSeq2Seq architecture
│   │   ├── spectral_dataset.py          ← Full MGF loader with binning
│   │   └── psm_dataset.py              ← PSM-paired training set loader
│   │
│   └── _archive/                        ← Old scripts, never deleted, never run
│       ├── train_lite.py
│       ├── train_v2.py
│       ├── train_production.py
│       └── predict_lite.py
│
├── configs/
│   ├── maxquant_config.xml              ← MaxQuant search parameters
│   ├── casanovo_config.yaml             ← Casanovo inference settings
│   └── sample_manifest.tsv             ← Patient metadata (HLA, RNA-seq paths)
│
├── results/
│   ├── YYYYMMDD_<descriptive_name>/     ← One dated directory per run
│   │   ├── run_config.json
│   │   ├── immunopeptidome_psms.tsv
│   │   ├── de_novo_candidates.tsv
│   │   ├── ranked_neoantigens.tsv
│   │   ├── run_report.md
│   │   └── model_checkpoint.pth
│   └── _archive/                        ← Existing ad-hoc results kept for reference
│
├── logs/
│   └── YYYYMMDD_HHMMSS_<step>.log       ← One log per step execution
│
├── tests/
│   ├── test_extract_psms.py
│   ├── test_cnnlstm_model.py
│   └── test_rank_candidates.py
│
├── docs/
│   ├── pipeline_manual.md               ← Full research manual (moved from root)
│   └── runbook.md                       ← Operator quick-start
│
└── archive/
    └── workplan/
        └── obj3.txt
```

> **Note on step numbering (01–09):** Prefixing scripts with step numbers is standard in published bioinformatics workflows (nf-core, ENCODE pipelines). It makes the pipeline order self-documenting — anyone reading `src/` immediately understands the sequence without reading docs.

### Restructuring Actions (in order)

**R1. Rename `src/dl_module/` → `src/dl_model/`** (more standard bioinformatics naming)

**R2. Create `src/data_prep/`, `src/db_search/`, `src/training/`, `src/inference/`, `src/postprocess/`**
- Move scripts to their correct home
- Archive old training scripts under `src/training/_archive/`

**R3. Create `data/psms/` with `sample_manifest.tsv`**
- Map each `production_results/msms_*.raw.txt` → patient ID using the Bassani-Sternberg 2016 sample table
- Rename files: `msms_TIL3_run1.tsv`, `msms_TIL3_run2.tsv` etc. (human-readable)

**R4. Create `results/YYYYMMDD_<run_name>/` convention**
- Move current `results/checkpoints/`, `results/de_novo_hits*` → `results/_archive/`
- All future runs get their own dated directory

**R5. Move `OBJECTIVE3_MANUAL.md` → `docs/OBJECTIVE3_MANUAL.md`**
- Root directory should only have `README.md` and `CHANGELOG.md`

**R6. Create `data/raw/` and `data/mgf/` as symlinks**
```bash
ln -s /home/amity/hla_data_raw /home/amity/Documents/experiments/data/raw
ln -s /home/amity/hla_data_mgf /home/amity/Documents/experiments/data/mgf
```

**R7. Create `configs/sample_metadata.tsv`** — the patient manifest (HLA alleles, RNA-seq paths)

---

## PHASE A — Fix Existing Pipeline Infrastructure (Est. 1–2 days)
*All work in `experiments/`*

**A1. Write `src/data_prep/msms_to_psm_tsv.py`**

Converts existing MaxQuant `production_results/` output into the pipeline contract format. This is fast (no re-compute) since MaxQuant already ran.

> **Why converter, not re-running MaxQuant?** The 61 RAW files have already been searched. Re-running MaxQuant takes 24–48h per batch on this hardware. The converter extracts what we need in minutes. For future datasets, MaxQuant can be run directly and the output fed in the same format.

Input format: MaxQuant `msms.txt` columns — `Sequence`, `PEP`, `Scan number`, `Raw file`, `Decoy`
Output contract:
```
sample_id   spectrum_id           peptide    q_value   source_file
TIL3_R1     scan:24382            AASAAAAEL  0.020344  msms_TIL3_R1.tsv
```

Filter rules:
- `Decoy != "+"` — remove decoy hits
- `PEP ≤ 0.01` — use PEP as proxy for q-value at 1% FDR
- Peptide length 8–11 AA (HLA-I canonical)
- Remove rows with empty `Sequence`

**A2. Fix HLA allele format in postprocess (1 hour)**

In `src/postprocess/postprocess.py` (ported from `obj3/`):
- `_normalize_hla()` currently strips all formatting, yielding `A0201`
- MHCflurry expects `HLA-A*02:01`
- Fix: detect if the input already has `*` → pass through unchanged; otherwise reconstruct

```python
# CURRENT (broken)
value = value.replace("HLA-", "").replace("*", "").replace(":", "")
# returns "A0201" → MHCflurry rejects

# FIX
# Keep original format; MHCflurry accepts "HLA-A*02:01" and "A*02:01"
# Only strip if the format is already 4-digit like "A0201" (no asterisk)
if "*" not in value:
    # convert A0201 → HLA-A*02:01
    value = f"HLA-{value[0]}*{value[1:3]}:{value[3:5]}"
```

**A3. Add HLA-I length filter (30 min)**
- After PSM subtraction, retain only peptides 8–11 AA
- `CANONICAL_HLA_I_PATTERN = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]{8,11}$")`
- Rejects: any de novo prediction outside this range → `rejected_candidates.tsv`

**A4. Add single missense mutation detection in `src/postprocess/proteome_filter.py` (4 hours)**

The proposal explicitly requires filtering on *"presence of single missense mutation"*. Current proteome filter does binary in/out matching. Upgrade:
- For each de novo candidate, find the nearest match in `uniprot_human_reviewed.fasta`
- Use Levenshtein distance = 1 (exactly one amino acid substitution)
- If match found at distance 1: annotate `missense_position`, `wt_aa`, `mut_aa`, `source_protein`
- If match found at distance 0: it's a self-peptide → subtract
- If no match within distance 1: Class C (may still be a neoepitope from frameshift/indel)

---

## PHASE B — Patient-Specific CNN-LSTM Training (Est. 3–5 days)
*Core proposal requirement: train on the patient's own immunopeptidome PSMs*

**B1. Build the PSM-paired training dataset**
- Input: `immunopeptidome_psms.tsv` (from A1) + `data/mgf/` files
- For each PSM row: look up the spectrum by scan number in the corresponding MGF file
- Output: paired list of `(spectrum_vector, peptide_sequence)` — the labelled training set
- Expected scale: ~15,000 pairs from 40 MGF files

**B2. Proper 70 / 15 / 15 split**
- Stratify by peptide length (9-mers, 10-mers, etc.) to avoid length bias
- Split is per-spectrum, NOT per-peptide (same peptide can appear in multiple spectra)
- Test set is **locked** — never used until final evaluation

**B3. Train `src/training/train_patient_specific.py`**
- Load all ~15,000 PSM-paired spectra
- Model: `dl_model/model.py` (existing CNN-LSTM) — no architecture change needed
- Training: 30 epochs, early stopping on validation loss (patience = 5)
- Save: `results/YYYYMMDD_patient_specific_v1/model_checkpoint.pth`
- Log: per-epoch train loss, val loss, per-amino-acid accuracy on val set

**B4. Optional: Fine-tune Casanovo alongside (parallel benchmark)**
- `casanovo train` with `--train_from_scratch False` on the same 15,000 PSM pairs
- This gives a direct comparison: proposal method (custom CNN-LSTM trained from scratch) vs. SOTA (fine-tuned Casanovo)
- Both are valid implementations of the proposal since Casanovo is authored by Tran et al. (same group cited in obj3.txt)

---

## PHASE C — De Novo Sequencing on All Spectra (Est. 2–3 days)
*Apply the trained patient-specific model to every spectrum*

**C1. Run model on labelled spectra (accuracy validation)**
- Apply trained model to the PSM-matched spectra
- Compare prediction vs. ground-truth peptide
- Metric: per-amino-acid accuracy (target: >85% per Tran et al. 2019 benchmark)
- This validates the model before applying to unknowns

**C2. Run model on unlabelled spectra (neoepitope candidates)**
- Unlabelled = spectra in MGF with no PSM match from database search
- Apply patient-specific model → sequence predictions + confidence scores
- Save: `de_novo_candidates.tsv`

**C3. FDR estimation via target-decoy approach**
- Generate decoy predictions by reverse-complementing the predicted sequences
- At each score threshold, compute: `FDR = (# decoy hits) / (# target hits)`
- Apply 5% FDR cutoff — retain only `q_value ≤ 0.05` predictions
- This is the "quality control procedure" the proposal mandates

---

## PHASE D — Shortlisting & Ranking (Est. 2–3 days)
*Apply proposal's exact filter criteria and ranking criteria*

**D1. Apply proposal-mandated filters (in order)**
1. FDR ≤ 0.05 (from Phase C3)
2. Length 8–11 AA (from Phase A3)
3. **Single missense mutation** — Levenshtein-1 from reference proteome (from Phase A4)
4. **Mutations at flanking positions** — check if missense is at positions 1 or last of the peptide
5. PSM support count ≥ 2 (seen in at least 2 spectra)

**D2. Rank by binding affinity (MHCflurry)**
- Fix A2 ensures allele format is correct
- `binding_rank ≤ 2.0` → binding-supported

**D3. Link RNA-seq expression data**
- Source: Bassani-Sternberg 2016 companion RNA-seq (GEO accession with the paper)
- Format: TSV with `gene`, `expression_tpm` per patient
- Update `configs/sample_metadata.tsv` → `rna_expr_path` column per sample
- `expression_tpm ≥ 1.0` → expression-supported

**D4. Assign Evidence Classes and produce final output**
- Class A: missense match + binding + expression ← what the paper needs
- Class B: binding + expression (no direct mutation)
- Class C: MS-only

---

## PHASE E — Evaluation with RNA-seq (Est. 2–3 days)
*Proposal Activity 4: "Final evaluation of top neo antigens shall be performed using the available RNA-seq dataset"*

**E1. Obtain Bassani-Sternberg 2016 validated neoantigen list**
- Supplementary Table from the 2016 paper — validated neoantigens per patient
- Cross-reference against our top-ranked candidates
- Compute precision, recall at top-10, top-25, top-50 cutoffs

**E2. RNA-seq expression validation for top candidates**
- For each Class A/B candidate: confirm TPM ≥ 1 in the patient's RNA-seq
- Flag any candidate where gene is unexpressed (TPM < 0.5) → likely false positive

**E3. Produce paper-ready output**
- Update `docs/OBJECTIVE3_MANUAL.md` with final results
- Create results table: top-20 candidates per patient with evidence class, scores, mutation annotation
- Figures: evidence class pie chart, Casanovo score distribution, precision-recall curve

---

## Prioritised Immediate Actions

| # | Action | Phase | Time | Blocks |
|:--|:---|:---|:---|:---|
| 1 | **Directory restructure** — create academic layout, move files | 0 | 2h | Everything |
| 2 | **Write `msms_to_psm_tsv.py`** — convert MaxQuant PSMs in `production_results/` | A1 | 3h | Training data + subtraction |
| 3 | **Fix HLA format bug** in postprocess.py | A2 | 1h | Binding ranking |
| 4 | **Add missense mutation detection** to `proteome_filter.py` | A4 | 4h | Proposal filter requirement |
| 5 | **Build PSM-paired training dataset** | B1 | 4h | Patient-specific training |
| 6 | **Train patient-specific model** | B2–B3 | 2 days | All downstream |
| 7 | **Link RNA-seq data** | D3 | 1 day | Class A/B evidence |
| 8 | **Evaluate vs. Bassani-Sternberg 2016** | E1–E2 | 2 days | Paper results |

---

## Datasets (per obj3.txt)

| Dataset | Paper | Status | Role |
|:---|:---|:---|:---|
| PXD005231 | Bassani-Sternberg 2016 | ✅ 61 RAW downloaded, 40 MGF converted | Primary training + testing |
| Carreno et al. 2015 | PNAS | 🔲 Not yet downloaded | Validation dataset |
| Laumont et al. 2018 | Nat Commun (PXD008084) | 🔲 Not yet downloaded | Additional training data |
| RNA-seq companion | Bassani-Sternberg 2016 | 🔲 GEO accession to be retrieved | Expression ranking + evaluation |
