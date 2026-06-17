#!/usr/bin/env bash
# run_pipeline.sh — End-to-end neoepitope pipeline runner
# Usage:
#   bash run_pipeline.sh <ACCESSION_ID> [--dry-run] [--patient PATIENT_ID] [--hla-override PATH] [--expression PATH] [--ega-credentials PATH]
#
# Runs: data_acquire → build_manifest → extract_psms → autotype_hla → link_expression → 
#       predict_denovo → filter_neoantigens → rank_candidates → evaluate_neoantigens

set -euo pipefail

# ─── Config ────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$SCRIPT_DIR"
SRC="$ROOT/src"
CONFIGS="$ROOT/configs"
DATA="$ROOT/data"
RESULTS="$ROOT/results"
LOGS="$ROOT/logs"

MANIFEST="$CONFIGS/sample_manifest.tsv"
FASTA="$DATA/reference/uniprot_human_reviewed.only_human.fasta"
VALIDATED="$DATA/reference/s2_dataset_extracted/Dataset1/Dataset1.txt"
CHECKPOINT="$SRC/results/checkpoints/neoepitope_production_best.pth"
EXPRESSION_MATRIX="$DATA/expression_matrix.tsv"
PSMS_FILE="$RESULTS/immunopeptidome_psms.tsv"

DATESTAMP=$(date +%Y%m%d_%H%M%S)
RUN_DIR="$RESULTS/${DATESTAMP}_denovo_run"
mkdir -p "$RUN_DIR" "$LOGS" "$CONFIGS" "$DATA/raw" "$DATA/psms" "$DATA/mgf" "$DATA/mgf_unlabeled" "$RESULTS"

LOG_FILE="$LOGS/${DATESTAMP}_pipeline.log"
DRY_RUN=false
PATIENT_FILTER=""
ACCESSION=""
HLA_OVERRIDE=""
EXPRESSION_OVERRIDE=""
EGA_CREDS=""

# ─── Arg parsing ────────────────────────────────────────────────────────────
if [[ $# -eq 0 ]]; then
    echo "Error: Must provide an accession ID (e.g. PXD005231) or --local"
    exit 1
fi

ACCESSION="$1"
shift

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=true ;;
        --patient) PATIENT_FILTER="$2"; shift ;;
        --hla-override) HLA_OVERRIDE="$2"; shift ;;
        --expression) EXPRESSION_OVERRIDE="$2"; shift ;;
        --ega-credentials) EGA_CREDS="$2"; shift ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
    shift
done

run_validation() {
    local step_name="$1"
    shift

    echo ">>> [validate:${step_name}] $*" | tee -a "$LOG_FILE"
    if [[ "$DRY_RUN" == "false" ]]; then
        python3 "$SRC/validation/preflight_validate.py" "$@" 2>&1 | tee -a "$LOG_FILE"
    else
        echo "[DRY-RUN] validation skipped" | tee -a "$LOG_FILE"
    fi
}

run_provenance_audit() {
    echo ">>> [provenance_audit] Writing provenance audit for this run" | tee -a "$LOG_FILE"
    if [[ "$DRY_RUN" == "false" ]]; then
        python3 "$SRC/validation/provenance_audit.py" \
            --accession "$ACCESSION" \
            --manifest "$MANIFEST" \
            --excluded-runs "$CONFIGS/excluded_runs.tsv" \
            --psm-dir "$DATA/psms" \
            --mgf-dir "$DATA/mgf" \
            --unlabeled-mgf-dir "$DATA/mgf_unlabeled" \
            --reference-fasta "$FASTA" \
            --validated "$VALIDATED" \
            --de-novo "$CANDIDATES" \
            --filtered "$FILTERED" \
            --ranked "$RANKED" \
            $(expected_provenance_args) \
            --output-md "$RUN_DIR/provenance_audit.md" \
            --output-json "$RUN_DIR/provenance_audit.json" 2>&1 | tee -a "$LOG_FILE"
    else
        echo "[DRY-RUN] provenance audit skipped" | tee -a "$LOG_FILE"
    fi
}

expected_run_args() {
    if [[ "$ACCESSION" == "PXD005231" ]]; then
        echo "--expected-active-runs 31"
    fi
}

expected_provenance_args() {
    if [[ "$ACCESSION" == "PXD005231" ]]; then
        echo "--expected-active-runs 31 --expected-excluded-runs 9"
    fi
}

run_step() {
    local step_name="$1"
    shift
    local marker="$RUN_DIR/.done_${step_name}"

    echo ">>> [$step_name] $*" | tee -a "$LOG_FILE"
    
    if [[ -f "$marker" ]]; then
        echo "[SKIPPED] Step $step_name already completed (marker found: $marker)" | tee -a "$LOG_FILE"
        return 0
    fi

    if [[ "$DRY_RUN" == "false" ]]; then
        python3 "$@" 2>&1 | tee -a "$LOG_FILE"
        # Create marker on success (set -e ensures we don't reach this on failure)
        touch "$marker"
    else
        echo "[DRY-RUN] skipped"
    fi
}

echo "========================================" | tee "$LOG_FILE"
echo " Neoepitope Pipeline — ${DATESTAMP}" | tee -a "$LOG_FILE"
echo " Accession: $ACCESSION" | tee -a "$LOG_FILE"
echo " Run dir: $RUN_DIR" | tee -a "$LOG_FILE"
echo " Dry run: $DRY_RUN" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"

# ─── Step 00: Data Acquisition ──────────────────────────────────────────────
echo "" | tee -a "$LOG_FILE"
echo "── STEP 00: Smart Data Acquisition ──" | tee -a "$LOG_FILE"

ACQUIRE_CMD=("$SRC/data_prep/00_acquire_data.py" "--accession" "$ACCESSION" "--raw-dir" "$DATA/raw" "--psm-dir" "$DATA/psms")
if [[ -n "$EGA_CREDS" ]]; then
    ACQUIRE_CMD+=("--ega-credentials" "$EGA_CREDS")
fi
PYTHONPATH="$SRC" run_step "00_acquire" "${ACQUIRE_CMD[@]}"

# ─── Step 00b: Build Manifest ───────────────────────────────────────────────
echo "" | tee -a "$LOG_FILE"
echo "── STEP 00b: Build Manifest ──" | tee -a "$LOG_FILE"

MANIFEST_CMD=("$SRC/data_prep/00b_build_manifest.py" "--accession" "$ACCESSION" "--raw-dir" "$DATA/raw" "--psm-dir" "$DATA/psms" "--output" "$MANIFEST")
if [[ -n "$HLA_OVERRIDE" ]]; then
    MANIFEST_CMD+=("--hla-override" "$HLA_OVERRIDE")
fi
PYTHONPATH="$SRC" run_step "00b_manifest" "${MANIFEST_CMD[@]}"
# ─── Step 04: Extract PSMs ──────────────────────────────────────────────────
echo "" | tee -a "$LOG_FILE"
echo "── STEP 04: Extract PSMs from Search Results ──" | tee -a "$LOG_FILE"
PYTHONPATH="$SRC" run_step "04_extract" "$SRC/data_prep/04_extract_psms.py" \
    --input-dir "$DATA/psms" \
    --output-file "$PSMS_FILE" \
    --manifest "$MANIFEST"
# ─── Step 00c: Auto-type HLA ────────────────────────────────────────────────
echo "" | tee -a "$LOG_FILE"
echo "── STEP 00c: HLA Auto-typing ──" | tee -a "$LOG_FILE"
PYTHONPATH="$SRC" run_step "00c_autotype" "$SRC/data_prep/00c_autotype_hla.py" \
    --manifest "$MANIFEST" \
    --psms "$PSMS_FILE" \
    --output "$MANIFEST"

# ─── Step 00d: Link Expression ──────────────────────────────────────────────
echo "" | tee -a "$LOG_FILE"
echo "── STEP 00d: Link Expression ──" | tee -a "$LOG_FILE"

EXPR_CMD=("$SRC/data_prep/00d_link_expression.py" "--accession" "$ACCESSION" "--manifest" "$MANIFEST" "--output" "$EXPRESSION_MATRIX")
if [[ -n "$EXPRESSION_OVERRIDE" ]]; then
    EXPR_CMD+=("--expression" "$EXPRESSION_OVERRIDE")
fi
PYTHONPATH="$SRC" run_step "00d_expression" "${EXPR_CMD[@]}"
run_validation "manifest" \
    --manifest "$MANIFEST" \
    --psm-dir "$DATA/psms" \
    --mgf-dir "$DATA/mgf" \
    --unlabeled-mgf-dir "$DATA/mgf_unlabeled" \
    --reference-fasta "$FASTA" \
    --psms "$PSMS_FILE" \
    $(expected_run_args) \
    --output "$RUN_DIR/preflight_manifest.md"

# ─── Step 05: Extract Unlabeled Spectra ──────────────────────────────────────────
echo "" | tee -a "$LOG_FILE"
echo "── STEP 05: Extract Unlabeled Spectra ──" | tee -a "$LOG_FILE"
PYTHONPATH="$SRC" run_step "05_extract_unlabeled" "$SRC/data_prep/05_extract_unlabeled_spectra.py" \
    --mgf-dir "$DATA/mgf" \
    --psms    "$PSMS_FILE" \
    --output-dir "$DATA/mgf_unlabeled" \
    --manifest "$MANIFEST"
run_validation "unlabeled" \
    --manifest "$MANIFEST" \
    --psm-dir "$DATA/psms" \
    --mgf-dir "$DATA/mgf" \
    --unlabeled-mgf-dir "$DATA/mgf_unlabeled" \
    --reference-fasta "$FASTA" \
    --psms "$PSMS_FILE" \
    $(expected_run_args) \
    --output "$RUN_DIR/preflight_unlabeled.md"

# ─── Step 06: De novo prediction ────────────────────────────────────────────
echo "" | tee -a "$LOG_FILE"
echo "── STEP 06: De novo prediction ──" | tee -a "$LOG_FILE"
CANDIDATES="$RUN_DIR/de_novo_candidates.tsv"

PYTHONPATH="$SRC" run_step "06_predict" "$SRC/inference/06_predict_denovo.py" \
    --model   "$CHECKPOINT" \
    --mgf_dir "$DATA/mgf_unlabeled" \
    --output  "$CANDIDATES"

# ─── Step 07: Filter neoantigens ────────────────────────────────────────────
echo "" | tee -a "$LOG_FILE"
echo "── STEP 07: Filter neoantigens ──" | tee -a "$LOG_FILE"
FILTERED="$RUN_DIR/filtered_neoantigens.tsv"

PYTHONPATH="$SRC" run_step "07_filter" "$SRC/postprocess/07_filter_neoantigens.py" \
    --input  "$CANDIDATES" \
    --psms   "$PSMS_FILE" \
    --fasta  "$FASTA" \
    --manifest "$MANIFEST" \
    --output "$FILTERED" \
    --score_cutoff -0.5 \
    --min_psm_support 2

# ─── Step 08: Rank candidates ───────────────────────────────────────────────
echo "" | tee -a "$LOG_FILE"
echo "── STEP 08: Rank candidates ──" | tee -a "$LOG_FILE"
RANKED="$RUN_DIR/ranked_neoantigens.tsv"

PYTHONPATH="$SRC" run_step "08_rank" "$SRC/postprocess/08_rank_candidates.py" \
    --input    "$FILTERED" \
    --manifest "$MANIFEST" \
    --output   "$RANKED" \
    --binding_rank_cutoff 2.0 \
    --tpm_cutoff 1.0
run_validation "ranked" \
    --manifest "$MANIFEST" \
    --psm-dir "$DATA/psms" \
    --mgf-dir "$DATA/mgf" \
    --unlabeled-mgf-dir "$DATA/mgf_unlabeled" \
    --reference-fasta "$FASTA" \
    --psms "$PSMS_FILE" \
    --de-novo "$CANDIDATES" \
    --filtered "$FILTERED" \
    --ranked "$RANKED" \
    $(expected_run_args) \
    --output "$RUN_DIR/preflight_ranked.md"

# ─── Step 09: Evaluate ──────────────────────────────────────────────────────
echo "" | tee -a "$LOG_FILE"
echo "── STEP 09: Evaluate ──" | tee -a "$LOG_FILE"
REPORT="$RUN_DIR/evaluation_report.md"

PYTHONPATH="$SRC" run_step "09_evaluate" "$SRC/evaluation/09_evaluate_neoantigens.py" \
    --input     "$RANKED" \
    --validated "$VALIDATED" \
    --manifest  "$MANIFEST" \
    --output    "$REPORT" \
    --top_n 10 25 50
run_provenance_audit

echo "" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo " Pipeline complete." | tee -a "$LOG_FILE"
echo " Results: $RUN_DIR" | tee -a "$LOG_FILE"
echo " Log:     $LOG_FILE" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
