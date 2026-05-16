"""Hybrid Objective 3 pipeline logic.

This module bridges the proposal-faithful custom MaxQuant + CNN/LSTM workflow
with the downstream de novo subtraction, ranking, and report shape used by the
standalone obj3 workspace.
"""

from __future__ import annotations

import json
import math
import re
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

from .io_utils import read_table, sha256_file, write_tsv
from .mgf_utils import iter_mgf

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
CANONICAL_HLA_I = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]{8,11}$")
DEFAULT_PSM_FDR_THRESHOLD = 0.01
DEFAULT_DE_NOVO_SCORE_CUTOFF = 0.70
DEFAULT_BINDING_RANK_SUPPORT_MAX = 2.0
DEFAULT_EXPRESSION_TPM_SUPPORT_MIN = 1.0
EVIDENCE_CLASS_ORDER = {"A": 0, "B": 1, "C": 2}

PSM_FIELDS = ["sample_id", "spectrum_id", "peptide", "q_value", "source_file"]
DE_NOVO_FIELDS = ["sample_id", "spectrum_id", "peptide", "length", "de_novo_score", "source_file"]
RANKED_FIELDS = [
    "sample_id",
    "spectrum_id",
    "peptide",
    "length",
    "source",
    "de_novo_score",
    "psm_support",
    "best_hla",
    "binding_score",
    "binding_rank",
    "expression_tpm",
    "variant_match",
    "evidence_class",
    "rank",
]


def normalize_peptide(raw_value: object) -> str:
    return "".join(str(raw_value or "").strip().upper().split())


def to_float(raw_value: object) -> float | None:
    value = "" if raw_value is None else str(raw_value).strip()
    if not value or value.lower() == "nan":
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    return None if math.isnan(parsed) else parsed


def normalize_maxquant_psms(results_dir: Path, output: Path, fdr_column: str | None = None) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for path in sorted(results_dir.glob("msms_*.txt")):
        sample_id = _sample_from_msms_name(path.name)
        for raw in read_table(path):
            peptide = normalize_peptide(raw.get("Sequence", ""))
            scan = str(raw.get("Scan number", "")).strip()
            if not peptide or not scan:
                continue
            q_value = _pick_confidence_value(raw, fdr_column)
            rows.append(
                {
                    "sample_id": sample_id,
                    "spectrum_id": scan,
                    "peptide": peptide,
                    "q_value": "" if q_value is None else q_value,
                    "source_file": str(path),
                }
            )
    write_tsv(output, rows, PSM_FIELDS)
    return rows


def dataset_status(mgf_dir: Path, results_dir: Path) -> dict[str, object]:
    mgf_files = sorted(mgf_dir.glob("*.mgf"))
    msms_files = sorted(results_dir.glob("msms_*.txt"))
    mgf_samples = {path.stem for path in mgf_files}
    msms_samples = {_sample_from_msms_name(path.name) for path in msms_files}
    total_psm_rows = 0
    for path in msms_files:
        total_psm_rows += max(0, sum(1 for _ in path.open("r", encoding="utf-8", errors="replace")) - 1)
    return {
        "mgf_count": len(mgf_files),
        "msms_count": len(msms_files),
        "matched_count": len(mgf_samples & msms_samples),
        "missing_msms_for_mgf": sorted(mgf_samples - msms_samples),
        "missing_mgf_for_msms": sorted(msms_samples - mgf_samples),
        "raw_psm_rows": total_psm_rows,
    }


def predict_de_novo(
    *,
    mgf_dir: Path,
    checkpoint: Path,
    output: Path,
    max_spectra_per_file: int | None = None,
    min_length: int = 8,
    max_length: int = 11,
    bin_size: float = 0.1,
    max_mz: float = 2000.0,
) -> list[dict[str, object]]:
    import sys
    import numpy as np
    import torch

    src_root = Path(__file__).resolve().parents[1]
    if str(src_root) not in sys.path:
        sys.path.insert(0, str(src_root))
    from dl_module.model import NeoepitopeSeq2Seq

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = NeoepitopeSeq2Seq().to(device)
    model.load_state_dict(torch.load(checkpoint, map_location=device))
    model.eval()

    rows: list[dict[str, object]] = []
    vector_size = int(max_mz / bin_size)
    for mgf_path in sorted(mgf_dir.glob("*.mgf")):
        sample_id = mgf_path.stem
        for index, spectrum in enumerate(iter_mgf(mgf_path), start=1):
            if max_spectra_per_file is not None and index > max_spectra_per_file:
                break
            vector = _bin_spectrum(spectrum.mz_array, spectrum.intensity_array, vector_size, bin_size, max_mz, np)
            tensor = torch.tensor(vector, dtype=torch.float32).unsqueeze(0).unsqueeze(0).to(device)
            with torch.no_grad():
                logits = model(tensor).squeeze(0)
                probs = torch.softmax(logits, dim=-1)
                confidence, indices = torch.max(probs, dim=-1)
            peptide, score = _decode_prediction(indices.cpu().tolist(), confidence.cpu().tolist())
            if not peptide or not (min_length <= len(peptide) <= max_length):
                continue
            if not CANONICAL_HLA_I.fullmatch(peptide):
                continue
            rows.append(
                {
                    "sample_id": sample_id,
                    "spectrum_id": spectrum.spectrum_id,
                    "peptide": peptide,
                    "length": len(peptide),
                    "de_novo_score": round(score, 6),
                    "source_file": str(mgf_path),
                }
            )
    write_tsv(output, rows, DE_NOVO_FIELDS)
    return rows


def postprocess(
    *,
    psm_path: Path,
    de_novo_path: Path,
    outdir: Path,
    metadata_path: Path | None = None,
    binding_input: Path | None = None,
    fasta_path: Path | None = None,
    psm_fdr_threshold: float = DEFAULT_PSM_FDR_THRESHOLD,
    de_novo_score_cutoff: float = DEFAULT_DE_NOVO_SCORE_CUTOFF,
    binding_rank_support_max: float = DEFAULT_BINDING_RANK_SUPPORT_MAX,
    expression_tpm_support_min: float = DEFAULT_EXPRESSION_TPM_SUPPORT_MIN,
) -> tuple[Path, Path, Path, Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    psm_rows = read_table(psm_path)
    de_novo_input = read_table(de_novo_path)
    metadata = _load_metadata(metadata_path)
    candidates = subtract_database_search(psm_rows, de_novo_input, psm_fdr_threshold, de_novo_score_cutoff)
    ranked, warnings = rank_candidates(
        candidates,
        metadata,
        binding_input,
        binding_rank_support_max,
        expression_tpm_support_min,
    )

    de_novo_out = outdir / "de_novo_candidates.tsv"
    ranked_out = outdir / "ranked_neoantigens.tsv"
    report_md = outdir / "run_report.md"
    report_json = outdir / "run_report.json"
    write_tsv(de_novo_out, candidates, DE_NOVO_FIELDS)
    write_tsv(ranked_out, ranked, RANKED_FIELDS)
    _write_report(
        report_md,
        report_json,
        psm_path=psm_path,
        de_novo_path=de_novo_path,
        metadata_path=metadata_path,
        fasta_path=fasta_path,
        psm_rows=psm_rows,
        de_novo_rows=candidates,
        ranked_rows=ranked,
        warnings=warnings,
        thresholds={
            "psm_fdr_threshold": psm_fdr_threshold,
            "de_novo_score_cutoff": de_novo_score_cutoff,
            "binding_rank_support_max": binding_rank_support_max,
            "expression_tpm_support_min": expression_tpm_support_min,
        },
    )
    return de_novo_out, ranked_out, report_md, report_json


def subtract_database_search(
    psm_rows: list[dict[str, object]],
    de_novo_rows: list[dict[str, object]],
    psm_fdr_threshold: float,
    de_novo_score_cutoff: float,
) -> list[dict[str, object]]:
    accepted: dict[str, set[str]] = defaultdict(set)
    for row in psm_rows:
        q_value = to_float(row.get("q_value"))
        if q_value is None or q_value <= psm_fdr_threshold:
            accepted[str(row["sample_id"])].add(str(row["peptide"]))

    candidates: list[dict[str, object]] = []
    for row in de_novo_rows:
        sample_id = str(row.get("sample_id", ""))
        peptide = normalize_peptide(row.get("peptide", ""))
        score = to_float(row.get("de_novo_score"))
        if not sample_id or not CANONICAL_HLA_I.fullmatch(peptide):
            continue
        if peptide in accepted.get(sample_id, set()):
            continue
        if score is not None and score < de_novo_score_cutoff:
            continue
        candidates.append(
            {
                "sample_id": sample_id,
                "spectrum_id": str(row.get("spectrum_id", "")),
                "peptide": peptide,
                "length": len(peptide),
                "de_novo_score": "" if score is None else score,
                "source_file": str(row.get("source_file", "")),
            }
        )
    return candidates


def rank_candidates(
    candidates: list[dict[str, object]],
    metadata: dict[str, dict[str, object]],
    binding_input: Path | None,
    binding_rank_support_max: float,
    expression_tpm_support_min: float,
) -> tuple[list[dict[str, object]], list[str]]:
    warnings: list[str] = []
    expression = _load_expression(metadata)
    variants = _load_variants(metadata)
    binding = _load_binding(binding_input)
    if not binding:
        scored, scoring_warnings = _score_with_mhcflurry(candidates, metadata)
        binding.update(scored)
        warnings.extend(scoring_warnings)

    support_counts = Counter((row["sample_id"], row["peptide"]) for row in candidates)
    ranked: list[dict[str, object]] = []
    for row in candidates:
        sample_id = str(row["sample_id"])
        peptide = str(row["peptide"])
        variant = variants.get(sample_id, {}).get(peptide, {})
        bind = binding.get((sample_id, peptide), {})
        binding_rank = to_float(bind.get("binding_rank"))
        binding_score = to_float(bind.get("binding_score"))
        expression_tpm = _expression_for(sample_id, peptide, str(variant.get("gene", "")), expression)
        has_binding = binding_rank is not None and binding_rank <= binding_rank_support_max
        has_expression = expression_tpm is not None and expression_tpm >= expression_tpm_support_min
        mutation = str(variant.get("mutation", ""))
        evidence_class = "A" if mutation and has_binding and has_expression else "B" if has_binding and has_expression else "C"
        ranked.append(
            {
                "sample_id": sample_id,
                "spectrum_id": row["spectrum_id"],
                "peptide": peptide,
                "length": row["length"],
                "source": "de_novo",
                "de_novo_score": row["de_novo_score"],
                "psm_support": support_counts[(sample_id, peptide)],
                "best_hla": bind.get("best_hla", ""),
                "binding_score": "" if binding_score is None else binding_score,
                "binding_rank": "" if binding_rank is None else binding_rank,
                "expression_tpm": "" if expression_tpm is None else expression_tpm,
                "variant_match": mutation,
                "evidence_class": evidence_class,
                "rank": 0,
            }
        )

    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in ranked:
        grouped[str(row["sample_id"])].append(row)
    final_rows: list[dict[str, object]] = []
    for _, rows in grouped.items():
        rows.sort(
            key=lambda row: (
                EVIDENCE_CLASS_ORDER[str(row["evidence_class"])],
                -int(row["psm_support"]),
                to_float(row["binding_rank"]) if to_float(row["binding_rank"]) is not None else float("inf"),
                -(to_float(row["expression_tpm"]) or 0.0),
                -(to_float(row["de_novo_score"]) or 0.0),
                str(row["spectrum_id"]),
            )
        )
        for rank, row in enumerate(rows, start=1):
            row["rank"] = rank
            final_rows.append(row)
    return final_rows, warnings


def _sample_from_msms_name(name: str) -> str:
    sample = name
    if sample.startswith("msms_"):
        sample = sample[len("msms_") :]
    if sample.endswith(".txt"):
        sample = sample[:-4]
    if sample.endswith(".raw"):
        sample = sample[:-4]
    return sample


def _pick_confidence_value(row: dict[str, str], preferred: str | None) -> float | None:
    candidates = [preferred] if preferred else []
    candidates.extend(["Q-value", "q_value", "Q Value", "PEP", "Posterior Error Probability"])
    lowered = {key.lower(): value for key, value in row.items()}
    for candidate in candidates:
        if candidate and candidate.lower() in lowered:
            value = to_float(lowered[candidate.lower()])
            if value is not None:
                return value
    return None


def _bin_spectrum(mz_array, intensity_array, vector_size, bin_size, max_mz, np):
    vector = np.zeros(vector_size, dtype=np.float32)
    for mz, intensity in zip(mz_array, intensity_array):
        if 0 <= mz < max_mz:
            idx = int(mz / bin_size)
            if idx < vector_size:
                vector[idx] += intensity
    max_intensity = float(np.max(vector)) if vector.size else 0.0
    if max_intensity > 0:
        vector = vector / max_intensity
    return vector


def _decode_prediction(indices: list[int], confidence: list[float]) -> tuple[str, float]:
    int_to_aa = {idx + 1: aa for idx, aa in enumerate(AMINO_ACIDS)}
    peptide: list[str] = []
    residue_conf: list[float] = []
    for idx, prob in zip(indices, confidence):
        if idx == 22:
            break
        aa = int_to_aa.get(idx)
        if aa:
            peptide.append(aa)
            residue_conf.append(float(prob))
    if not peptide:
        return "", 0.0
    return "".join(peptide), sum(residue_conf) / len(residue_conf)


def _load_metadata(path: Path | None) -> dict[str, dict[str, object]]:
    if not path or not path.exists():
        return {}
    base = path.parent
    metadata: dict[str, dict[str, object]] = {}
    for row in read_table(path):
        sample_id = str(row.get("sample_id", "")).strip()
        if not sample_id:
            continue
        row = dict(row)
        for key in ["rna_expr_path", "variant_candidates_path"]:
            value = str(row.get(key, "")).strip()
            row[key] = str((base / value).resolve()) if value and not Path(value).is_absolute() else value
        metadata[sample_id] = row
    return metadata


def _load_expression(metadata: dict[str, dict[str, object]]) -> dict[str, dict[str, dict[str, float]]]:
    output: dict[str, dict[str, dict[str, float]]] = {}
    for sample_id, row in metadata.items():
        peptide_map: dict[str, float] = {}
        gene_map: dict[str, float] = {}
        path = Path(str(row.get("rna_expr_path", "")))
        if path.exists():
            for expr in read_table(path):
                tpm = to_float(expr.get("expression_tpm") or expr.get("TPM") or expr.get("tpm"))
                if tpm is None:
                    continue
                peptide = normalize_peptide(expr.get("peptide") or expr.get("Peptide"))
                gene = str(expr.get("gene") or expr.get("Gene") or "").strip()
                if peptide:
                    peptide_map[peptide] = tpm
                if gene:
                    gene_map[gene] = tpm
        output[sample_id] = {"peptide": peptide_map, "gene": gene_map}
    return output


def _load_variants(metadata: dict[str, dict[str, object]]) -> dict[str, dict[str, dict[str, str]]]:
    output: dict[str, dict[str, dict[str, str]]] = {}
    for sample_id, row in metadata.items():
        variant_map: dict[str, dict[str, str]] = {}
        path = Path(str(row.get("variant_candidates_path", "")))
        if path.exists():
            for variant in read_table(path):
                peptide = normalize_peptide(variant.get("peptide") or variant.get("Peptide"))
                if not peptide:
                    continue
                variant_map[peptide] = {
                    "gene": str(variant.get("gene") or variant.get("Gene") or ""),
                    "mutation": str(
                        variant.get("mutation")
                        or variant.get("Mutation")
                        or variant.get("amino_acid_change")
                        or variant.get("variant")
                        or ""
                    ),
                }
        output[sample_id] = variant_map
    return output


def _load_binding(path: Path | None) -> dict[tuple[str, str], dict[str, object]]:
    if not path or not path.exists():
        return {}
    best: dict[tuple[str, str], dict[str, object]] = {}
    for row in read_table(path):
        sample_id = str(row.get("sample_id") or row.get("Sample") or "").strip()
        peptide = normalize_peptide(row.get("peptide") or row.get("Peptide"))
        if not sample_id or not peptide:
            continue
        score = to_float(row.get("binding_score") or row.get("presentation_score") or row.get("affinity"))
        rank = to_float(row.get("binding_rank") or row.get("presentation_percentile") or row.get("affinity_percentile"))
        current = best.get((sample_id, peptide))
        if current is None or _binding_sort_key(rank, score) < _binding_sort_key(to_float(current.get("binding_rank")), to_float(current.get("binding_score"))):
            best[(sample_id, peptide)] = {
                "best_hla": str(row.get("best_hla") or row.get("allele") or ""),
                "binding_score": "" if score is None else score,
                "binding_rank": "" if rank is None else rank,
            }
    return best


def _score_with_mhcflurry(candidates: list[dict[str, object]], metadata: dict[str, dict[str, object]]):
    predictor = shutil.which("mhcflurry-predict")
    if not predictor:
        return {}, ["mhcflurry-predict not found; binding columns were left blank."]
    results: dict[tuple[str, str], dict[str, object]] = {}
    grouped: dict[str, set[str]] = defaultdict(set)
    for row in candidates:
        grouped[str(row["sample_id"])].add(str(row["peptide"]))
    with tempfile.TemporaryDirectory(prefix="objective3_mhcflurry_") as tmpdir:
        tmp = Path(tmpdir)
        for sample_id, peptides in grouped.items():
            alleles = [
                _normalize_hla(str(metadata.get(sample_id, {}).get(column, "")))
                for column in ["hla_a1", "hla_a2", "hla_b1", "hla_b2", "hla_c1", "hla_c2"]
            ]
            alleles = [allele for allele in alleles if allele]
            if not alleles:
                continue
            input_path = tmp / f"{sample_id}_binding_input.csv"
            output_path = tmp / f"{sample_id}_binding_output.csv"
            input_path.write_text(
                "allele,peptide\n"
                + "".join(f"{allele},{peptide}\n" for allele in alleles for peptide in sorted(peptides)),
                encoding="utf-8",
            )
            try:
                subprocess.run([predictor, str(input_path), "--out", str(output_path)], check=True, capture_output=True, text=True)
            except subprocess.CalledProcessError as exc:
                return results, [f"MHCflurry failed for sample {sample_id}: {exc.stderr.strip() or exc.stdout.strip()}"]
            for row in _read_csv_as_table(output_path):
                peptide = normalize_peptide(row.get("peptide"))
                allele = str(row.get("allele") or row.get("Allele") or "")
                score = to_float(row.get("presentation_score") or row.get("binding_score") or row.get("affinity"))
                rank = to_float(row.get("presentation_percentile") or row.get("affinity_percentile") or row.get("binding_rank"))
                key = (sample_id, peptide)
                current = results.get(key)
                if current is None or _binding_sort_key(rank, score) < _binding_sort_key(to_float(current.get("binding_rank")), to_float(current.get("binding_score"))):
                    results[key] = {"best_hla": allele, "binding_score": "" if score is None else score, "binding_rank": "" if rank is None else rank}
    return results, []


def _read_csv_as_table(path: Path) -> list[dict[str, str]]:
    import csv

    with path.open("r", encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _normalize_hla(value: str) -> str:
    value = value.strip()
    if not value:
        return ""
    return value.replace("HLA-", "").replace("*", "").replace(":", "")


def _binding_sort_key(rank: float | None, score: float | None) -> tuple[float, float]:
    return (rank if rank is not None else float("inf"), -(score if score is not None else float("-inf")))


def _expression_for(sample_id: str, peptide: str, gene: str, expression) -> float | None:
    sample = expression.get(sample_id, {})
    peptide_map = sample.get("peptide", {})
    gene_map = sample.get("gene", {})
    if peptide in peptide_map:
        return peptide_map[peptide]
    if gene and gene in gene_map:
        return gene_map[gene]
    return None


def _write_report(
    markdown_path: Path,
    json_path: Path,
    *,
    psm_path: Path,
    de_novo_path: Path,
    metadata_path: Path | None,
    fasta_path: Path | None,
    psm_rows: list[dict[str, object]],
    de_novo_rows: list[dict[str, object]],
    ranked_rows: list[dict[str, object]],
    warnings: list[str],
    thresholds: dict[str, object],
) -> None:
    class_counts = Counter(row.get("evidence_class", "") for row in ranked_rows)
    payload = {
        "summary": {
            "psm_input": str(psm_path),
            "de_novo_input": str(de_novo_path),
            "sample_metadata": str(metadata_path) if metadata_path else "",
            "num_psms": len(psm_rows),
            "num_de_novo_candidates": len(de_novo_rows),
            "num_ranked_candidates": len(ranked_rows),
            "evidence_class_counts": dict(class_counts),
        },
        "thresholds": thresholds,
        "warnings": warnings,
        "plain_language_explanation": (
            "This run combines database-search supported immunopeptidome peptides with "
            "custom CNN/LSTM de novo candidates, removes exact overlaps, and ranks the "
            "remaining candidates using available binding, expression, and variant evidence."
        ),
    }
    if fasta_path and fasta_path.exists():
        payload["reference_fasta"] = {"path": str(fasta_path), "sha256": sha256_file(fasta_path)}

    markdown = f"""# Objective 3 Run Report

## Summary

- PSM input: `{payload["summary"]["psm_input"]}`
- De novo input: `{payload["summary"]["de_novo_input"]}`
- Sample metadata: `{payload["summary"]["sample_metadata"]}`
- Immunopeptidome PSM rows: `{payload["summary"]["num_psms"]}`
- De novo candidate rows: `{payload["summary"]["num_de_novo_candidates"]}`
- Ranked neoantigen rows: `{payload["summary"]["num_ranked_candidates"]}`

## Evidence class counts

{json.dumps(payload["summary"]["evidence_class_counts"], indent=2)}

## Thresholds

{json.dumps(thresholds, indent=2)}

## Plain-language explanation

{payload["plain_language_explanation"]}

## Warnings

{json.dumps(warnings, indent=2)}
"""
    markdown_path.write_text(markdown, encoding="utf-8")
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
