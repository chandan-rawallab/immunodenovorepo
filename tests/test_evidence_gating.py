import importlib.util
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).parent.parent
RANK_SCRIPT = ROOT / "src" / "postprocess" / "08_rank_candidates.py"


def load_rank_module():
    spec = importlib.util.spec_from_file_location("rank_candidates", RANK_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def run_rank(tmp_path: Path, rna_source: str) -> pd.DataFrame:
    module = load_rank_module()
    filtered = tmp_path / "filtered.tsv"
    manifest = tmp_path / "manifest.tsv"
    expr = tmp_path / "expr.tsv"
    output = tmp_path / "ranked.tsv"

    filtered.write_text(
        "sample_id\trun_id\tpeptide\tscore\tmutation_type\tsource_protein\n"
        "S1\tRUN1\tACDEFGHIK\t0.9\tmissense\tP00001\n"
    )
    expr.write_text("gene\texpression_tpm\nP00001\t25.0\n")
    manifest.write_text(
        "study_id\trun_id\tpatient_id\thla_alleles\thla_source\trna_expr_path\trna_source\tinclude_in_pipeline\n"
        f"STUDY1\tRUN1\tS1\tHLA-A*02:01\tclinical_ngs\t{expr}\t{rna_source}\tTrue\n"
    )

    def fake_mhcflurry(peptides, alleles, output_path):
        Path(output_path).write_text(
            "allele,peptide,mhcflurry_presentation_percentile\n"
            "HLA-A*02:01,ACDEFGHIK,0.5\n"
        )
        return True

    module.run_mhcflurry = fake_mhcflurry
    argv = sys.argv
    try:
        sys.argv = [
            "08_rank_candidates.py",
            "--input", str(filtered),
            "--manifest", str(manifest),
            "--output", str(output),
            "--binding_rank_cutoff", "2.0",
            "--tpm_cutoff", "1.0",
        ]
        module.main()
    finally:
        sys.argv = argv

    return pd.read_csv(output, sep="\t")


def test_mock_rna_cannot_create_strong_evidence_class(tmp_path: Path):
    ranked = run_rank(tmp_path, "mock_debug")

    assert ranked.loc[0, "evidence_class"] == "C"
    assert ranked.loc[0, "expression_evidence_status"] == "debug_or_non_patient_matched"
    assert not bool(ranked.loc[0, "expression_supports_biology"])
    assert "cannot support" in ranked.loc[0, "evidence_limitations"]


def test_patient_matched_rna_can_support_strong_evidence_class(tmp_path: Path):
    ranked = run_rank(tmp_path, "real_patient_matched")

    assert ranked.loc[0, "evidence_class"] == "A"
    assert ranked.loc[0, "expression_evidence_status"] == "patient_matched"
    assert bool(ranked.loc[0, "expression_supports_biology"])
    assert pd.isna(ranked.loc[0, "evidence_limitations"]) or ranked.loc[0, "evidence_limitations"] == ""
