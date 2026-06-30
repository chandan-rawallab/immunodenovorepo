import os
import subprocess
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).parent.parent
SCRIPT = ROOT / "src" / "data_prep" / "00d_link_expression.py"


def run_linker(tmp_path: Path, extra_args: list[str] | None = None) -> pd.DataFrame:
    manifest = tmp_path / "manifest.tsv"
    output = tmp_path / "expression_matrix.tsv"
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(
        ">sp|P00001|TEST_HUMAN Test protein OS=Homo sapiens OX=9606\n"
        "MPEPTIDESEQ\n",
        encoding="utf-8",
    )
    manifest.write_text(
        """study_id	run_id	patient_id	validation_id	filename	cohort	sample_role	hla_alleles	hla_source	rna_expr_path	rna_source	raw_source	psm_source	include_in_pipeline	notes\nSTUDY1	RUN1	P1	P1	msms_RUN1.raw.txt	STUDY1	hla_peptidome	HLA-A*02:01	manual			local_raw	local_psm	True	\n""",
        encoding="utf-8",
    )

    cmd = [
        sys.executable,
        str(SCRIPT),
        "--manifest",
        str(manifest),
        "--fasta",
        str(fasta),
        "--output",
        str(output),
    ]
    if extra_args:
        cmd.extend(extra_args)

    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    result = subprocess.run(cmd, cwd=tmp_path, capture_output=True, text=True, env=env)
    assert result.returncode == 0, result.stderr
    return pd.read_csv(manifest, sep="\t", dtype=str).fillna("")


def test_linker_is_strict_by_default(tmp_path: Path):
    manifest = run_linker(tmp_path)
    row = manifest.iloc[0]
    assert row["rna_expr_path"] == ""
    assert row["rna_source"] == "missing"
    assert not (tmp_path / "data" / "expression" / "P1_tpm.tsv").exists()


def test_linker_generates_mock_only_in_debug_mode(tmp_path: Path):
    manifest = run_linker(tmp_path, ["--debug-expression"])
    row = manifest.iloc[0]
    assert row["rna_source"] == "mock_debug"
    assert row["rna_expr_path"].endswith("P1_tpm.tsv")
    assert (tmp_path / "data" / "expression" / "P1_tpm.tsv").exists()
