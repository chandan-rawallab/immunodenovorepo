from pathlib import Path

from objective3.io_utils import read_table, write_tsv
from objective3.pipeline import normalize_maxquant_psms, postprocess


def test_normalize_maxquant_psms_uses_pep_as_q_value(tmp_path: Path):
    results_dir = tmp_path / "production_results"
    results_dir.mkdir()
    msms = results_dir / "msms_SAMPLE_01.raw.txt"
    msms.write_text(
        "Raw file\tSequence\tPEP\tScan number\n"
        "SAMPLE_01\tACDEFGHIK\t0.002\t42\n",
        encoding="utf-8",
    )
    output = tmp_path / "immunopeptidome_psms.tsv"

    rows = normalize_maxquant_psms(results_dir, output)

    assert rows == [
        {
            "sample_id": "SAMPLE_01",
            "spectrum_id": "42",
            "peptide": "ACDEFGHIK",
            "q_value": 0.002,
            "source_file": str(msms),
        }
    ]
    written = read_table(output)
    assert written[0]["sample_id"] == "SAMPLE_01"
    assert written[0]["q_value"] == "0.002"


def test_postprocess_subtracts_database_hits_and_ranks_with_evidence(tmp_path: Path):
    psms = tmp_path / "psms.tsv"
    denovo = tmp_path / "denovo.tsv"
    metadata = tmp_path / "metadata.tsv"
    expression = tmp_path / "expr.tsv"
    variants = tmp_path / "variants.tsv"
    binding = tmp_path / "binding.tsv"
    outdir = tmp_path / "final"

    write_tsv(
        psms,
        [{"sample_id": "S1", "spectrum_id": "1", "peptide": "ACDEFGHIK", "q_value": "0.001", "source_file": "x"}],
        ["sample_id", "spectrum_id", "peptide", "q_value", "source_file"],
    )
    write_tsv(
        denovo,
        [
            {"sample_id": "S1", "spectrum_id": "1", "peptide": "ACDEFGHIK", "length": 9, "de_novo_score": 0.95, "source_file": "x"},
            {"sample_id": "S1", "spectrum_id": "2", "peptide": "KLGGALQAK", "length": 9, "de_novo_score": 0.91, "source_file": "x"},
        ],
        ["sample_id", "spectrum_id", "peptide", "length", "de_novo_score", "source_file"],
    )
    write_tsv(
        expression,
        [{"peptide": "KLGGALQAK", "gene": "GENE1", "expression_tpm": 8.0}],
        ["peptide", "gene", "expression_tpm"],
    )
    write_tsv(
        variants,
        [{"peptide": "KLGGALQAK", "gene": "GENE1", "mutation": "R175H"}],
        ["peptide", "gene", "mutation"],
    )
    write_tsv(
        binding,
        [{"sample_id": "S1", "peptide": "KLGGALQAK", "best_hla": "HLA-A*02:01", "binding_score": 0.9, "binding_rank": 0.5}],
        ["sample_id", "peptide", "best_hla", "binding_score", "binding_rank"],
    )
    write_tsv(
        metadata,
        [
            {
                "sample_id": "S1",
                "hla_a1": "HLA-A*02:01",
                "hla_a2": "",
                "hla_b1": "",
                "hla_b2": "",
                "hla_c1": "",
                "hla_c2": "",
                "rna_expr_path": str(expression),
                "variant_candidates_path": str(variants),
            }
        ],
        ["sample_id", "hla_a1", "hla_a2", "hla_b1", "hla_b2", "hla_c1", "hla_c2", "rna_expr_path", "variant_candidates_path"],
    )

    _, ranked_path, report_md, report_json = postprocess(
        psm_path=psms,
        de_novo_path=denovo,
        outdir=outdir,
        metadata_path=metadata,
        binding_input=binding,
    )

    ranked = read_table(ranked_path)
    assert len(ranked) == 1
    assert ranked[0]["peptide"] == "KLGGALQAK"
    assert ranked[0]["evidence_class"] == "A"
    assert report_md.exists()
    assert report_json.exists()
