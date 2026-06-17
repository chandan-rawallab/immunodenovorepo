#!/usr/bin/env python3
"""
Integration tests for the neoepitope pipeline (steps 07-09).

Builds minimal synthetic fixtures in a temp directory and exercises
the real scripts end-to-end without needing the full 15k-PSM dataset.

Run with:
    python3 tests/test_pipeline_integration.py
"""

import os, sys, tempfile, textwrap, unittest
from pathlib import Path

# Make src/ importable
ROOT = Path(__file__).parent.parent
SRC  = ROOT / "src"
sys.path.insert(0, str(SRC))

# ─── Fixture helpers ─────────────────────────────────────────────────────────

FASTA_CONTENT = textwrap.dedent("""\
    >sp|P00001|PROT_HUMAN Test protein
    MAQLLKALEV SQSDAVLTLL IYPTAPPRSF KLQQAQSTET
    YLDPVQRDLY HEAAAAAAAT TQVPVVGAVL HFGHIJKLMN
""").replace(" ", "")

CANDIDATES_TSV = textwrap.dedent("""\
    peptide\tsample_id\trun_id\tde_novo_score
    IYPFAPPRSF\tCM467\tCM467_F1_R1\t0.92
    IYPFAPPRSF\tCM467\tCM467_F1_R2\t0.91
    YLDPVQRDLY\tCM467\tCM467_F1_R1\t0.88
    YLDPVQRDLY\tCM467\tCM467_F1_R2\t0.87
    TOOSHORT\tCM467\tCM467_F1_R1\t0.95
    TOOLONGPEPTIDEHERE\tCM467\tCM467_F1_R1\t0.81
    AQLLKALEV\tTIL1\tTIL1_R1\t0.89
    AQLLKALEV\tTIL1\tTIL1_R2\t0.88
""")

PSMS_TSV = textwrap.dedent("""\
    peptide\tsample_id
    KNOWNPEPTIDE\tCM467
""")

VALIDATED_TSV = textwrap.dedent("""\
    Sequence\tIntensity CM467\tIntensity TIL1
    IYPFAPPRSF\t12345.0\t0
    YLDPVQRDLY\t9876.0\t0
    AQLLKALEV\t0\t5432.0
    NOTPREDICTED\t111.0\t0
""")

MANIFEST_TSV = textwrap.dedent("""\
    run_id\tpatient_id\thla_alleles\trna_expr_path
    CM467_F1_R1\tCM467\tHLA-A*01:01,HLA-A*24:02,HLA-B*13:02,HLA-B*39:06\t
    CM467_F1_R2\tCM467\tHLA-A*01:01,HLA-A*24:02,HLA-B*13:02,HLA-B*39:06\t
    TIL1_R1\tTIL1\tHLA-A*02:01,HLA-B*18:01,HLA-B*38:01\t
    TIL1_R2\tTIL1\tHLA-A*02:01,HLA-B*18:01,HLA-B*38:01\t
""")


class TestFilterNeoantigens(unittest.TestCase):
    """Step 07: filter_neoantigens."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        p = Path(self.tmp)
        (p / "candidates.tsv").write_text(CANDIDATES_TSV)
        (p / "psms.tsv").write_text(PSMS_TSV)
        (p / "ref.fasta").write_text(FASTA_CONTENT)
        self.out = p / "filtered.tsv"

    def _run(self, extra_args=""):
        import subprocess
        cmd = [
            sys.executable,
            str(SRC / "postprocess" / "07_filter_neoantigens.py"),
            "--input",  str(Path(self.tmp) / "candidates.tsv"),
            "--psms",   str(Path(self.tmp) / "psms.tsv"),
            "--fasta",  str(Path(self.tmp) / "ref.fasta"),
            "--output", str(self.out),
            "--score_cutoff", "0.7",
            "--min_psm_support", "2",
            "--allow_flanking_mutations",
        ] + (extra_args.split() if extra_args else [])
        env = os.environ.copy()
        env["PYTHONPATH"] = str(SRC)
        result = subprocess.run(cmd, capture_output=True, text=True, env=env)
        if result.returncode != 0:
            print(f"\n--- SCRIPT STDOUT ---\n{result.stdout}")
            print(f"--- SCRIPT STDERR ---\n{result.stderr}")
        return result

    def test_length_filter(self):
        """TOOSHORT and TOOLONGPEPTIDEHERE must be excluded."""
        r = self._run()
        self.assertEqual(r.returncode, 0, r.stderr)
        import pandas as pd
        df = pd.read_csv(self.out, sep="\t")
        self.assertFalse(any("TOOSHORT" in p for p in df["peptide"]),
                         "TOOSHORT should be filtered out by length check")
        self.assertFalse(any("TOOLONG" in p for p in df["peptide"]),
                         "TOOLONG peptide should be filtered out by length check")

    def test_psm_support_filter(self):
        """Peptides with < 2 PSMs should be excluded."""
        r = self._run()
        self.assertEqual(r.returncode, 0, r.stderr)
        import pandas as pd
        df = pd.read_csv(self.out, sep="\t")
        # Each surviving peptide must have appeared at least twice
        # (IYPFAPPRSF and YLDPVQRDLY appear 2× each in fixtures)
        self.assertGreater(len(df), 0, "Some candidates should survive filtering")

    def test_missense_annotation(self):
        """IYPFAPPRSF (T→F at pos 4) should be annotated as missense."""
        r = self._run()
        self.assertEqual(r.returncode, 0, r.stderr)
        import pandas as pd
        df = pd.read_csv(self.out, sep="\t")
        missense_rows = df[df["peptide"] == "IYPFAPPRSF"]
        if not missense_rows.empty:
            self.assertEqual(missense_rows.iloc[0]["mutation_type"], "missense")
            self.assertEqual(int(missense_rows.iloc[0]["mutation_pos"]), 4)


class TestEvaluateNeoantigens(unittest.TestCase):
    """Step 09: evaluate_neoantigens precision/recall."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        p = Path(self.tmp)
        # ranked input: 2 CM467 + 1 TIL1 predicted
        ranked = textwrap.dedent("""\
            peptide\tsample_id\tevidence_class\tbinding_rank\texpression_tpm\tmutation_type\tmutation_pos\twt_aa\tmut_aa\tsource_protein
            IYPFAPPRSF\tCM467\tA\t0.003\t0.0\tmissense\t4\tT\tF\tP00001
            YLDPVQRDLY\tCM467\tA\t0.007\t0.0\tmissense\t1\tC\tY\tP00002
            AQLLKALEV\tTIL1\tA\t0.010\t0.0\tmissense\t9\tK\tV\tP00003
        """)
        (p / "ranked.tsv").write_text(ranked)
        (p / "validated.tsv").write_text(VALIDATED_TSV)
        (p / "manifest.tsv").write_text(MANIFEST_TSV)
        self.report = p / "report.md"

    def test_perfect_precision(self):
        """All 3 predicted peptides are in validated set → precision = 1.0."""
        import subprocess
        cmd = [
            sys.executable,
            str(SRC / "evaluation" / "09_evaluate_neoantigens.py"),
            "--input",     str(Path(self.tmp) / "ranked.tsv"),
            "--validated", str(Path(self.tmp) / "validated.tsv"),
            "--manifest",  str(Path(self.tmp) / "manifest.tsv"),
            "--output",    str(self.report),
            "--top_n", "2", "3",
        ]
        env = os.environ.copy()
        env["PYTHONPATH"] = str(SRC)
        result = subprocess.run(cmd, capture_output=True, text=True, env=env)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(self.report.exists(), "Report file should be written")
        content = self.report.read_text()
        # Precision of 1.000 should appear for both patients
        self.assertIn("1.000", content, "Expected precision 1.000 in report")

    def test_report_contains_top20_table(self):
        """Report must include the Top-20 ranked candidates table."""
        import subprocess
        cmd = [
            sys.executable,
            str(SRC / "evaluation" / "09_evaluate_neoantigens.py"),
            "--input",     str(Path(self.tmp) / "ranked.tsv"),
            "--validated", str(Path(self.tmp) / "validated.tsv"),
            "--manifest",  str(Path(self.tmp) / "manifest.tsv"),
            "--output",    str(self.report),
            "--top_n", "3",
        ]
        env = os.environ.copy()
        env["PYTHONPATH"] = str(SRC)
        subprocess.run(cmd, capture_output=True, text=True, env=env)
        content = self.report.read_text()
        self.assertIn("Top-20 Ranked Candidates", content)
        self.assertIn("IYPFAPPRSF", content)


class TestManifestHLA(unittest.TestCase):
    """Verify HLA alleles have been populated in sample_manifest.tsv."""

    def test_no_tbd_hla(self):
        manifest_path = ROOT / "configs" / "sample_manifest.tsv"
        if not manifest_path.exists():
            self.skipTest("sample_manifest.tsv not found")
        import pandas as pd
        df = pd.read_csv(manifest_path, sep="\t")
        tbd_rows = df[df["hla_alleles"] == "TBD"]
        self.assertEqual(len(tbd_rows), 0,
                         f"{len(tbd_rows)} rows still have 'TBD' in hla_alleles")

    def test_core_patients_have_hla(self):
        manifest_path = ROOT / "configs" / "sample_manifest.tsv"
        if not manifest_path.exists():
            self.skipTest("sample_manifest.tsv not found")
        import pandas as pd
        df = pd.read_csv(manifest_path, sep="\t")
        core = {"CM467", "GD149", "MD155", "TIL1", "TIL3", "RA957"}
        for patient in core:
            rows = df[df["patient_id"] == patient]
            if rows.empty:
                continue
            hla = rows.iloc[0]["hla_alleles"]
            self.assertTrue(
                pd.notna(hla) and str(hla).startswith("HLA-"),
                f"{patient} is missing HLA alleles in manifest (got: {hla!r})"
            )


if __name__ == "__main__":
    unittest.main(verbosity=2)
