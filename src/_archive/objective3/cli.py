"""Command-line entrypoint for the integrated Objective 3 workflow."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .pipeline import (
    DEFAULT_BINDING_RANK_SUPPORT_MAX,
    DEFAULT_DE_NOVO_SCORE_CUTOFF,
    DEFAULT_EXPRESSION_TPM_SUPPORT_MIN,
    DEFAULT_PSM_FDR_THRESHOLD,
    dataset_status,
    normalize_maxquant_psms,
    postprocess,
    predict_de_novo,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Integrated Objective 3 helper CLI")
    sub = parser.add_subparsers(dest="command", required=True)

    status = sub.add_parser("status", help="Summarize MGF and MaxQuant result coverage")
    status.add_argument("--mgf-dir", default="/home/amity/hla_data_mgf")
    status.add_argument("--results-dir", default="production_results")

    norm = sub.add_parser("normalize-maxquant-psms", help="Normalize MaxQuant msms files")
    norm.add_argument("--results-dir", default="production_results")
    norm.add_argument("--output", default="results/objective3/immunopeptidome_psms.tsv")
    norm.add_argument("--fdr-column")

    predict = sub.add_parser("predict-denovo", help="Run the custom CNN/LSTM model on MGF spectra")
    predict.add_argument("--mgf-dir", default="/home/amity/hla_data_mgf")
    predict.add_argument("--checkpoint", default="neoepitope_model_epoch_50.pth")
    predict.add_argument("--output", default="results/objective3/de_novo_predictions.tsv")
    predict.add_argument("--max-spectra-per-file", type=int)
    predict.add_argument("--min-length", type=int, default=8)
    predict.add_argument("--max-length", type=int, default=11)

    post = sub.add_parser("postprocess", help="Subtract database hits and rank de novo candidates")
    post.add_argument("--psms", default="results/objective3/immunopeptidome_psms.tsv")
    post.add_argument("--denovo", default="results/objective3/de_novo_predictions.tsv")
    post.add_argument("--outdir", default="results/objective3/final")
    post.add_argument("--sample-metadata")
    post.add_argument("--binding-input")
    post.add_argument("--fasta")
    post.add_argument("--psm-fdr-threshold", type=float, default=DEFAULT_PSM_FDR_THRESHOLD)
    post.add_argument("--de-novo-score-cutoff", type=float, default=DEFAULT_DE_NOVO_SCORE_CUTOFF)
    post.add_argument("--binding-rank-support-max", type=float, default=DEFAULT_BINDING_RANK_SUPPORT_MAX)
    post.add_argument("--expression-tpm-support-min", type=float, default=DEFAULT_EXPRESSION_TPM_SUPPORT_MIN)

    run = sub.add_parser("run-objective3", help="Run normalize, optional prediction, and postprocess")
    run.add_argument("--mgf-dir", default="/home/amity/hla_data_mgf")
    run.add_argument("--results-dir", default="production_results")
    run.add_argument("--outdir", default="results/objective3")
    run.add_argument("--checkpoint", default="neoepitope_model_epoch_50.pth")
    run.add_argument("--skip-predict", action="store_true")
    run.add_argument("--existing-denovo")
    run.add_argument("--max-spectra-per-file", type=int)
    run.add_argument("--sample-metadata")
    run.add_argument("--binding-input")
    run.add_argument("--fasta")

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "status":
        summary = dataset_status(Path(args.mgf_dir), Path(args.results_dir))
        print(json.dumps(summary, indent=2))
        return 0

    if args.command == "normalize-maxquant-psms":
        rows = normalize_maxquant_psms(Path(args.results_dir), Path(args.output), fdr_column=args.fdr_column)
        print(f"Wrote {len(rows)} PSM rows to {args.output}")
        return 0

    if args.command == "predict-denovo":
        rows = predict_de_novo(
            mgf_dir=Path(args.mgf_dir),
            checkpoint=Path(args.checkpoint),
            output=Path(args.output),
            max_spectra_per_file=args.max_spectra_per_file,
            min_length=args.min_length,
            max_length=args.max_length,
        )
        print(f"Wrote {len(rows)} de novo rows to {args.output}")
        return 0

    if args.command == "postprocess":
        outputs = postprocess(
            psm_path=Path(args.psms),
            de_novo_path=Path(args.denovo),
            outdir=Path(args.outdir),
            metadata_path=Path(args.sample_metadata) if args.sample_metadata else None,
            binding_input=Path(args.binding_input) if args.binding_input else None,
            fasta_path=Path(args.fasta) if args.fasta else None,
            psm_fdr_threshold=args.psm_fdr_threshold,
            de_novo_score_cutoff=args.de_novo_score_cutoff,
            binding_rank_support_max=args.binding_rank_support_max,
            expression_tpm_support_min=args.expression_tpm_support_min,
        )
        print("\n".join(str(path) for path in outputs))
        return 0

    if args.command == "run-objective3":
        outdir = Path(args.outdir)
        psm_path = outdir / "immunopeptidome_psms.tsv"
        prediction_path = Path(args.existing_denovo) if args.existing_denovo else outdir / "de_novo_predictions.tsv"
        normalize_maxquant_psms(Path(args.results_dir), psm_path)
        if not args.skip_predict and not args.existing_denovo:
            predict_de_novo(
                mgf_dir=Path(args.mgf_dir),
                checkpoint=Path(args.checkpoint),
                output=prediction_path,
                max_spectra_per_file=args.max_spectra_per_file,
            )
        if not prediction_path.exists():
            parser.error("No de novo prediction file available. Provide --existing-denovo or omit --skip-predict.")
        outputs = postprocess(
            psm_path=psm_path,
            de_novo_path=prediction_path,
            outdir=outdir / "final",
            metadata_path=Path(args.sample_metadata) if args.sample_metadata else None,
            binding_input=Path(args.binding_input) if args.binding_input else None,
            fasta_path=Path(args.fasta) if args.fasta else None,
        )
        print("\n".join(str(path) for path in outputs))
        return 0

    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
