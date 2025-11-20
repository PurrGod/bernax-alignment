#!/usr/bin/env python

"""
align.py

User-facing entrypoint for the RNA alignment + quantification + DESeq2 pipeline.

This script:
- Parses command-line arguments.
- Calls run_align_pipeline() to:
  - Load sample sheet and reference config.
  - Run QC and trimming (optional).
  - Run STAR alignments.
  - Run featureCounts.
  - Split reads into sequenceA (assigned) and sequenceUa (unassigned).
  - Optionally run DESeq2 via R.

Intended audience: lab members with basic terminal familiarity.
"""

from pathlib import Path
import argparse

from rna_pipeline import (
    samplesheet,
    qc,
    star_runner,
    featurecounts,
    assign_split,
    deseq2_wrapper,
    utils,
)
from rna_pipeline.cli_common import build_align_argparser


def run_align_pipeline(args: argparse.Namespace) -> None:
    """
    Main orchestrator for the alignment pipeline.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments from build_align_argparser().
    """
    outdir = Path(args.outdir)
    utils.ensure_dir(outdir)

    # Load reference configuration (STAR index, GTF, etc.)
    ref_cfg = utils.load_reference_config(Path(args.reference_config))

    # Parse and validate samplesheet
    sample_list = samplesheet.parse_samplesheet(Path(args.samples))
    samplesheet.validate_samples(sample_list)

    # QC and trimming
    if not args.skip_qc:
        qc.run_fastqc(sample_list, args, ref_cfg)
    if args.trim:
        sample_list = qc.run_trimming(sample_list, args, ref_cfg)

    # STAR alignment
    star_outputs = star_runner.run_star_batch(sample_list, args, ref_cfg)

    # featureCounts quantification
    fc_result = featurecounts.run_featurecounts(
        star_outputs=star_outputs,
        args=args,
        ref_cfg=ref_cfg,
        outdir=outdir / "featureCounts",
    )

    assignments = featurecounts.parse_assignments(fc_result)

    # Split into assigned and unassigned read sets
    assign_split.build_sequence_fastqs(
        star_outputs=star_outputs,
        assignments=assignments,
        outdir=outdir / "assign_split",
        args=args,
        ref_cfg=ref_cfg,
    )

    # Optional DESeq2
    if args.run_deseq2:
        deseq2_wrapper.run_deseq2(
            counts_file=fc_result.counts_file,
            samplesheet_path=Path(args.samples),
            organism=args.organism,
            outdir=outdir / "deseq2",
            ref_cfg=ref_cfg,
        )


def main() -> None:
    parser = build_align_argparser()
    args = parser.parse_args()
    run_align_pipeline(args)


if __name__ == "__main__":
    main()

