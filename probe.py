#!/usr/bin/env python

"""
probe.py

User-facing entrypoint for BLASTing unassigned reads (sequenceUa).

This script:
- Parses command-line arguments.
- Calls run_probe_pipeline() to:
  - Combine all sequenceUa FASTQs into one FASTA.
  - Run BLAST+ against a specified database.
  - Filter and summarize BLAST hits.
  - Write matched_sequences.tsv and summary tables.
"""

from pathlib import Path
import argparse

from rna_pipeline import blast_runner, blast_parser, utils
from rna_pipeline.cli_common import build_probe_argparser


def run_probe_pipeline(args: argparse.Namespace) -> None:
    """
    Orchestrate the BLAST pipeline on unassigned reads.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments from build_probe_argparser().
    """
    outdir = Path(args.outdir)
    utils.ensure_dir(outdir)

    # 1) Build combined FASTA from sequenceUa FASTQs
    fasta_path = blast_runner.build_unassigned_fasta(
        input_pattern=args.input_sequences,
        outdir=outdir,
        sample_size=args.sample_size
    )

    # 2) Run BLAST
    blast_tab = blast_runner.run_blast(
        fasta_path=fasta_path,
        db=args.blast_db,
        threads=args.threads,
        outdir=outdir,
    )

    # 3) Parse, filter, and summarize hits
    blast_parser.filter_and_summarize(
        blast_tab=blast_tab,
        min_pident=args.min_pident,
        min_qcov=args.min_qcov,
        max_evalue=args.max_evalue,
        out_dir=outdir,
    )


def main() -> None:
    parser = build_probe_argparser()
    args = parser.parse_args()
    run_probe_pipeline(args)


if __name__ == "__main__":
    main()
