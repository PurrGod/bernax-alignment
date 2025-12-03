"""
cli_common.py

Shared command-line argument parsers and helpers for align.py and probe.py.
"""

import argparse


def build_align_argparser() -> argparse.ArgumentParser:
    """
    Build an ArgumentParser for align.py.

    Returns
    -------
    argparse.ArgumentParser
        Configured parser for the alignment pipeline.
    """
    p = argparse.ArgumentParser(
        description="Run RNA-seq alignment + quantification + DESeq2 pipeline."
    )
    p.add_argument("--samples", required=True, help="Path to samplesheet TSV.")
    p.add_argument("--sample-size", type=int, default=10000, help="Number of reads to subsample for BLAST.")
    p.add_argument("--genome-index", required=True, help="STAR genome index directory.")
    p.add_argument("--gtf", required=True, help="Gene annotation GTF file.")
    p.add_argument("--outdir", required=True, help="Output directory.")
    p.add_argument("--organism", required=True, help="Organism key (e.g. mus_musculus).")
    p.add_argument("--reference-config", required=False, default=None,
                   help="YAML file with reference configuration (STAR index, GTF, BLAST DB).")
    p.add_argument("--threads", type=int, default=4, help="Number of threads to use.")
    p.add_argument("--skip-qc", action="store_true", help="Skip FastQC.")
    p.add_argument("--trim", action="store_true", help="Enable read trimming step.")
    p.add_argument("--run-deseq2", action="store_true", help="Run DESeq2 analysis.")
    return p


def build_probe_argparser() -> argparse.ArgumentParser:
    """
    Build an ArgumentParser for probe.py.

    Returns
    -------
    argparse.ArgumentParser
        Configured parser for the BLAST pipeline.
    """
    p = argparse.ArgumentParser(
        description="Run BLAST on unassigned reads (sequenceUa)."
    )
    p.add_argument(
        "--input-sequences",
        required=True,
        help='Glob pattern or comma-separated list of sequenceUa FASTQs, e.g. "results/assign_split/sequenceUa_*.fastq".',
    )
    p.add_argument("--blast-db", required=True, help="BLAST+ database path/prefix.")
    p.add_argument("--outdir", required=True, help="Output directory for BLAST results.")
    p.add_argument("--threads", type=int, default=4, help="Number of BLAST threads.")
    p.add_argument("--min-pident", type=float, default=90.0, help="Minimum percent identity.")
    p.add_argument("--min-qcov", type=float, default=0.7, help="Minimum query coverage (fraction).")
    p.add_argument("--max-evalue", type=float, default=1e-5, help="Maximum E-value.")
    p.add_argument("--sample-size", type=int, default=10000, help="Number of reads to subsample for BLAST.")
    return p