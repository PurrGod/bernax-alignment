"""
blast_runner.py

Build combined FASTA from sequenceUa fastqs and run BLAST+.
"""

from pathlib import Path
from typing import Any

from . import utils


def build_unassigned_fasta(input_pattern: str, outdir: Path) -> Path:
    """
    Convert all sequenceUa FASTQs matching a pattern into a single FASTA file.

    Parameters
    ----------
    input_pattern : str
        Glob pattern for sequenceUa FASTQs.
    outdir : Path

    Returns
    -------
    Path
        Path to combined FASTA file.
    """
    # TODO: implement (e.g. using seqtk or Biopython)
    raise NotImplementedError


def run_blast(
    fasta_path: Path,
    db: str,
    threads: int,
    outdir: Path,
) -> Path:
    """
    Run BLAST+ on the combined unassigned FASTA.

    Parameters
    ----------
    fasta_path : Path
    db : str
        BLAST database path or prefix.
    threads : int
    outdir : Path

    Returns
    -------
    Path
        Path to raw BLAST tabular output file.
    """
    # TODO: build and run blastn command via utils.run_cmd
    raise NotImplementedError
