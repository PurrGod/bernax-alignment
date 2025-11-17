"""
blast_parser.py

Parse BLAST tabular output, filter hits, and produce summary tables.
"""

from pathlib import Path
from typing import Tuple

import pandas as pd  # ensure added to requirements


def filter_and_summarize(
    blast_tab: Path,
    min_pident: float,
    min_qcov: float,
    max_evalue: float,
    outdir: Path,
) -> Tuple[Path, Path]:
    """
    Filter BLAST hits and write matched_sequences.tsv and summary_per_sample.tsv.

    Parameters
    ----------
    blast_tab : Path
        Path to raw BLAST tabular output.
    min_pident : float
    min_qcov : float
    max_evalue : float
    outdir : Path

    Returns
    -------
    (Path, Path)
        Paths to matched_sequences.tsv and summary_per_sample.tsv.
    """
    # TODO: implement using pandas
    raise NotImplementedError
