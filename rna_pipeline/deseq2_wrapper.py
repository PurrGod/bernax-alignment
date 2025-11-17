"""
deseq2_wrapper.py

Thin wrapper for running DESeq2 via an R script.
"""

from pathlib import Path
from typing import Dict

from . import utils


def run_deseq2(
    counts_file: Path,
    samplesheet_path: Path,
    organism: str,
    outdir: Path,
    ref_cfg: Dict,
) -> None:
    """
    Run DESeq2 and annotation in R using run_deseq2.R.

    Parameters
    ----------
    counts_file : Path
        featureCounts counts.txt.
    samplesheet_path : Path
        Samplesheet TSV file.
    organism : str
        Organism key (e.g. 'mus_musculus').
    outdir : Path
        Output directory for DESeq2 results.
    ref_cfg : dict
        Reference configuration (for organism-specific settings if needed).
    """
    utils.ensure_dir(outdir)
    # TODO: build Rscript command and call utils.run_cmd
    raise NotImplementedError
