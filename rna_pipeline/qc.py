"""
qc.py

Quality control and trimming for FASTQ files.

This module runs:
- FastQC (before/after trimming)
- Optional trimming (cutadapt or fastp)
"""

from typing import List, Any, Dict
from pathlib import Path

from .samplesheet import Sample


def run_fastqc(samples: List[Sample], args: Any, ref_cfg: Dict) -> None:
    """
    Run FastQC on input FASTQs.

    Parameters
    ----------
    samples : list of Sample
        Samples to process.
    args : argparse.Namespace
        Command-line arguments.
    ref_cfg : dict
        Reference configuration (not heavily used here, but included for symmetry).
    """
    # TODO: implement FastQC calls (likely using utils.run_cmd)
    raise NotImplementedError


def run_trimming(samples: List[Sample], args: Any, ref_cfg: Dict) -> List[Sample]:
    """
    Run trimming (cutadapt/fastp) and return updated samples.

    Parameters
    ----------
    samples : list of Sample
        Original samples.
    args : argparse.Namespace
        Command-line arguments (e.g. trimming options).
    ref_cfg : dict
        Reference configuration.

    Returns
    -------
    List[Sample]
        New Sample objects pointing to trimmed FASTQs.
    """
    # TODO: implement
    raise NotImplementedError