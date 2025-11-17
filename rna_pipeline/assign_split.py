"""
assign_split.py

Create sequenceA (assigned reads) and sequenceUa (unassigned reads) FASTQ files
based on featureCounts assignments and STAR outputs.
"""

from pathlib import Path
from typing import Dict, Any

from .star_runner import StarBatchOutputs
from .featurecounts import SampleAssignments
from . import utils


def build_sequence_fastqs(
    star_outputs: StarBatchOutputs,
    assignments: Dict[str, SampleAssignments],
    outdir: Path,
    args: Any,
    ref_cfg: dict,
) -> None:
    """
    Top-level function to build sequenceA and sequenceUa FASTQs for each sample.

    Parameters
    ----------
    star_outputs : StarBatchOutputs
    assignments : dict
        Mapping sample_id -> SampleAssignments.
    outdir : Path
        Output directory for split FASTQs and summary.
    args : argparse.Namespace
    ref_cfg : dict
    """
    utils.ensure_dir(outdir)
    # TODO: iterate over samples, call per-sample helpers, then write summary
    raise NotImplementedError


def build_sequenceA_for_sample(
    sample_id: str,
    sample_outputs,
    sample_assignments: SampleAssignments,
    outdir: Path,
) -> Path:
    """
    Create sequenceA FASTQ (assigned reads) for a single sample.

    Returns
    -------
    Path
        Path to sequenceA FASTQ.
    """
    # TODO: implement
    raise NotImplementedError


def build_sequenceUa_for_sample(
    sample_id: str,
    sample_outputs,
    sample_assignments: SampleAssignments,
    outdir: Path,
) -> Path:
    """
    Create sequenceUa FASTQ (unassigned reads) for a single sample.

    Returns
    -------
    Path
        Path to sequenceUa FASTQ.
    """
    # TODO: implement
    raise NotImplementedError


def write_assignment_summary(
    assignments: Dict[str, SampleAssignments],
    outpath: Path,
) -> None:
    """
    Write a summary table of assignment categories per sample.

    Parameters
    ----------
    assignments : dict
    outpath : Path
    """
    # TODO: implement
    raise NotImplementedError
