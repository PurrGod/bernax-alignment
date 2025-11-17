"""
star_runner.py

Wrapper around the STAR aligner.

Responsibilities:
- Build and run STAR commands for each sample.
- Collect paths to BAM and unmapped FASTQs.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Any

from .samplesheet import Sample
from . import utils


@dataclass
class StarSampleOutputs:
    """Paths to STAR outputs for a single sample."""
    bam: Path
    unmapped_fastq1: Path
    unmapped_fastq2: Path | None = None


@dataclass
class StarBatchOutputs:
    """STAR outputs for all samples in the run."""
    per_sample: Dict[str, StarSampleOutputs]

    @property
    def bam_files(self) -> List[Path]:
        """Return a list of all BAM files."""
        return [o.bam for o in self.per_sample.values()]


def run_star(sample: Sample, args: Any, ref_cfg: dict, outdir: Path) -> StarSampleOutputs:
    """
    Run STAR for a single sample.

    Parameters
    ----------
    sample : Sample
        Sample to align.
    args : argparse.Namespace
        Command-line arguments.
    ref_cfg : dict
        Reference configuration (includes STAR genome index).
    outdir : Path
        Output directory for STAR results for this sample.

    Returns
    -------
    StarSampleOutputs
    """
    # TODO: build STAR command and call utils.run_cmd
    raise NotImplementedError


def run_star_batch(samples: List[Sample], args: Any, ref_cfg: dict) -> StarBatchOutputs:
    """
    Run STAR alignment for all samples.

    Parameters
    ----------
    samples : list of Sample
    args : argparse.Namespace
    ref_cfg : dict

    Returns
    -------
    StarBatchOutputs
    """
    # TODO: implement loop over samples, call run_star
    raise NotImplementedError
