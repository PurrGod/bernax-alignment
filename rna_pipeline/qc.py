# qc.py

"""
QC and trimming wrappers.

For now, these are minimal stubs so the pipeline runs without doing real QC.
You can later replace these with real FastQC / trimming steps.
"""

from typing import List, Any
from pathlib import Path

from .samplesheet import Sample
from . import utils


def run_fastqc(samples: List[Sample], args: Any, ref_cfg: dict) -> None:
    """
    Run FastQC (stub).

    Currently a no-op; implement real FastQC here later.
    """
    # Example placeholder: create a qc/ directory so users see something.
    qc_dir = utils.subdir(Path(args.outdir), "qc")
    # TODO: implement real FastQC if desired.
    return None


def run_trimming(samples: List[Sample], args: Any, ref_cfg: dict) -> List[Sample]:
    """
    Run read trimming (stub).

    Currently just returns the input samples unchanged.
    """
    # TODO: implement trimming with e.g. fastp, trimmomatic, etc.
    return samples
