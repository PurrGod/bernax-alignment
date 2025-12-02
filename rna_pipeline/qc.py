"""
qc.py

Quality control (FastQC) and read trimming wrappers used by the RNA-seq pipeline.

These functions are intentionally lightweight stubs so that the pipeline can run
end-to-end even if no real QC / trimming tools are installed. You can later
replace the internals with calls to FastQC, fastp, trimmomatic, etc., without
changing the public interface.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Any, Dict
from argparse import Namespace

from .samplesheet import Sample
from . import utils

# Type alias for reference configuration dictionaries used across the pipeline.
RefConfig = Dict[str, Any]

logger = utils.get_logger(__name__)


def run_fastqc(samples: List[Sample], args: Namespace, ref_cfg: RefConfig) -> None:
    """
    Run FastQC on all samples.

    Parameters
    ----------
    samples : list[Sample]
        Parsed sample metadata for this run.
    args : argparse.Namespace
        Parsed command-line arguments. Must contain ``outdir`` at minimum.
    ref_cfg : RefConfig
        Reference configuration dictionary (currently unused by this stub).

    Notes
    -----
    This implementation is a no-op placeholder. It only ensures that a ``qc/``
    subdirectory exists under the output directory so that downstream steps and
    users have a visible location for QC results.

    Replace the body of this function with real FastQC invocation(s) if desired.
    """
    qc_dir = utils.subdir(Path(args.outdir), "qc")
    logger.info("QC stub: created/verified QC directory at %s", qc_dir)
    logger.debug("QC stub: FastQC not executed; %d samples provided", len(samples))
    # TODO: implement real FastQC logic here.
    return None


def run_trimming(
    samples: List[Sample],
    args: Namespace,
    ref_cfg: RefConfig,
) -> List[Sample]:
    """
    Run read trimming on all samples.

    Parameters
    ----------
    samples : list[Sample]
        Parsed sample metadata for this run.
    args : argparse.Namespace
        Parsed command-line arguments for this pipeline run.
    ref_cfg : RefConfig
        Reference configuration dictionary (currently unused by this stub).

    Returns
    -------
    list[Sample]
        The list of samples to be passed to downstream steps. For this stub,
        the input list is returned unchanged.

    Notes
    -----
    This is a no-op placeholder so the rest of the pipeline can run without
    requiring a trimming tool. To enable real trimming, update this function to
    call a tool such as fastp or trimmomatic, and adjust the returned
    ``Sample.fastq1`` / ``Sample.fastq2`` paths to point at trimmed FASTQs.
    """
    logger.info("Trimming stub: no reads were trimmed; passing through samples.")
    logger.debug("Trimming stub: %d samples provided", len(samples))
    # TODO: implement real trimming and update FASTQ paths in `samples`.
    return samples

