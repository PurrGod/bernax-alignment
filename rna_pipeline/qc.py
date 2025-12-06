"""
qc.py

Purpose: This file is used for quality control (FastQC) and trimming for the pipeline
Output: Creates a directory for QC and logs messages, however it is not implemented yet, 
        so it is just a placeholder for now.
"""

#Imports all the useful tools for QC and trimming as well as the Samples from the samplesheet, created earlier
from __future__ import annotations
from pathlib import Path
from typing import List, Any, Dict
from argparse import Namespace
from .samplesheet import Sample
from . import utils

#Types alias for reference configuration dictionaries used across the pipeline.
RefConfig = Dict[str, Any]

#Sets up a logger, so log messages go into the same place as the rest of the pipeline
logger = utils.getLogger(__name__)


def runFastqc(samples: List[Sample], args: Namespace, ref_cfg: RefConfig) -> None:
    """
    Function: runFastqc
    Purpose: Runs FastQC on all samples
    - Right now it just sets up the folder, but it will run quality check on all the
        fastq1 and fastq2 data. Eventually this would output per-sample QC reports into
        the qc/ directory. For now it is only a placeholder and does not run FastQC yet.
    """

    #Takes the main output directory and makes sure "qc" subdirectory exists inside
    #The path will go towards that directory
    qc_dir = utils.subDir(Path(args.outdir), "qc")
    
    #Creates information messages about QC
    logger.info("QC stub: created/verified QC directory at %s", qc_dir)
    logger.debug("QC stub: FastQC not executed; %d samples provided", len(samples))
    # TODO: implement real FastQC logic here.
    
    #Doesnt return anything right now
    return None


def runTrimming(    
    samples: List[Sample],
    args: Namespace,
    ref_cfg: RefConfig,
) -> List[Sample]:
    """
    Function: runTrimming
    Purpose: A trimming tool for all the samples
    - Right now it returns the original samples unchanged
    """
   
    #Creates information messages for the user about trimming information
    logger.info("Trimming stub: no reads were trimmed; passing through samples.")
    logger.debug("Trimming stub: %d samples provided", len(samples))
    # TODO: implement real trimming and update FASTQ paths in `samples`.
    
    #Retuns original samples for now
    return samples

