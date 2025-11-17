"""
samplesheet.py

Parsing and validation of the samplesheet TSV into Sample objects.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional


@dataclass
class Sample:
    """
    Representation of a single RNA-seq sample.

    Attributes
    ----------
    id : str
        Unique sample ID.
    condition : str
        Experimental condition (e.g. control, treated).
    fastq1 : Path
        Path to read 1 FASTQ(.gz).
    fastq2 : Optional[Path]
        Path to read 2 FASTQ(.gz) for paired-end data, or None if single-end.
    """
    id: str
    condition: str
    fastq1: Path
    fastq2: Optional[Path] = None


def parse_samplesheet(path: Path) -> List[Sample]:
    """
    Parse a tab-delimited samplesheet into a list of Sample objects.

    Expected columns: sample_id, condition, fastq1, [fastq2]

    Parameters
    ----------
    path : Path
        Path to the samplesheet TSV.

    Returns
    -------
    List[Sample]
        List of samples.
    """
    # TODO: implement TSV parsing
    raise NotImplementedError


def validate_samples(samples: List[Sample]) -> None:
    """
    Validate sample list:
    - Ensure sample IDs are unique.
    - Ensure all FASTQ files exist.

    Parameters
    ----------
    samples : List[Sample]
        Samples to validate.

    Raises
    ------
    ValueError
        If validation fails.
    """
    # TODO: implement validation
    raise NotImplementedError