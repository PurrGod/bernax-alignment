"""
featurecounts.py

Wrapper for featureCounts.

Responsibilities:
- Run featureCounts on all BAM files.
- Parse counts.txt and counts.txt.summary.
- Parse per-read assignment files when -R NAME is used.
"""

from dataclasses import dataclass
from enum import Enum, auto
from pathlib import Path
from typing import Dict, Any

from .star_runner import StarBatchOutputs


class AssignmentCategory(Enum):
    """Categorization of read assignment outcomes."""
    ASSIGNED = auto()
    UNASSIGNED_UNMAPPED = auto()
    UNASSIGNED_NO_FEATURES = auto()
    UNASSIGNED_MAPPING_QUALITY = auto()
    UNASSIGNED_AMBIGUITY = auto()
    # TODO: extend as needed


@dataclass
class FeatureCountsResult:
    """Paths to key featureCounts output files."""
    counts_file: Path
    summary_file: Path
    per_sample_assignment_files: Dict[str, Path]


@dataclass
class SampleAssignments:
    """Per-sample assignment information."""
    assigned_ids: set[str]
    unassigned_ids: set[str]
    category_counts: Dict[AssignmentCategory, int]


def run_featurecounts(
    star_outputs: StarBatchOutputs,
    args: Any,
    ref_cfg: dict,
    outdir: Path,
) -> FeatureCountsResult:
    """
    Run featureCounts on all BAM files.

    Parameters
    ----------
    star_outputs : StarBatchOutputs
        BAM file locations.
    args : argparse.Namespace
    ref_cfg : dict
    outdir : Path

    Returns
    -------
    FeatureCountsResult
    """
    # TODO: build featureCounts command and run via utils.run_cmd
    raise NotImplementedError


def parse_assignments(fc_result: FeatureCountsResult) -> Dict[str, SampleAssignments]:
    """
    Parse read assignment files into SampleAssignments.

    Parameters
    ----------
    fc_result : FeatureCountsResult

    Returns
    -------
    dict
        Mapping from sample ID to SampleAssignments.
    """
    # TODO: implement parsing logic
    raise NotImplementedError
