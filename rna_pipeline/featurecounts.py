# featurecounts.py

"""
Wrapper for featureCounts.

Responsibilities:
- Run featureCounts on all BAM files.
- Parse counts.txt and per-read assignment BAMs (-R BAM).
"""

from dataclasses import dataclass
from enum import Enum, auto
from pathlib import Path
from typing import Dict, Any

from .star_runner import StarBatchOutputs
from . import utils
import pysam


class AssignmentCategory(Enum):
    """Categorization of read assignment outcomes."""
    ASSIGNED = auto()
    UNASSIGNED_UNMAPPED = auto()
    UNASSIGNED_NO_FEATURES = auto()
    UNASSIGNED_MAPPING_QUALITY = auto()
    UNASSIGNED_AMBIGUITY = auto()


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
    utils.ensure_dir(outdir)

    counts_file = outdir / "counts.txt"
    summary_file = outdir / "counts.txt.summary"

    # Determine GTF from ref_cfg or args
    gtf = ref_cfg.get("gtf") or getattr(args, "gtf", None)
    if gtf is None:
        raise ValueError(
            "No GTF file provided. Use --gtf or include 'gtf' in reference-config."
        )

    # Collect BAM files from STAR
    bam_files = [o.bam for o in star_outputs.per_sample.values()]
    if not bam_files:
        raise ValueError("No BAM files found in StarBatchOutputs.")

    cmd = [
        "featureCounts",
        "-T", str(getattr(args, "threads", 4)),
        "-a", str(gtf),
        "-o", str(counts_file),
        "-R", "BAM",
        *(str(b) for b in bam_files),
    ]

    utils.run_cmd(cmd)

    # featureCounts with -R BAM creates <input>.featureCounts.bam
    per_sample_assignment_files: Dict[str, Path] = {}
    for sample_id, star_out in star_outputs.per_sample.items():
        output_bam_name = star_out.bam.name + ".featureCounts.bam"
        per_sample_assignment_files[sample_id] = outdir / output_bam_name


    return FeatureCountsResult(
        counts_file=counts_file,
        summary_file=summary_file,
        per_sample_assignment_files=per_sample_assignment_files,
    )


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
    assignments: Dict[str, SampleAssignments] = {}

    for sample_id, bam_path in fc_result.per_sample_assignment_files.items():
        assigned_ids: set[str] = set()
        unassigned_ids: set[str] = set()
        category_counts: Dict[AssignmentCategory, int] = {
            cat: 0 for cat in AssignmentCategory
        }

        if not bam_path.exists():
            raise FileNotFoundError(f"Expected featureCounts output BAM not found at {bam_path}.")

        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for read in bam:
                read_id = read.query_name
                tag = read.get_tag("XS") if read.has_tag("XS") else None

                if tag == "Assigned":
                    assigned_ids.add(read_id)
                    category_counts[AssignmentCategory.ASSIGNED] += 1
                else:
                    unassigned_ids.add(read_id)
                    if tag == "Unassigned_Unmapped":
                        category_counts[AssignmentCategory.UNASSIGNED_UNMAPPED] += 1
                    elif tag == "Unassigned_NoFeatures":
                        category_counts[AssignmentCategory.UNASSIGNED_NO_FEATURES] += 1
                    elif tag == "Unassigned_MappingQuality":
                        category_counts[AssignmentCategory.UNASSIGNED_MAPPING_QUALITY] += 1
                    elif tag == "Unassigned_Ambiguous":
                        category_counts[AssignmentCategory.UNASSIGNED_AMBIGUITY] += 1

        assignments[sample_id] = SampleAssignments(
            assigned_ids=assigned_ids,
            unassigned_ids=unassigned_ids,
            category_counts=category_counts,
        )

    return assignments

