# assign_split.py

"""
Create sequenceA (assigned reads) and sequenceUa (unassigned reads) FASTQ files
based on featureCounts assignments and STAR outputs.
"""

from pathlib import Path
from typing import Dict, Any, Optional

from .star_runner import StarBatchOutputs, StarSampleOutputs
from .featurecounts import SampleAssignments
from . import utils
import gzip


def _get_input_fastq_from_sample_outputs(sample_outputs: StarSampleOutputs | None) -> Optional[Path]:
    """
    Locate the unmapped FASTQ for a sample from StarSampleOutputs.

    Prefer unmapped_fastq1; fall back to unmapped_fastq2; return None if missing.
    """
    if sample_outputs is None:
        return None

    # Primary: unmapped mate 1
    if getattr(sample_outputs, "unmapped_fastq1", None):
        return Path(sample_outputs.unmapped_fastq1)

    # Fallback: unmapped mate 2
    if getattr(sample_outputs, "unmapped_fastq2", None):
        return Path(sample_outputs.unmapped_fastq2)

    return None


def _open_maybe_gz(path: Path):
    """Open text-mode for plain or gzip FASTQ files."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "r")


def _extract_read_id_from_header(header: str) -> str:
    """
    Given a FASTQ header line (starting with '@'), return the read id token
    that's typically used to match assignments.
    """
    header = header.strip()
    if header.startswith("@"):
        header = header[1:]
    # Read id is usually the first whitespace-separated token
    return header.split()[0]


def _get_assigned_sets(sample_assignments: SampleAssignments):
    """
    Extract assigned/unassigned read id sets and category counts
    from the SampleAssignments dataclass.
    """
    if sample_assignments is None:
        return None, None, None

    assigned = sample_assignments.assigned_ids
    unassigned = sample_assignments.unassigned_ids
    counts = sample_assignments.category_counts
    return assigned, unassigned, counts


def _write_filtered_fastq(
    in_fastq: Path,
    out_fastq: Path,
    keep_ids: Optional[set],
) -> int:
    """
    Stream through in_fastq and write records whose read-id is in keep_ids.
    If keep_ids is None, write no records and return 0.
    Returns the number of records written.
    """
    written = 0
    if keep_ids is None:
        # nothing to write
        open(out_fastq, "w").close()
        return 0

    try:
        with _open_maybe_gz(in_fastq) as inh, open(out_fastq, "w") as outh:
            while True:
                h = inh.readline()
                if not h:
                    break
                s = inh.readline()
                p = inh.readline()
                q = inh.readline()
                if not q:
                    break  # malformed but stop
                rid = _extract_read_id_from_header(h)
                if rid in keep_ids:
                    outh.write(h)
                    outh.write(s)
                    outh.write(p)
                    outh.write(q)
                    written += 1
    except FileNotFoundError:
        # input missing -> create empty output
        open(out_fastq, "w").close()
        written = 0

    return written


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

    for sample_id, sample_assignments in assignments.items():
        # find corresponding sample outputs in star_outputs
        sample_outputs = star_outputs.per_sample.get(sample_id, None)

        seqA_path = build_sequenceA_for_sample(sample_id, sample_outputs, sample_assignments, outdir)
        seqUa_path = build_sequenceUa_for_sample(sample_id, sample_outputs, sample_assignments, outdir)

        # optional logging
        assigned_set, unassigned_set, _ = _get_assigned_sets(sample_assignments)
        try:
            if hasattr(utils, "log"):
                utils.log(
                    f"{sample_id}: sequenceA={seqA_path} "
                    f"(n={len(assigned_set) if assigned_set is not None else 'NA'}), "
                    f"sequenceUa={seqUa_path} "
                    f"(n={len(unassigned_set) if unassigned_set is not None else 'NA'})"
                )
        except Exception:
            pass

    # write a simple summary TSV
    summary_path = outdir / "assignment_summary.tsv"
    write_assignment_summary(assignments, summary_path)


def build_sequenceA_for_sample(
    sample_id: str,
    sample_outputs: StarSampleOutputs | None,
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
    utils.ensure_dir(outdir)
    assigned_set, _, _ = _get_assigned_sets(sample_assignments)

    outpath = outdir / f"{sample_id}.sequenceA.fastq"
    in_fastq = _get_input_fastq_from_sample_outputs(sample_outputs)

    if in_fastq is None:
        # no input FASTQ found; create empty output
        open(outpath, "w").close()
        return outpath

    written = _write_filtered_fastq(in_fastq, outpath, assigned_set)
    try:
        if hasattr(utils, "log"):
            utils.log(f"Wrote {written} assigned reads for sample {sample_id} -> {outpath}")
    except Exception:
        pass

    return outpath


def build_sequenceUa_for_sample(
    sample_id: str,
    sample_outputs: StarSampleOutputs | None,
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
    utils.ensure_dir(outdir)
    _, unassigned_set, _ = _get_assigned_sets(sample_assignments)

    outpath = outdir / f"{sample_id}.sequenceUa.fastq"
    in_fastq = _get_input_fastq_from_sample_outputs(sample_outputs)

    if in_fastq is None:
        open(outpath, "w").close()
        return outpath

    written = _write_filtered_fastq(in_fastq, outpath, unassigned_set)
    try:
        if hasattr(utils, "log"):
            utils.log(f"Wrote {written} unassigned reads for sample {sample_id} -> {outpath}")
    except Exception:
        pass

    return outpath


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
    all_categories = set()
    sample_stats: Dict[str, Dict[str, int]] = {}

    for sample_id, samp in assignments.items():
        assigned_set, unassigned_set, counts = _get_assigned_sets(samp)

        stats: Dict[str, int] = {}

        if assigned_set is not None:
            stats["Assigned"] = len(assigned_set)
        if unassigned_set is not None:
            stats["Unassigned"] = len(unassigned_set)

        if counts:
            # Use human-readable category names
            for cat, v in counts.items():
                stats[cat.name] = v

        sample_stats[sample_id] = stats
        all_categories.update(stats.keys())

    # write TSV
    with open(outpath, "w") as outfh:
        header = ["sample_id"] + sorted(all_categories)
        outfh.write("\t".join(header) + "\n")
        for sample_id in sorted(sample_stats.keys()):
            row = [sample_id]
            stats = sample_stats[sample_id]
            for c in sorted(all_categories):
                v = stats.get(c)
                row.append(str(v) if v is not None else "NA")
            outfh.write("\t".join(row) + "\n")

