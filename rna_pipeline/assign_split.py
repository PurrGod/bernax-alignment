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
import gzip
from typing import Optional, Iterable # Iterable actually isn't used but might be useful later

#


def _get_input_fastq_from_sample_outputs(sample_outputs) -> Optional[Path]:
    """
    Try a variety of common attribute/key names to locate the original/unmapped FASTQ
    for a sample in the STAR outputs object. Return a Path or None.
    """
    candidates = [
        "fastq", "fastq1", "reads", "input_fastq", "input_reads",
        "unmapped_fastq", "unmapped_out_mate1", "Unmapped.out.mate1",
        "Unmapped", "unmapped",
    ]

    # sample_outputs might be a dict-like or an object with attributes
    for name in candidates:
        # dict-like access
        try:
            if isinstance(sample_outputs, dict) and name in sample_outputs:
                p = sample_outputs[name]
                if p:
                    return Path(p)
        except Exception:
            pass

        # attribute access
        try:
            if hasattr(sample_outputs, name):
                p = getattr(sample_outputs, name)
                if p:
                    return Path(p)
        except Exception:
            pass

    # if sample_outputs itself looks like a path
    try:
        if isinstance(sample_outputs, (str, Path)):
            return Path(sample_outputs)
    except Exception:
        pass

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


def _get_assigned_sets(sample_assignments) -> (Optional[set], Optional[set], Optional[dict]):
    """
    Try to extract assigned/unassigned read name sets and an optional counts/categories mapping
    from the SampleAssignments object. Returns (assigned_set, unassigned_set, counts_dict)
    where any element may be None if not available.
    """
    assigned = None
    unassigned = None
    counts = None

    # If sample_assignments is a dict-like structure
    try:
        if isinstance(sample_assignments, dict):
            # common keys
            if "assigned" in sample_assignments:
                assigned = set(sample_assignments["assigned"])
            if "assigned_reads" in sample_assignments:
                assigned = set(sample_assignments["assigned_reads"])
            if "unassigned" in sample_assignments:
                unassigned = set(sample_assignments["unassigned"])
            if "unassigned_reads" in sample_assignments:
                unassigned = set(sample_assignments["unassigned_reads"])
            if "assignments" in sample_assignments and isinstance(sample_assignments["assignments"], dict):
                # assignments: read -> category
                mapping = sample_assignments["assignments"]
                assigned = set(k for k, v in mapping.items() if v and v != "unassigned")
                unassigned = set(k for k, v in mapping.items() if not v or v == "unassigned")
            if "counts" in sample_assignments and isinstance(sample_assignments["counts"], dict):
                counts = dict(sample_assignments["counts"])
    except Exception:
        pass

    # Try attribute-based access
    try:
        if hasattr(sample_assignments, "assigned"):
            val = getattr(sample_assignments, "assigned")
            if val is not None:
                assigned = set(val)
    except Exception:
        pass

    try:
        if hasattr(sample_assignments, "assigned_reads"):
            val = getattr(sample_assignments, "assigned_reads")
            if val is not None:
                assigned = set(val)
    except Exception:
        pass

    try:
        if hasattr(sample_assignments, "unassigned"):
            val = getattr(sample_assignments, "unassigned")
            if val is not None:
                unassigned = set(val)
    except Exception:
        pass

    try:
        if hasattr(sample_assignments, "unassigned_reads"):
            val = getattr(sample_assignments, "unassigned_reads")
            if val is not None:
                unassigned = set(val)
    except Exception:
        pass

    try:
        if hasattr(sample_assignments, "assignments") and isinstance(getattr(sample_assignments, "assignments"), dict):
            mapping = getattr(sample_assignments, "assignments")
            assigned = set(k for k, v in mapping.items() if v and v != "unassigned")
            unassigned = set(k for k, v in mapping.items() if not v or v == "unassigned")
    except Exception:
        pass

    try:
        if hasattr(sample_assignments, "counts"):
            val = getattr(sample_assignments, "counts")
            if isinstance(val, dict):
                counts = dict(val)
    except Exception:
        pass

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

    results = {}
    for sample_id, sample_assignments in assignments.items():
        # find corresponding sample outputs in star_outputs
        sample_outputs = None
        try:
            if isinstance(star_outputs, dict) and sample_id in star_outputs:
                sample_outputs = star_outputs[sample_id]
            elif hasattr(star_outputs, "outputs") and sample_id in getattr(star_outputs, "outputs"):
                sample_outputs = getattr(star_outputs, "outputs")[sample_id]
            elif hasattr(star_outputs, sample_id):
                sample_outputs = getattr(star_outputs, sample_id)
            elif hasattr(star_outputs, "get"):
                sample_outputs = star_outputs.get(sample_id)
        except Exception:
            sample_outputs = None

        seqA_path = build_sequenceA_for_sample(sample_id, sample_outputs, sample_assignments, outdir)
        seqUa_path = build_sequenceUa_for_sample(sample_id, sample_outputs, sample_assignments, outdir)

        # collect counts if possible
        assigned_set, unassigned_set, counts = _get_assigned_sets(sample_assignments)
        stats = {
            "sequenceA": seqA_path,
            "sequenceUa": seqUa_path,
            "assigned_count": len(assigned_set) if assigned_set is not None else None,
            "unassigned_count": len(unassigned_set) if unassigned_set is not None else None,
        }
        if counts:
            stats["counts"] = counts
        results[sample_id] = stats

    # write a simple summary TSV
    summary_path = outdir / "assignment_summary.tsv"
    write_assignment_summary(assignments, summary_path)

    return None


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
    utils.ensure_dir(outdir)
    assigned_set, _, _ = _get_assigned_sets(sample_assignments)

    outpath = outdir / f"{sample_id}.sequenceA.fastq"
    in_fastq = _get_input_fastq_from_sample_outputs(sample_outputs)

    if in_fastq is None:
        # no input FASTQ found; create empty output
        open(outpath, "w").close()
        return outpath

    written = _write_filtered_fastq(in_fastq, outpath, assigned_set)
    # optional: annotate via utils if available
    try:
        if hasattr(utils, "log"):
            utils.log(f"Wrote {written} assigned reads for sample {sample_id} -> {outpath}")
    except Exception:
        pass

    return outpath


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
    # prepare header columns
    # try to detect category names from any counts dict present
    all_categories = set()
    sample_stats = {}

    for sample_id, samp in assignments.items():
        _, _, counts = _get_assigned_sets(samp)
        if counts:
            sample_stats[sample_id] = dict(counts)
            all_categories.update(counts.keys())
        else:
            # fallback to Assigned/Unassigned counts if available
            assigned_set, unassigned_set, _ = _get_assigned_sets(samp)
            stats = {}
            if assigned_set is not None:
                stats["Assigned"] = len(assigned_set)
            if unassigned_set is not None:
                stats["Unassigned"] = len(unassigned_set)
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
