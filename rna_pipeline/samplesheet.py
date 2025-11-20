"""
samplesheet.py

Parsing and validation of the samplesheet TSV into Sample objects.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional
import csv


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
    samples: List[Sample] = []
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"Samplesheet not found: {path}")

    with path.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Empty samplesheet")

        # allow case-insensitive header names by mapping lowercase -> actual name
        header_map = {h.strip().lower(): h for h in reader.fieldnames}

        for req in ("sample_id", "condition", "fastq1"):
            if req not in header_map:
                raise ValueError(f"Missing required column '{req}' in samplesheet")

        for row in reader:
            # skip fully empty rows
            if not any((val or "").strip() for val in row.values()):
                continue

            sid = (row[header_map["sample_id"]] or "").strip()
            cond = (row[header_map["condition"]] or "").strip()
            fq1_s = (row[header_map["fastq1"]] or "").strip()

            if not sid:
                raise ValueError("Found a row with empty sample_id")
            if not cond:
                raise ValueError(f"Missing condition for sample '{sid}'")
            if not fq1_s:
                raise ValueError(f"Missing fastq1 for sample '{sid}'")

            def _to_path(pstr: str) -> Path:
                p = Path(pstr)
                if not p.is_absolute():
                    # interpret relative paths relative to the samplesheet location
                    p = (path.parent / p).resolve(strict=False)
                return p

            fq1 = _to_path(fq1_s)

            fq2 = None
            if "fastq2" in header_map:
                fq2_s = (row[header_map["fastq2"]] or "").strip()
                if fq2_s:
                    fq2 = _to_path(fq2_s)

            samples.append(Sample(id=sid, condition=cond, fastq1=fq1, fastq2=fq2))

    return samples


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
    errors: List[str] = []
    seen_ids = set()

    for s in samples:
        # check unique IDs
        if s.id in seen_ids:
            errors.append(f"Duplicate sample_id '{s.id}'")
        else:
            seen_ids.add(s.id)

        # check fastq1 exists
        if not s.fastq1.exists():
            errors.append(f"fastq1 for sample '{s.id}' not found: {s.fastq1}")

        # check fastq2 exists when provided
        if s.fastq2 is not None and not s.fastq2.exists():
            errors.append(f"fastq2 for sample '{s.id}' not found: {s.fastq2}")

    if errors:
        raise ValueError("Sample validation failed:\n" + "\n".join(errors))