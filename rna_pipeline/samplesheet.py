"""
samplesheet.py

Purpose: Parse and validate the sample sheet TSV that is provided by the user
Output: A structured Python object that contains the sample ID, FASTQ paths, and metadata columns
-This is important because it tells other modules which samples to run.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional
import csv


@dataclass
class Sample:
    """
    Represents one entry from the samplesheet, a single RNA-seq.
    It stores attributes for the Sample ID, condition, and paths to its fastq files. 
    Keeps things bundled in a small object.
    """
    id: str
    condition: str
    fastq1: Path
    fastq2: Optional[Path] = None
    # These are the specific attributes (ID, condition, paths)


def parse_samplesheet(path: Path) -> List[Sample]:
    """
    Function: parse_samplesheet
    Purpose: Reads the samplesheet TSV and converts each row into a Sample object.
    - To do this, it validates required columns, handles header capitalization, cleans
      and verifies FASTQ paths and interprets the paths
    """
    
    #Creates an empty list where samples will be put into
    samples: List[Sample] = []
    path = Path(path)

    #Makes sure the file exists first, if not raises an error
    if not path.exists():
        raise FileNotFoundError(f"Samplesheet not found: {path}")
    
    #Opens the file for reading, if there's no header it raises an error
    #Uses a csv.DictReader to check to read rows as dictionaries
    with path.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Empty samplesheet")

        #Allow case-insensitive header names by mapping lowercase -> actual name
        #This standardizes the process and makes it more forgiving
        header_map = {h.strip().lower(): h for h in reader.fieldnames}

        #Enforces the 3 key pieces of information without which the pipeline will not run properly
        for req in ("sample_id", "condition", "fastq1"):
            if req not in header_map:
                raise ValueError(f"Missing required column '{req}' in samplesheet")

        #This loop will make each line of the TSV into a Sample object
        for row in reader:
            
            #Skips fully empty rows, without breaking the script
            if not any((val or "").strip() for val in row.values()):
                continue
            
            #Unpacks values with the header map created earlier, makes the strings cleaner
            sid = (row[header_map["sample_id"]] or "").strip()
            cond = (row[header_map["condition"]] or "").strip()
            fq1_s = (row[header_map["fastq1"]] or "").strip()

            #This just checks for presence of the values, and raises an error if they are missing
            if not sid:
                raise ValueError("Found a row with empty sample_id")
            if not cond:
                raise ValueError(f"Missing condition for sample '{sid}'")
            if not fq1_s:
                raise ValueError(f"Missing fastq1 for sample '{sid}'")

            #Helper that interprets the FastQ paths
            def _to_path(pstr: str) -> Path:
                p = Path(pstr)
                
                # interpret relative paths relative to the samplesheet location
                if not p.is_absolute():
                    p = (path.parent / p).resolve(strict=False)
                return p
            
            #Used for converting FastQ1 and the optional FastQ2
            #If FastQ2 is present it treats it as a paired-end
            fq1 = _to_path(fq1_s)
            fq2 = None
            if "fastq2" in header_map:
                fq2_s = (row[header_map["fastq2"]] or "").strip()
                if fq2_s:
                    fq2 = _to_path(fq2_s)

            #This just collects the cleaned input, that will be returned by the function
            samples.append(Sample(id=sid, condition=cond, fastq1=fq1, fastq2=fq2))
    
    return samples


def validate_samples(samples: List[Sample]) -> None:
    """
    Function: Validates the data after parsing 
    Purpose: Check all Sample IDs are unique, and that FastQ files actually exists on the disk
    - This function basically just checks if the actual data makes sense, and raises one collected
        error if something is off.
    """
    
    #Creates a list that will collect all posiblle errors, and tracks IDs it already saw
    errors: List[str] = []
    seen_ids = set()

    #Loops through the objects created earlier
    for s in samples:
       
        #This checks unique IDs, a samplesheet shouldn't contain same sample twice
        if s.id in seen_ids:
            errors.append(f"Duplicate sample_id '{s.id}'")
        else:
            seen_ids.add(s.id)

        #This checks that the actual FastQ1 file exists
        if not s.fastq1.exists():
            errors.append(f"fastq1 for sample '{s.id}' not found: {s.fastq1}")

        #For the samples that have FastQ2, it checks they actually exist
        if s.fastq2 is not None and not s.fastq2.exists():
            errors.append(f"fastq2 for sample '{s.id}' not found: {s.fastq2}")

    #If there were any errors/problems discovered, it raises an error with a list of everything thats wrong
    if errors:
        raise ValueError("Sample validation failed:\n" + "\n".join(errors))