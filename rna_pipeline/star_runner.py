"""
star_runner.py

Purpose: STAR aligner wrapper, that builds and runs STAR commands for each sample. 
        It collects paths to BAM and unmapped FastQs. For each sample the module outputs the aligned BAM file 
        and the unmapped reads, and wraps them in a clean Python object. When running multiple samples, 
        it returns a dictionary mapping each sample to its output files.
"""

#Imports useful tools as well as mainly the samples from the samplesheet
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Any
from .samplesheet import Sample
from . import utils

@dataclass
class StarSampleOutputs:
    """Paths to STAR outputs for a single sample."""
    bam: Path
    unmapped_fastq1: Path
    unmapped_fastq2: Path | None = None

@dataclass
class StarBatchOutputs:
    """STAR outputs for all samples in the run."""
    per_sample: Dict[str, StarSampleOutputs]

    @property
    def bam_files(self) -> List[Path]:
        """Return a list of all BAM files."""
        return [o.bam for o in self.per_sample.values()]

def _resolve_genome_index(args: Any, ref_cfg: dict) -> Path:
    """
    Determine the STAR genome index directory from ref_cfg or args.
    """
    genome_dir = (
        ref_cfg.get("star_index")
        or ref_cfg.get("genome_index")
        or getattr(args, "genome_index", None)
    )
    if genome_dir is None:
        raise ValueError(
            "No STAR genome index provided. "
            "Set --genome-index or provide 'star_index' in reference-config."
        )
    return Path(genome_dir)

def run_star(sample: Sample, args: Any, ref_cfg: dict, outdir: Path) -> StarSampleOutputs:
    """
    Function: run_star
    Purpose: Run STAR for a single sample
    - Returns a StarSampleOutputs file containing the BAM path and unmapped FastQ files.
    """
    
    utils.ensure_dir(outdir)
    genome_dir = _resolve_genome_index(args, ref_cfg)

    # Input reads
    read_files = [sample.fastq1]
    if sample.fastq2 is not None:
        read_files.append(sample.fastq2)

    threads = getattr(args, "threads", 4)

    # STAR output prefix: outdir/sampleid_
    prefix = outdir / f"{sample.id}_"

    cmd: list[str] = [
        "STAR",
        "--genomeDir", str(genome_dir),
        "--runThreadN", str(threads),
        "--readFilesIn", *(str(p) for p in read_files),
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--outFileNamePrefix", str(prefix),
        "--outReadsUnmapped", "Fastx",
        "--twopassMode", "Basic",
    ]

    # Use zcat for gzipped FASTQ
    if any(str(p).endswith((".gz", ".gzip")) for p in read_files):
        cmd.extend(["--readFilesCommand", "zcat"])

    utils.run_cmd(cmd)

    bam = Path(str(prefix) + "Aligned.sortedByCoord.out.bam")
    unmapped1 = Path(str(prefix) + "Unmapped.out.mate1")
    unmapped2 = Path(str(prefix) + "Unmapped.out.mate2") if sample.fastq2 is not None else None

    return StarSampleOutputs(
        bam=bam,
        unmapped_fastq1=unmapped1,
        unmapped_fastq2=unmapped2,
    )

def run_star_batch(samples: List[Sample], args: Any, ref_cfg: dict) -> StarBatchOutputs:
    """
    Function: run_star_batch
    Purpose: Run STAR alignment for all samples.
    - Loops over all samples, creates per-sample STAR output directories, 
      calls run_star for each, and returns a batch object mapping sample IDs to their STAR outputs.
    """
  
    # Use a 'star' subdir under the main outdir
    base_outdir = utils.subdir(Path(args.outdir), "star")

    per_sample: Dict[str, StarSampleOutputs] = {}
    for sample in samples:
        sample_outdir = utils.subdir(base_outdir, sample.id)
        outputs = run_star(sample, args, ref_cfg, sample_outdir)
        per_sample[sample.id] = outputs

    return StarBatchOutputs(per_sample=per_sample)
