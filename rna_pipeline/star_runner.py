# star_runner.py

"""
Wrapper around the STAR aligner.

Responsibilities:
- Build and run STAR commands for each sample.
- Collect paths to BAM and unmapped FASTQs.
"""

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
    Run STAR for a single sample.

    Parameters
    ----------
    sample : Sample
        Sample to align.
    args : argparse.Namespace
        Command-line arguments.
    ref_cfg : dict
        Reference configuration (includes STAR genome index).
    outdir : Path
        Output directory for STAR results for this sample.

    Returns
    -------
    StarSampleOutputs
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
    Run STAR alignment for all samples.

    Parameters
    ----------
    samples : list of Sample
    args : argparse.Namespace
    ref_cfg : dict

    Returns
    -------
    StarBatchOutputs
    """
    # Use a 'star' subdir under the main outdir
    base_outdir = utils.subdir(Path(args.outdir), "star")

    per_sample: Dict[str, StarSampleOutputs] = {}
    for sample in samples:
        sample_outdir = utils.subdir(base_outdir, sample.id)
        outputs = run_star(sample, args, ref_cfg, sample_outdir)
        per_sample[sample.id] = outputs

    return StarBatchOutputs(per_sample=per_sample)
