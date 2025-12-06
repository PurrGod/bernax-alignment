# star_runner.py
# Group Members: Moe Sithu Maung Maung Lay, Akhilesh Nidamanuri, David Jiricek, Evan Fitzhugh

'''
star_runner.py

Purpose: STAR aligner wrapper, that builds and runs STAR commands for each sample. 
        It collects paths to BAM and unmapped FastQs. For each sample the module outputs the aligned BAM file 
        and the unmapped reads, and wraps them in a clean Python object. When running multiple samples, 
        it returns a dictionary mapping each sample to its output files.
'''

#Imports useful tools as well as mainly the samples from the samplesheet
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Any
from .samplesheet import Sample
from . import utils

@dataclass
class StarSampleOutputs :
    '''Paths to STAR outputs for a single sample.'''
    bam: Path
    unmappedFastq1: Path
    unmappedFastq2: Path | None = None

@dataclass
class StarBatchOutputs :
    '''STAR outputs for all samples in the run.'''
    perSample: Dict[str, StarSampleOutputs]

    @property
    def bamFiles (self) -> List[Path]:
        '''Return a list of all BAM files.'''
        return [o.bam for o in self.perSample.values()]

def _resolveGenomeIndex (args: Any, refCfg: dict) -> Path:
    '''
    Determine the STAR genome index directory from refCfg or args.
    '''
    genomeDir = (
        refCfg.get("star_index")
        or refCfg.get("genome_index")
        or getattr(args, "genome_index", None)
    )
    if genomeDir is None:
        raise ValueError(
            "No STAR genome index provided. "
            "Set --genome-index or provide 'star_index' in reference-config."
        )
    return Path(genomeDir)

def runStar (sample: Sample, args: Any, refCfg: dict, outDir: Path) -> StarSampleOutputs:
    '''
    Function: runStar
    Purpose: Run STAR for a single sample
    - Returns a StarSampleOutputs file containing the BAM path and unmapped FastQ files.
    '''
    
    utils.ensureDir(outDir)
    genomeDir = _resolveGenomeIndex(args, refCfg)

    # Input reads
    readFiles = [sample.fastq1]
    if sample.fastq2 is not None:
        readFiles.append(sample.fastq2)

    threads = getattr(args, "threads", 4)

    # STAR output prefix: outDir/sampleid_
    prefix = outDir / f"{sample.id}_"

    cmd = [
        "STAR",
        "--genomeDir", str(genomeDir),
        "--runThreadN", str(threads),
        "--readFilesIn", *(str(p) for p in readFiles),
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--outFileNamePrefix", str(prefix),
        "--outReadsUnmapped", "Fastx",
        "--twopassMode", "Basic",
        # Added to ensure stability on small genomes
        "--outFilterScoreMinOverLread", "0",
        "--outFilterMatchNminOverLread", "0",
        "--outFilterMatchNmin", "0"
    ]

    # Use zcat for gzipped FASTQ
    if any(str(p).endswith((".gz", ".gzip")) for p in readFiles):
        cmd.extend(["--readFilesCommand", "zcat"])

    utils.runCmd(cmd)

    bam = Path(str(prefix) + "Aligned.sortedByCoord.out.bam")
    unmapped1 = Path(str(prefix) + "Unmapped.out.mate1")
    unmapped2 = Path(str(prefix) + "Unmapped.out.mate2") if sample.fastq2 is not None else None

    return StarSampleOutputs(
        bam=bam,
        unmappedFastq1=unmapped1,
        unmappedFastq2=unmapped2,
    )

def runStarBatch (samples: List[Sample], args: Any, refCfg: dict) -> StarBatchOutputs:
    '''
    Function: runStarBatch
    Purpose: Run STAR alignment for all samples.
    - Loops over all samples, creates per-sample STAR output directories, 
      calls runStar for each, and returns StarBatchOutputs object mapping sample IDs to their STAR outputs.
    - StarBatchOutputs object contains the map in .perSample
    '''

    # Use a 'star' subDir under the main outdir
    baseOutDir = utils.subDir(Path(args.outdir), "star")

    perSample: Dict[str, StarSampleOutputs] = {}
    for sample in samples:
        sampleOutDir = utils.subDir(baseOutDir, sample.id)
        outputs = runStar(sample, args, refCfg, sampleOutDir)
        perSample[sample.id] = outputs

    return StarBatchOutputs(perSample=perSample)