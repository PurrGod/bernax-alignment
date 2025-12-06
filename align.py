# align.py
# Group Members: Moe Sithu Maung Maung Lay, Akhilesh Nidamanuri, David Jiricek, Evan Fitzhugh

'''
align.py

User-facing entrypoint for the RNA alignment + quantification pipeline.

This script:
- Parses command-line arguments.
- Calls runAlignPipeline() to:
  - Load sample sheet and reference config.
  - Run STAR alignments.
  - Run featureCounts.

Intended audience: lab members with basic terminal familiarity.
'''

from pathlib import Path
import argparse

from rna_pipeline import (
    samplesheet,
    qc, # optional
    star_runner,
    featurecounts,
    deseq2_wrapper,
    utils,
)
from rna_pipeline.cli_common import buildAlignArgparser


def runAlignPipeline (args):
    '''
    Main orchestrator for the alignment pipeline.
    Inputs: args (Namespace)
    Outputs: None
    '''
    outDir = Path(args.outdir)
    utils.ensureDir(outDir)

    # Load reference configuration (STAR index, GTF, etc.)
    configPath = Path(args.reference_config) if args.reference_config else None
    refCfg = utils.loadReferenceConfig(configPath)

    # Parse and validate samplesheet
    sampleList = samplesheet.parseSamplesheet(Path(args.samples))
    samplesheet.validateSamples(sampleList)

    # QC and trimming -> optional
    if not args.skip_qc:
        qc.runFastqc(sampleList, args, refCfg)
    if args.trim:
        sampleList = qc.runTrimming(sampleList, args, refCfg)

    # STAR alignment
    starOutputs = star_runner.runStarBatch(sampleList, args, refCfg)

    # featureCounts quantification
    fcResult = featurecounts.runFeatureCounts(
        starOutputs=starOutputs,
        args=args,
        refCfg=refCfg,
        outDir=outDir / "featureCounts",
    )

    # Parsing assignments (Useful for debugging, but splitting handled by STAR)
    assignments = featurecounts.parseAssignments(fcResult)

    # Optional DESeq2
    if args.run_deseq2:
        deseq2_wrapper.run_deseq2(
            counts_file=fcResult.countsFile,
            samplesheet_path=Path(args.samples),
            organism=args.organism,
            outdir=outDir / "deseq2",
            ref_cfg=refCfg,
        )

def main ():
    '''
    Entry point for the script.
    '''
    parser = buildAlignArgparser()
    args = parser.parse_args()
    runAlignPipeline(args)

if __name__ == "__main__":
    main()
