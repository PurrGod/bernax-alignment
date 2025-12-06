# probe.py
# Group Members: Moe Sithu Maung Maung Lay, Akhilesh Nidamanuri, David Jiricek, Evan Fitzhugh


'''
probe.py

User-facing entrypoint for BLASTing unassigned reads (sequenceUa).

This script:
- Parses command-line arguments.
- Calls run_probe_pipeline() to:
  - Combine all sequenceUa FASTQs into one FASTA.
  - Run BLAST+ against a specified database.
  - Filter and summarize BLAST hits.
  - Write matched_sequences.tsv and summary tables.
'''


from pathlib import Path
import argparse

from rna_pipeline import blast_runner, blast_parser, utils
from rna_pipeline.cli_common import buildProbeArgparser


def runProbePipeline (args):
    '''
    Orchestrate the BLAST pipeline on unassigned reads.
    Inputs: args (Namespace)
    Outputs: None
    '''
    outDir = Path(args.outdir)
    utils.ensureDir(outDir)

    # 1) Build combined FASTA from sequenceUa FASTQs
    fastaPath = blast_runner.buildUnassignedFasta(
        inputPattern=args.input_sequences,
        outDir=outDir,
        sampleSize=args.sample_size
    )

    # 2) Run BLAST
    blastTab = blast_runner.runBlast(
        fastaPath=fastaPath,
        db=args.blast_db,
        threads=args.threads,
        outDir=outDir,
    )

    # 3) Parse, filter, and summarize hits
    blast_parser.filterAndSummarize(
        blastTab=blastTab,
        minPident=args.min_pident,
        minQcov=args.min_qcov,
        maxEvalue=args.max_evalue,
        outDir=outDir,
    )


def main ():
    '''
    Entry point for the script.
    '''
    parser = buildProbeArgparser()
    args = parser.parse_args()
    runProbePipeline(args)

if __name__ == "__main__":
    main()