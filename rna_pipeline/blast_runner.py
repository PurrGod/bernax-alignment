"""
blast_runner.py

Build combined FASTA from sequenceUa fastqs and run BLAST+.
"""

from pathlib import Path
from typing import Any

from . import utils


def build_unassigned_fasta(input_pattern: str, outdir: Path) -> Path:
    """
    Convert all sequenceUa FASTQs matching a pattern into a single FASTA file.

    Parameters
    ----------
    input_pattern : str
        Glob pattern for sequenceUa FASTQs.
    outdir : Path

    Returns
    -------
    Path
        Path to combined FASTA file.
     """
    
    outdir.mkdir(parents=True, exist_ok=True)
    #^^Checks if the output directory exists,if not creates it, if it does it skips creation and continues.

    combined_fasta = outdir / "sequenceUa_combined.fasta"
    #^^Name of the output combined FASTA file, that will be created in here.

    fastq_files = sorted(Path().glob(input_pattern))
    #^^Find all the fastq files matching the input pattern and sort them., no specific names, just patterns.

    with combined_fasta.open('w') as fasta_out:
    #^^Opens the output FASTA file for writing.

        for fq in fastq_files:
        #^^Starts a loop and goes thru each fastq file found
           
           with fq.open('r') as fastq_in:
            #^^Opens the current fastq file for reading
               while True:
                   header = fastq_in.readline()
                   if not header:
                       break  
                   seq = fastq_in.readline()
                   plus = fastq_in.readline()
                   qual = fastq_in.readline()
                #^^This loop should read 4 lines at a time from the fastq file
                #^^Unless it is corrupted, it should follow, header, sequence, plus, and quality for each line
                #^^ Once it reaches the end of the file it break the loop

                   header = header.strip()
                   seq = seq.strip()
                    #^^ Since we only need the header and sequence for FASTA, we strip them
                    
                   if header.startswith('@'):
                       fasta_header = '>' + header[1:]
                    #^^ Should change the fastq header into the fasta format
                   else:
                       fasta_header = '>' + header
                    #^^ If for some reason it doesn't start with @, this adds the > anyway
                    
                   fasta_out.write(f"{fasta_header}\n{seq}\n")
                   #^^ Writes the fasta with only the header and sequence into the output file
    
    return combined_fasta
    #^^ Finally, returns the path to combined FASTA file


def run_blast(
    fasta_path: Path,
    db: str,
    threads: int,
    outdir: Path,
) -> Path:
    """
    Run BLAST+ on the combined unassigned FASTA.

    Parameters
    ----------
    fasta_path : Path
    db : str
        BLAST database path or prefix.
    threads : int
    outdir : Path

    Returns
    -------
    Path
        Path to raw BLAST tabular output file.
    """
    
    outdir.mkdir(parents=True, exist_ok=True)
    #^^ Same as before, makes sure the output directory exists
    
    blast_out = outdir / f"{fasta_path.stem}.blast.tsv"
    #^^ Name the blast output file according to the input fasta file name

    outfmt = ("6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs stitle")
    #^^ Blast output format, with the columns described in blast_parser.py

    cmd = [
        "blastn",                       
        "-query", str(fasta_path),
        "-db", db,
        "-out", str(blast_out),
        "-evalue", "1e-5",
        "-num_threads", str(threads),
        "-outfmt", outfmt
    ]
    #^^ This is the list of command line arguments that can be used to run blast in the terminal
    
    utils.run_cmd(cmd)
    #^^ This is the thing that actually runs the command

    return blast_out
    #^^ Finally returns the path to the blast output file
    
