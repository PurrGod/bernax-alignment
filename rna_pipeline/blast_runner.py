# blast_runner.py
# Group Members: Moe Sithu Maung Maung Lay, Akhilesh Nidamanuri, David Jiricek, Evan Fitzhugh

'''
blast_runner.py

Purpose: Build combined FASTA from sequenceUA FastQs and run BLASTN/BLAST+ on them.
'''

from pathlib import Path
import gzip
from . import utils

def _openMaybeGz (path):
    '''
    Open text-mode for plain or gzip FASTQ files.
    Inputs: path (Path)
    Outputs: File handle
    '''
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "r")

def buildUnassignedFasta (inputPattern, outDir, sampleSize=10000):
    '''
    Convert sequenceUa FASTQs into a single FASTA file.
    Inputs: inputPattern (str), outDir (Path), sampleSize (int)
    Outputs: Path to combined FASTA file
    '''
    
    # Checks that the output directory exists, if not it creates one
    outDir.mkdir(parents=True, exist_ok=True)
    combinedFasta = outDir / "sequenceUa_combined.fasta"
    fastqFiles = sorted(Path().glob(inputPattern))
    
    count = 0

    print(f"Building BLAST input from {len(fastqFiles)} files (Subsampling {sampleSize} reads)...")

    with combinedFasta.open('w') as fastaOut:
        for fq in fastqFiles:
            if count >= sampleSize: break # Stop if global limit reached
            
            try:
                with _openMaybeGz(fq) as fastqIn:
                    while True:
                        if count >= sampleSize: break
                        
                        # FASTQ Record = 4 lines
                        header = fastqIn.readline()
                        if not header: break # End of file
                        seq = fastqIn.readline()
                        fastqIn.readline() # Plus line
                        fastqIn.readline() # Quality line
                        
                        # Convert to FASTA (Header + Sequence)
                        # Ensure header starts with '>'
                        headerClean = header.strip()
                        if headerClean.startswith("@"):
                            headerClean = ">" + headerClean[1:]
                        else:
                            headerClean = ">" + headerClean
                       
                        # Use the FASTQ file name as the sample prefix
                        sampleName = fq.name.split('.')[0]
                        
                        # Put sample name FIRST to prevent BLAST from truncating
                        fastaOut.write(f">{sampleName}_{headerClean[1:]}\n{seq.strip()}\n")
                         
                        count += 1
            
            except Exception as e:
                print(f"Warning: Could not read {fq}: {e}")

    # Prints a message with information about the created FASTA file.            
    print(f"Created {combinedFasta} with {count} sequences.")
    return combinedFasta

def runBlast (fastaPath, db, threads, outDir):
    '''
    Run BLAST+ on the combined unassigned FASTA.
    Inputs: fastaPath (Path), db (str), threads (int), outDir (Path)
    Outputs: Path to BLAST output file
    '''
    outDir.mkdir(parents=True, exist_ok=True)
    blastOut = outDir / f"{fastaPath.stem}.blast.tsv"
    
    # Custom output format matching blast_parser expectations
    outFmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs stitle"

    cmd = [
        "blastn",
        "-query", str(fastaPath),
        "-db", db,
        "-out", str(blastOut),
        "-evalue", "1e-5",
        "-num_threads", str(threads),
        "-outfmt", outFmt
    ]
    
    # Delegate to the shared command runner
    utils.runCmd(cmd)
    return blastOut