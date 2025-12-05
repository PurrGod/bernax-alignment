"""
blast_runner.py

Build combined FASTA from sequenceUa fastqs (Subsampled) and run BLAST+.
"""

from pathlib import Path
import gzip
from . import utils

def _open_maybe_gz(path: Path):
    """Open text-mode for plain or gzip FASTQ files."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "r")

def build_unassigned_fasta(input_pattern: str, outdir: Path, sample_size: int = 10000) -> Path:
    """
    Convert sequenceUa FASTQs into a single FASTA file.
    
    Parameters
    ----------
    input_pattern : str
        Glob pattern for input files.
    outdir : Path
        Output directory.
    sample_size : int
        Max number of reads to process (Subsampling).
    """
    outdir.mkdir(parents=True, exist_ok=True)
    combined_fasta = outdir / "sequenceUa_combined.fasta"
    fastq_files = sorted(Path().glob(input_pattern))
    
    
    count = 0

    print(f"Building BLAST input from {len(fastq_files)} files (Subsampling {sample_size} reads)...")

    with combined_fasta.open('w') as fasta_out:
        for fq in fastq_files:
            if count >= sample_size: break # Stop if global limit reached
            
            try:
                with _open_maybe_gz(fq) as fastq_in:
                    while True:
                        if count >= sample_size: break
                        
                        # FASTQ Record = 4 lines
                        header = fastq_in.readline()
                        if not header: break # End of file
                        seq = fastq_in.readline()
                        fastq_in.readline() # Plus line
                        fastq_in.readline() # Quality line
                        
                        # Convert to FASTA (Header + Sequence)
                        # Ensure header starts with '>'
                        header_clean = header.strip()
                        if header_clean.startswith("@"):
                            header_clean = ">" + header_clean[1:]
                        else:
                            header_clean = ">" + header_clean
                            
                        sample_name = fq.name.split('.')[0]
                        fasta_out.write(f">{sample_name}_{header_clean[1:]}\n{seq.strip()}\n")
                        #Use the FASTQ file name as the sample prefix so that the files are more easily traced
                        #Format according to how normal FASTA files are formatted
                        
                        count += 1
            except Exception as e:
                print(f"Warning: Could not read {fq}: {e}")
                #This ensures that a bad file does not disrupt the entire BLAST preparation step, sends out a warning

    print(f"Created {combined_fasta} with {count} sequences.")
    return combined_fasta


def run_blast(fasta_path: Path, db: str, threads: int, outdir: Path) -> Path:
    """
    Run BLAST+ on the combined unassigned FASTA.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    blast_out = outdir / f"{fasta_path.stem}.blast.tsv"
    
    # Custom output format matching blast_parser expectations
    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs stitle"

    cmd = [
        "blastn",
        "-query", str(fasta_path),
        "-db", db,
        "-out", str(blast_out),
        "-evalue", "1e-5",
        "-num_threads", str(threads),
        "-outfmt", outfmt
    ]
    #Delegate to the shared command runner so logging or error handling is consistent across the pipeline
    utils.run_cmd(cmd)
    return blast_out
    