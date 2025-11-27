"""
blast_parser.py

Main module for parsing BLAST output files and filtering results.
"""

import pandas as pd
from pathlib import Path
from typing import Tuple

# BLAST standard 6 columns + custom columns we requested (qcovs, stitle)
BLAST_COLUMNS = [
    'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 
    'qcovs', 'stitle'
]

def getSampleFromReadId(readId):
    """
    Extracts the Sample ID from the Read ID.
    Assumes the format: 'SampleName_ReadID'
    """
    # Split by the first underscore only
    if "_" in readId:
        return readId.split("_", 1)[0]
    return "Unknown"

def filterAndSummarize(
    blastTab: Path,
    minPident: float,
    minQcov: float,
    maxEvalue: float,
    outDir: Path,
) -> Tuple[Path, Path]:
    """
    Filter BLAST hits and write matchedSequences.tsv and summaryPerSample.tsv.
    """

    # 1. Setup Output Directory
    outputFolder = Path(outDir)
    outputFolder.mkdir(parents=True, exist_ok=True)

    matchOut = outputFolder / "matchedSequences.tsv"
    summaryOut = outputFolder / "summaryPerSample.tsv"

    # 2. Check if input empty
    # If the file doesn't exist or is empty, make empty outputs and stop.
    if not blastTab.exists() or blastTab.stat().st_size == 0:
        print("Warning: {} is empty. Creating empty output files.".format(blastTab))
        pd.DataFrame(columns=BLAST_COLS).to_csv(matchOut, sep='\t', index=False)
        pd.DataFrame(columns=['sampleId', 'stitle', 'count']).to_csv(summaryOut, sep='\t', index=False)
        return matchOut, summaryOut
    
    # 3. Load Data
    print("Loading BLAST results from {}...".format(blastTab))
    try:
        # We use header=None because BLAST output usually has no header row
        df = pd.read_csv(blastTab, sep='\t', names=BLAST_COLUMNS, comment='#')
    except pd.errors.EmptyDataError:
        # Handle files that exist but have no data rows
        pd.DataFrame(columns=BLAST_COLS).to_csv(matchOut, sep='\t', index=False)
        pd.DataFrame(columns=['sampleId', 'stitle', 'count']).to_csv(summaryOut, sep='\t', index=False)
        return matchOut, summaryOut
    
    # 4. Filter the Hits 
    print("Filtering: Identity >= {}%, Coverage >= {}%, E-value <= {}".format(minPident, minQcov, maxEvalue))
    
    # Create a filter mask (True/False for each row)
    goodHitsMask = (
        (df['pident'] >= minPident) &
        (df['qcovs'] >= minQcov) & 
        (df['evalue'] <= maxEvalue)
    )

    # Keep only the rows that match our criteria
    cleanDf = df[goodHitsMask].copy()

    # 5. Find best hit per read
    # Sort by Bitscore (Highest first) and E-value (Lowest first)
    cleanDf = cleanDf.sort_values(by=["bitscore", "evalue"], ascending=[False, True])

    # Drop duplicates on 'qseqid', keeping the first one (which is the best one)
    bestHitsDf = cleanDf.drop_duplicates(subset='qseqid', keep='first')

    # 6. Extract Sample IDs
    # We use the helper function to pull the sample name from the read ID
    if not bestHitsDf.empty:
        bestHitsDf['sampleId'] = bestHitsDf['qseqid'].apply(getSampleFromReadId)
    else:
        bestHitsDf['sampleId'] = []

    # 7. Save the Detailed Report 
    print("Saving detailed matches to {}...".format(matchOut))
    bestHitsDf.to_csv(matchOut, sep='\t', index=False)

    # 8. Create and Save Summary
    # Group by Sample ID and Organism Name (stitle), then count the rows
    if not bestHitsDf.empty:
        summaryDf = bestHitsDf.groupby(['sampleId', 'stitle']).size().reset_index(name='count')
        # Sort so the most abundant organisms are at the top
        summaryDf = summaryDf.sort_values(by=['sampleId', 'count'], ascending=[True, False])
    else:
        summaryDf = pd.DataFrame(columns=['sampleId', 'stitle', 'count'])

    print("Saving summary report to {}...".format(summaryOut))
    summaryDf.to_csv(summaryOut, sep='\t', index=False)

    return matchOut, summaryOut