"""
blast_parser.py 
This module parses BLAST output files and extracts relevant information.
"""

import os

# Created indices first for the columns for easier access and understanding

IDX_QSEQID = 0
IDX_SSEQID = 1
IDX_PIDENT = 2
IDX_LENGTH = 3
IDX_MISMATCH = 4
IDX_GAPOPEN = 5
IDX_QSTART = 6
IDX_QEND = 7
IDX_SSTART = 8
IDX_SEND = 9
IDX_EVALUE = 10
IDX_BITSCORE = 11
IDX_QCOVS = 12
IDX_STITLE = 13


def getSampleFromReadID(readID):
    """Extracts sample ID (everything before the first _)"""
    if "_" in readID:
        return readID.split("_", 1)[0]
    return "Unknown Sample"

def filterAndSummarize(blastTab, minPident, minQcov, maxEvalue, outDir):
    # Setup output files 
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    
    matchOut = os.path.join(outDir, "matchedSequences.tsv")
    summaryOut = os.path.join(outDir, "summaryPerSample.tsv")

    # bestHits dictionary to store qseqid -> (line, bitscore) 

    bestHits = {}

    print("Reading {}...".format(blastTab))

    # Read the file line-by-line
    with open(blastTab, 'r') as infile: 
        for line in infile: 
            if line.startswith('#') or not line.strip():
                continue 

            # Split the line into list of columns
            columns = line.strip().split('\t')

            # Check if line meets enough columns
            if len(columns) < 14:
                continue

            # Data extraction 
            pident = float(columns[IDX_PIDENT])
            qcovs = float(columns[IDX_QCOVS])
            evalue = float(columns[IDX_EVALUE])
            bitscore = float(columns[IDX_BITSCORE])
            qseqid = columns[IDX_QSEQID]

            # Applying the filters 
            if pident >= minPident and qcovs >= minQcov and evalue <= maxEvalue:
                # Checks if current iteration is the best hit for the readID
                if qseqid not in bestHits:
                    # First time -> save it 
                    bestHits[qseqid] = (line.strip(), bitscore)
                else:
                    # Compare existing bitscore 
                    currentBest = bestHits[qseqid][1] 
                    if bitscore > currentBest:
                        # Higher score is kept 
                        bestHits[qseqid] = (line.strip(), bitscore)

    if not bestHits:
        print("No hits passed the filtering crtiera.")

        # Write empty output tsv file with header
        with open(matchOut, 'w') as outputFile:
            header = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                      'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
                      'qcovs', 'stitle']
            outputFile.write('\t'.join(header) + '\n')
            
        with open(summaryOut, 'w') as outputFile:
            outputFile.write("sampleID\tstitle\tcount\n")

        return matchOut, summaryOut 

    # Write best hits to the output tsv file
    print("Saving matched sequences to {}...".format(matchOut))
    with open(matchOut, 'w') as outputFile: 
        # Header
        header = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
                  'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 
                  'qcovs', 'stitle']
        outputFile.write('\t'.join(header) + '\n')
        
        for hit in bestHits.values():
            outputFile.write(hit[0] + '\n')

    # Generate Summary Counts per sample using dictionary
    summaryCounts = {} # Key: (sampleID, stitle), Value: Count in int

    print("Generating summary per sample...")
    # Loop through best matches again
    for hit in bestHits.values():
        lineStr = hit[0]
        columns = lineStr.split('\t')

        # Get sampleID 
        sID = getSampleFromReadID(columns[IDX_QSEQID])
        stitle = columns[IDX_STITLE]

        # Add to summary counts
        key = (sID, stitle) 
        if key not in summaryCounts:
            summaryCounts[key] = 0
        summaryCounts[key] += 1

    # Writing Summary to output file
    print("Saving summary to {}...".format(summaryOut))
    with open(summaryOut, 'w') as outputFile:
        outputFile.write("sampleID\tstitle\tcount\n")

        # Sort alphabetically by sampleID and descending count
        # Dictionary converted to list for sorting
        sortedSummary = sorted(summaryCounts.items(), key=lambda item: (item[0][0], -item[1]))

        for (sampleId, stitle), count in sortedSummary:
            outputFile.write("{}\t{}\t{}\n".format(sampleId, stitle, count))

    return matchOut, summaryOut
