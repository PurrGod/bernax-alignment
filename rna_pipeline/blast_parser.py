# blast_parser.py
# Group Members: Moe Sithu Maung Maung Lay, Akhilesh Nidamanuri, David Jiricek, Evan Fitzhugh


'''
blast_parser.py 
Purpose: This module parses BLAST output files and extracts relevant information.
'''

import os

# Column Indices
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

def getSampleFromReadId (readId):
    '''
    Extracts sample ID (everything before the first _).
    Inputs: readId (str)
    Outputs: Sample ID (str)
    '''
    if "_" in readId:
        return readId.split("_", 1)[0]
    return "Unknown Sample"

def filterAndSummarize (blastTab, minPident, minQcov, maxEvalue, outDir):
    '''
    Parses BLAST results, filters them, and writes summary reports.
    Inputs: blastTab (Path), minPident (float), minQcov (float), maxEvalue (float), outDir (Path)
    Outputs: Tuple of paths (matchOut, summaryOut)
    '''
    # Setup output files 
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    
    matchOut = os.path.join(outDir, "matchedSequences.tsv")
    summaryOut = os.path.join(outDir, "summaryPerSample.tsv")

    bestHits = {}

    print(f"Reading {blastTab}...")

    with open(blastTab, 'r') as infile: 
        for line in infile: 
            if line.startswith('#') or not line.strip():
                continue 

            columns = line.strip().split('\t')
            if len(columns) < 14:
                continue

            pident = float(columns[IDX_PIDENT])
            qcovs = float(columns[IDX_QCOVS])
            evalue = float(columns[IDX_EVALUE])
            bitscore = float(columns[IDX_BITSCORE])
            qseqid = columns[IDX_QSEQID]

            # Filters 
            if pident >= minPident and qcovs >= minQcov and evalue <= maxEvalue:
                # Check if current iteration is the best hit for the readID
                if qseqid not in bestHits:
                    # save it for the first hit
                    bestHits[qseqid] = (line.strip(), bitscore)
                else:
                    # Compare existing bitscore, higher score is saved
                    currentBest = bestHits[qseqid][1] 
                    if bitscore > currentBest:
                        bestHits[qseqid] = (line.strip(), bitscore)

    if not bestHits:
        print("No hits passed the filtering criteria.")

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
    print(f"Saving matched sequences to {matchOut}...")
    with open(matchOut, 'w') as outputFile: 
        header = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
                  'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 
                  'qcovs', 'stitle']
        outputFile.write('\t'.join(header) + '\n')
        
        for hit in bestHits.values():
            outputFile.write(hit[0] + '\n')

    # Generate Summary Counts per sample using dictionary
    # Key: (sampleID, stitle), Value: Count in int
    summaryCounts = {} 

    print("Generating summary per sample...")
    # Loop through best matches again
    for hit in bestHits.values():
        lineStr = hit[0]
        columns = lineStr.split('\t')

        sId = getSampleFromReadId(columns[IDX_QSEQID])
        stitle = columns[IDX_STITLE]

        key = (sId, stitle) 
        if key not in summaryCounts:
            summaryCounts[key] = 0
        summaryCounts[key] += 1

    # Writing Summary to output file
    print(f"Saving summary to {summaryOut}...")
    with open(summaryOut, 'w') as outputFile:
        outputFile.write("sampleID\tstitle\tcount\n")

        # Sort alphabetically by sampleID and descending count
        # Dictionary converted to list for sorting
        sortedSummary = sorted(summaryCounts.items(), key=lambda item: (item[0][0], -item[1]))

        for (sampleId, stitle), count in sortedSummary:
            outputFile.write(f"{sampleId}\t{stitle}\t{count}\n")

    return matchOut, summaryOut