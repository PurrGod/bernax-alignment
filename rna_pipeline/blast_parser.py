"""
blast_parser.py 
This module parses BLAST output files and extracts relevant information.
"""

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

def get_sample_from_read_id(read_id):
    """Extracts sample ID (everything before the first _)"""
    if "_" in read_id:
        return read_id.split("_", 1)[0]
    return "Unknown Sample"

def filter_and_summarize(blast_tab, min_pident, min_qcov, max_evalue, out_dir):
    # Setup output files 
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    match_out = os.path.join(out_dir, "matchedSequences.tsv")
    summary_out = os.path.join(out_dir, "summaryPerSample.tsv")

    best_hits = {}

    print(f"Reading {blast_tab}...")

    # Read the file line-by-line
    with open(blast_tab, 'r') as infile: 
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
            if pident >= min_pident and qcovs >= min_qcov and evalue <= max_evalue:
                if qseqid not in best_hits:
                    best_hits[qseqid] = (line.strip(), bitscore)
                else:
                    current_best = best_hits[qseqid][1] 
                    if bitscore > current_best:
                        best_hits[qseqid] = (line.strip(), bitscore)

    if not best_hits:
        print("No hits passed the filtering criteria.")
        # Create empty files so pipeline doesn't crash later
        with open(match_out, 'w') as f:
            f.write("qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcovs\nstitle\n")
        with open(summary_out, 'w') as f:
            f.write("sampleID\tstitle\tcount\n")
        return match_out, summary_out 

    # Write best hits
    print(f"Saving matched sequences to {match_out}...")
    with open(match_out, 'w') as f: 
        f.write("qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcovs\nstitle\n")
        for hit in best_hits.values():
            f.write(hit[0] + '\n')

    # Generate Summary
    summary_counts = {} 
    print("Generating summary per sample...")
    
    for hit in best_hits.values():
        line_str = hit[0]
        columns = line_str.split('\t')
        s_id = get_sample_from_read_id(columns[IDX_QSEQID])
        stitle = columns[IDX_STITLE]

        key = (s_id, stitle) 
        if key not in summary_counts:
            summary_counts[key] = 0
        summary_counts[key] += 1

    print(f"Saving summary to {summary_out}...")
    with open(summary_out, 'w') as f:
        f.write("sampleID\tstitle\tcount\n")
        sorted_summary = sorted(summary_counts.items(), key=lambda item: (item[0][0], -item[1]))
        for (sample_id, stitle), count in sorted_summary:
            f.write(f"{sample_id}\t{stitle}\t{count}\n")

    return match_out, summary_out