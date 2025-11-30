# Pipeline Overview

## Problem we are investigating
Given the dataset acquired from the Chamorro-Garcia Lab has unmapped reads of >60%, we believe it can be a useful marker on these two hypotheses:
1. Purple Bar = Evidence of chromatin disruption
2. Black Bar = Evidence of contamination

## Solution we came up with 
We have split the entire project into two parts, creating two main scripts: align.py and probe.py

### align.py
- The script will align the genes the mouse genome using STAR
    - It will then split the data into two areas: Aligned and Unaligned
        - Aligned reads that hit the mouse genome
        - Unaligned reads that do not match with the mouse genome


### probe.py 
- The script will take the unaligned reads from align.py
    - The reads will be run through BLAST against a known database 
        - The output will demonstate what reads are being matched with which species 


## Summary of the flowchart
[ Raw Data ] --> [ align.py ] --(Aligned) --> [ Further Analysis ]
                       |
                   (Unmapped)
                       |
                       V
                 [ probe.py ] --> [ Further Analysis ]