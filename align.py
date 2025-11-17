"""
align.py

User-facing entrypoint for the RNA alignment + quantification + DESeq2 pipeline.

This script:
- Parses command-line arguments.
- Calls run_align_pipeline() to:
  - Load sample sheet and reference config.
  - Run QC and trimming (optional).
  - Run STAR alignments.
  - Run featureCounts.
  - Split reads into sequenceA (assigned) and sequenceUa (unassigned).
  - Optionally run DESeq2 via R.

Intended audience: lab members with basic terminal familiarity.
"""

