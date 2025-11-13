# bernax-alignment
bme160 project desc
.

📘 README.md — RNA-Probe Project

Version: 0.1
Audience: Wet-lab biologists, undergraduate researchers, graduate students, rotation students, and PIs with minimal computational background.

🔬 Project Summary

The RNA-Probe Project is a reproducible, Python-based pipeline that mirrors the Galaxy ref-based RNA-seq workflow used in transcriptomics research. It aligns RNA sequencing reads, counts gene expression, performs DESeq2 analysis, and then identifies unknown or unassigned reads using BLAST.

This repository provides two high-level scripts:

📂 Repository Layout
rna_probe_project/
├── align.py                 # Run complete RNA-seq → DESeq2 pipeline
├── probe.py                 # BLAST the unassigned reads
├── README.md
├── docs/
│   ├── install_guide.md
│   ├── usage_align.md
│   ├── usage_probe.md
│   └── pipeline_overview.md
├── config/
│   ├── reference_mouse.yml
│   ├── reference_human.yml
│   └── samples_example.tsv
├── envs/
│   ├── alignment.yml        # STAR, samtools, featureCounts, FastQC, BLAST+
│   └── r_deseq2.yml         # R + DESeq2 + annotation libraries
├── r_scripts/
│   └── run_deseq2.R
└── rna_pipeline/
    ├── samplesheet.py
    ├── qc.py
    ├── star_runner.py
    ├── featurecounts.py
    ├── assign_split.py
    ├── deseq2_wrapper.py
    ├── blast_runner.py
    ├── blast_parser.py
    ├── cli_common.py
    └── utils.py

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
🧬 Support

This repository is designed for use inside a research lab. If something breaks:

Ask a computationally-inclined member of the lab

Open a GitHub issue

Send an email/slack to the pipeline maintainer
