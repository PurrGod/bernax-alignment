# bernax-alignment

Version: 0.1

Audience: Wet-lab biologists, undergraduate researchers, graduate students, rotation students, and PIs with minimal computational background.

## Overview

The RNA-Probe Project is a reproducible, Python-based pipeline that mirrors a reference-based RNA-seq workflow. It performs the following high-level tasks:

- Align RNA sequencing reads
- Count gene expression
- Run DESeq2 differential expression analysis
- Identify unknown or unassigned reads using BLAST

The repository exposes two top-level scripts for common workflows:

- `align.py` — Run the complete RNA-seq → DESeq2 pipeline
- `probe.py` — BLAST the unassigned reads

## Repository layout

The important files and directories are:

```
bernax-alignment/
├── align.py                 # Run complete RNA-seq → DESeq2 pipeline
├── probe.py                 # BLAST the unassigned reads
├── README.md
├── docs/                    # Usage and pipeline documentation
│   ├── install_guide.md
│   ├── usage_align.md
│   ├── usage_probe.md
│   └── pipeline_overview.md
├── config/                  # Example reference and sample config files
│   ├── reference_mouse.yml
│   ├── reference_human.yml
│   └── samples_example.tsv
├── envs/                    # Conda environment YAMLs (recommended)
│   ├── alignment.yml        # STAR, samtools, featureCounts, FastQC, BLAST+
│   └── r_deseq2.yml         # R + DESeq2 + annotation libraries
├── r_scripts/               # R helper scripts for DESeq2
│   └── run_deseq2.R
└── rna_pipeline/            # Core Python modules used by the two scripts
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
```

## Documentation

Detailed usage and installation instructions are in the `docs/` folder:

- `docs/install_guide.md` — Environment setup and dependencies
- `docs/usage_align.md` — How to run `align.py`
- `docs/usage_probe.md` — How to run `probe.py`
- `docs/pipeline_overview.md` — High-level pipeline description

## Quick notes

- This project is intended to be run in a reproducible environment — using the Conda YAMLs in `envs/` is recommended.
- The pipeline expects a properly formatted samplesheet and reference YAML (see `config/` for examples).
- For DESeq2 steps, R and the necessary Bioconductor libraries are required (see `envs/r_deseq2.yml`).

## Support

If you run into problems, try the following:

1. Ask a computationally-inclined member of your lab.
2. Open a GitHub issue in this repository.
3. Contact the pipeline maintainer via email or Slack.

## License

This repository includes a `LICENSE` file at the project root. Check it for license details.

## Contributing

Contributions are welcome. Please open an issue first to discuss larger changes.

## Credits

This project was developed as part of the BME160 course materials at UC Santa Cruz (SHOUTOUT DAVID BERNICK). Code Base developed in cooperation by Evan F, Moe SLL, David J, Akhilesh N.
