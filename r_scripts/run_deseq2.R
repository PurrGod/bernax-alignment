#!/usr/bin/env Rscript

'''run_deseq2.R

Run DESeq2 on a featureCounts counts matrix and annotate results.

Intended to mirror the Galaxy ref-based RNA-seq tutorial behavior.'''

suppressPackageStartupMessages({
  library("optparse")
  library("DESeq2")
  library("readr")
  library("dplyr")
  # TODO: add organism-specific annotation libraries (e.g. org.Mm.eg.db)
})

option_list <- list(
  make_option(c("--counts"), type = "character", help = "featureCounts counts.txt"),
  make_option(c("--samples"), type = "character", help = "Samplesheet TSV"),
  make_option(c("--organism"), type = "character", help = "Organism key, e.g. mus_musculus"),
  make_option(c("--outdir"), type = "character", help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$counts) || is.null(opt$samples) || is.null(opt$outdir)) {
  stop("Missing required arguments: --counts, --samples, --outdir")
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# TODO:
# 1) Read counts table
# 2) Read samplesheet and build colData
# 3) Create DESeqDataSet, run DESeq()
# 4) Extract results, write full and filtered tables
# 5) Add annotation columns
# 6) Generate basic plots (MA, PCA)

quit(save = "no")


