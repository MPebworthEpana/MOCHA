#!/usr/bin/env Rscript
# Build pkgdown site for MOCHA (Bioconductor Docker / CI).
args <- commandArgs(trailingOnly = TRUE)
clean <- !("--no-clean" %in% args)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

imports <- c(
  "data.table", "plyranges", "dplyr", "GenomicRanges", "RaggedExperiment",
  "MultiAssayExperiment", "SummarizedExperiment", "stringr", "ggbio", "wCorr",
  "magrittr", "rlang", "BiocGenerics", "GenomeInfoDb", "GenomicFeatures",
  "IRanges", "S4Vectors", "assertthat", "ensembldb", "ggplot2", "ggrepel",
  "matrixStats", "qvalue", "scales", "tidyr", "ggridges", "pbapply", "BSgenome",
  "tidyselect", "lifecycle", "BiocStyle", "knitr", "rmarkdown"
)
BiocManager::install(imports, ask = FALSE, update = FALSE)

if (!requireNamespace("pkgdown", quietly = TRUE)) {
  BiocManager::install("pkgdown", ask = FALSE, update = FALSE)
}

if (!dir.exists("docs") || clean) {
  if (dir.exists("docs")) {
    unlink("docs", recursive = TRUE)
  }
  dir.create("docs", showWarnings = FALSE)
}

status <- system2("R", c("CMD", "INSTALL", "."))
if (status != 0) {
  stop("MOCHA install failed (exit ", status, ").")
}

pkgdown::build_site_github_pages(
  clean = clean,
  new_process = FALSE,
  install = FALSE
)
