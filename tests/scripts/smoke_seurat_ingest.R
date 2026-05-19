#!/usr/bin/env Rscript
# Smoke test: Seurat/Signac -> seuratToMOCHAInputs -> callOpenTiles(Seurat)
# Run from repo root with mocha-test activated:
#   NOT_CRAN=true Rscript tests/scripts/smoke_seurat_ingest.R

suppressPackageStartupMessages({
  if (requireNamespace("ggbio", quietly = TRUE)) library(ggbio)
  if (requireNamespace("ensembldb", quietly = TRUE)) library(ensembldb)
})

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools is required")
}
devtools::load_all()

helpers_path <- file.path("tests", "scripts", "smoke_seurat_helpers.R")
if (!file.exists(helpers_path)) {
  stop("Run from MOCHA repo root; helper not found: ", helpers_path)
}
source(helpers_path, local = TRUE)

if (!requireNamespace("Seurat", quietly = TRUE) ||
    !requireNamespace("Signac", quietly = TRUE)) {
  stop("Install Seurat and Signac (see tests/scripts/install-seurat-test-deps.R)")
}

rds_path <- Sys.getenv("MOCHA_SEURAT_RDS", unset = "")
if (nzchar(rds_path) && file.exists(rds_path)) {
  message("Loading Seurat object from MOCHA_SEURAT_RDS: ", rds_path)
  seurat_obj <- readRDS(rds_path)
  cell_pop_label <- Sys.getenv("MOCHA_CELL_POP_LABEL", "predicted.id")
} else {
  message("Building minimal test Seurat object ...")
  helper_path <- file.path("tests", "testthat", "helper-seurat_fragments.R")
  if (!file.exists(helper_path)) {
    stop("Run from MOCHA repo root; helper not found: ", helper_path)
  }
  source(helper_path, local = TRUE)
  meta <- data.frame(
    Sample = c("sample1", "sample1", "sample2"),
    cellPop = c("TypeA", "TypeA", "TypeA"),
    row.names = c("c1", "c2", "c3"),
    stringsAsFactors = FALSE
  )
  lines <- c(
    "chr1\t760101\t760110\tAAAC\t1",
    "chr1\t760111\t760120\tAAAG\t1",
    "chr1\t760121\t760130\tTTTC\t1"
  )
  seurat_obj <- make_test_seurat_chromatin(fragment_lines = lines, meta = meta)
  if (is.null(seurat_obj)) {
    stop("Failed to build test Seurat object")
  }
  cell_pop_label <- "cellPop"
}

message("Running seuratToMOCHAInputs ...")
inputs <- MOCHA::seuratToMOCHAInputs(
  seuratObj = seurat_obj,
  cellPopLabel = cell_pop_label,
  verbose = TRUE
)

stopifnot(is.list(inputs))
stopifnot("ATACFragments" %in% names(inputs))
stopifnot("cellColData" %in% names(inputs))
stopifnot(methods::is(inputs$ATACFragments, "GRangesList"))
stopifnot(length(inputs$ATACFragments) >= 1L)
stopifnot("RG" %in% colnames(GenomicRanges::mcols(inputs$ATACFragments[[1]])))
message("seuratToMOCHAInputs OK: ", length(inputs$ATACFragments), " sample(s)")

missing_pkgs <- smoke_callOpenTiles_missing_pkgs()
if (length(missing_pkgs) > 0) {
  stop(
    "callOpenTiles smoke requires: ",
    paste(missing_pkgs, collapse = ", "),
    ". Install via conda (environment-mocha-test.yml) or BiocManager."
  )
}

blacklist <- smoke_load_blacklist()
cell_pops <- unique(as.character(inputs$cellColData[[cell_pop_label]]))
out_dir <- tempfile("MOCHA_smoke_")

message("Running callOpenTiles(Seurat) ...")
message("  outDir: ", out_dir)
tiles <- MOCHA::callOpenTiles(
  ATACFragments = seurat_obj,
  blackList = blacklist,
  genome = "BSgenome.Hsapiens.UCSC.hg19",
  cellColData = inputs$cellColData,
  cellPopLabel = cell_pop_label,
  cellPopulations = cell_pops,
  studySignal = 100,
  TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
  OrgDb = "org.Hs.eg.db",
  outDir = out_dir,
  numCores = 1,
  verbose = TRUE
)

stopifnot(methods::is(tiles, "MultiAssayExperiment"))
stopifnot(length(tiles) >= 1L)
message(
  "callOpenTiles(Seurat) OK: ",
  length(tiles),
  " cell population assay(s); names: ",
  paste(names(tiles), collapse = ", ")
)
message("Smoke test passed.")
