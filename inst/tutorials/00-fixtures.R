# Build shared MOCHA objects for extended tutorial scripts (Tier A/B).
# Run after loading MOCHA; sources bundled example data through the workflow.

shared_path <- file.path(
  Sys.getenv("MOCHA_TUTORIAL_DIR", unset = ""),
  "_shared.R"
)
if (!nzchar(Sys.getenv("MOCHA_TUTORIAL_DIR", unset = ""))) {
  shared_path <- system.file("tutorials", "_shared.R", package = "MOCHA")
  if (!nzchar(shared_path)) {
    shared_path <- file.path(getwd(), "inst", "tutorials", "_shared.R")
  }
}
if (file.exists(shared_path)) {
  source(shared_path, local = FALSE)
}

if (!requireNamespace("MOCHA", quietly = TRUE)) {
  stop("MOCHA must be loaded before sourcing 00-fixtures.R")
}

tutorial_check_deps()

# ---- fixtures-call-open-tiles ----
tileResults <- MOCHA::callOpenTiles(
  ATACFragments = MOCHA::exampleFragments,
  cellColData = MOCHA::exampleCellColData,
  blackList = MOCHA::exampleBlackList,
  genome = "BSgenome.Hsapiens.UCSC.hg19",
  TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
  OrgDb = "org.Hs.eg.db",
  outDir = tempdir(),
  cellPopLabel = "Clusters",
  cellPopulations = c("C2", "C5"),
  numCores = 1,
  verbose = FALSE
)

# ---- fixtures-sample-tile-matrix ----
sampleTileMatrices <- MOCHA::getSampleTileMatrix(
  tileResults,
  cellPopulations = c("C2", "C5"),
  reproducibilityThreshold = 0
)

# ---- fixtures-subset ----
c2Matrix <- MOCHA::subsetMOCHAObject(
  sampleTileMatrices,
  subsetBy = "celltype",
  groupList = "C2"
)

# ---- fixtures-annotate ----
annotated <- MOCHA::annotateTiles(c2Matrix)

# ---- fixtures-differentials ----
samples <- unique(SummarizedExperiment::colData(c2Matrix)$Sample)
diffs <- NULL
if (length(samples) >= 2L) {
  diffs <- suppressWarnings(
    MOCHA::getDifferentialAccessibleTiles(
      SampleTileObj = c2Matrix,
      cellPopulations = "C2",
      groupColumn = "Sample",
      foreground = samples[1],
      background = samples[2],
      numCores = 1,
      verbose = FALSE
    )
  )
}

# Synthetic multisample objects for downstream / modeling tutorials
tutorial_stm_multisample <- tutorial_synthetic_stm()
tutorial_dropout_stm <- tutorial_synthetic_dropout_stm()

# Alias for alt-TSS vignette variable names
differentials <- diffs
if (is.null(differentials) || length(differentials) == 0L) {
  tile_gr <- SummarizedExperiment::rowRanges(c2Matrix)[seq_len(min(20L, length(c2Matrix)))]
  differentials <- tile_gr
  differentials$FDR <- stats::runif(length(differentials), 0.01, 0.5)
  differentials$Log2FC_C <- stats::rnorm(length(differentials), 0, 1.5)
}

fixture_dir <- file.path(tempdir(), "mocha_tutorial_fixtures")
dir.create(fixture_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(tileResults, file.path(fixture_dir, "tileResults.rds"))
saveRDS(sampleTileMatrices, file.path(fixture_dir, "sampleTileMatrices.rds"))
saveRDS(c2Matrix, file.path(fixture_dir, "c2Matrix.rds"))
saveRDS(differentials, file.path(fixture_dir, "differentials.rds"))
saveRDS(tutorial_stm_multisample, file.path(fixture_dir, "tutorial_stm_multisample.rds"))
saveRDS(tutorial_dropout_stm, file.path(fixture_dir, "tutorial_dropout_stm.rds"))

invisible(list(
  tileResults = tileResults,
  sampleTileMatrices = sampleTileMatrices,
  c2Matrix = c2Matrix,
  annotated = annotated,
  diffs = diffs,
  differentials = differentials,
  tutorial_stm_multisample = tutorial_stm_multisample,
  tutorial_dropout_stm = tutorial_dropout_stm,
  fixture_dir = fixture_dir
))
