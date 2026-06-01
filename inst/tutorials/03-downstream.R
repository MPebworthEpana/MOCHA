# ---- libraries ----
library(MOCHA)

# ---- script-init ----
if (!exists("tutorial_try_run", mode = "function")) {
  shared_path <- system.file("tutorials", "_shared.R", package = "MOCHA")
  if (!nzchar(shared_path)) {
    shared_path <- file.path(getwd(), "inst", "tutorials", "_shared.R")
  }
  if (file.exists(shared_path)) {
    source(shared_path, local = FALSE)
  }
}
if (!exists("tutorial_stm_multisample")) {
  fixture_dir <- file.path(tempdir(), "mocha_tutorial_fixtures")
  if (!file.exists(file.path(fixture_dir, "sampleTileMatrices.rds"))) {
    stop("Source 00-fixtures.R before 03-downstream.R")
  }
  tileResults <- readRDS(file.path(fixture_dir, "tileResults.rds"))
  sampleTileMatrices <- readRDS(file.path(fixture_dir, "sampleTileMatrices.rds"))
  c2Matrix <- readRDS(file.path(fixture_dir, "c2Matrix.rds"))
  tutorial_stm_multisample <- readRDS(file.path(fixture_dir, "tutorial_stm_multisample.rds"))
  tutorial_dropout_stm <- readRDS(file.path(fixture_dir, "tutorial_dropout_stm.rds"))
}
stm_multi <- tutorial_stm_multisample
region_name <- rownames(stm_multi)[1L]

# ---- co-access-links ----
regions <- StringsToGRanges(region_name)

links <- getCoAccessibleLinks(
  SampleTileObj = stm_multi,
  cellPopulation = "C2",
  regions = regions,
  verbose = FALSE
)
head(links)

# ---- filter-co-access ----
filtered <- filterCoAccessibleLinks(
  links,
  threshold = 0.5
)

# ---- test-co-access ----
if (nrow(filtered) > 0L) {
  coAccessTests <- testCoAccessibility(
    SampleTileObj = stm_multi,
    tile1 = filtered$Tile1,
    tile2 = filtered$Tile2,
    numCores = 1,
    verbose = FALSE
  )
}

# ---- co-access-not-cran ----
# Optional local run with objects from the workflow tutorial:
# links <- getCoAccessibleLinks(sampleTileMatrices, "C2", regions, numCores = 1)

# ---- assess-dropout ----
stm_with_dropout <- assessDropout(
  TSAM_Object = tutorial_dropout_stm,
  cellPopulation = "C2",
  verbose = FALSE
)

# ---- dropout-model ----
model <- estimateDropoutModel(stm_with_dropout, cellPopulation = "C2")
probs <- getDropoutProb(stm_with_dropout, "C2")
cls <- classifyZeros(stm_with_dropout, "C2")

# ---- dropout-plots ----
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plotDropoutDiagnostics(stm_with_dropout, "C2", type = "probHist")
  plotDropoutDiagnostics(stm_with_dropout, "C2", type = "calibration")
}

# ---- bulk-dr ----
if (requireNamespace("uwot", quietly = TRUE) &&
    requireNamespace("irlba", quietly = TRUE)) {
  lse <- bulkDimReduction(
    SampleTileObj = stm_multi,
    cellType = "All",
    componentNumber = 2
  )
  umap_coords <- bulkUMAP(lse, components = 1:2, nNeighbors = 4, verbose = FALSE)
}

# ---- merge-tiles ----
# tileResults2 from a second cohort (same genome/annotation)
# merged <- mergeTileResults(
#   tileResultsList = list(tileResults, tileResults2),
#   cellPopulations = c("C2", "C5")
# )

# ---- combine-stm ----
combined <- tutorial_try_run(
  combineSampleTileMatrix(sampleTileMatrices),
  label = "combine-stm"
)
if (!is.null(combined)) {
  dim(combined)
}

# ---- getters ----
getCellTypes(sampleTileMatrices)
getCellTypeTiles(sampleTileMatrices, cellType = "C2")
getSampleCellTypeMetadata(sampleTileMatrices)

# ---- rename-add-subset ----
# renamed <- renameCellTypes(sampleTileMatrices, oldName = "C2", newName = "Mono")
# withCol <- addCellColData(tileResults, newColData = data.frame(...))
c2_only <- subsetMOCHAObject(
  sampleTileMatrices,
  subsetBy = "celltype",
  groupList = "C2"
)

# ---- session-info ----
sessionInfo()
