# Export and sharing tutorial — executable chunks (temp directories only).

# ---- libraries ----
library(MOCHA)

if (!exists("tileResults")) {
  fixture_dir <- file.path(tempdir(), "mocha_tutorial_fixtures")
  if (!file.exists(file.path(fixture_dir, "tileResults.rds"))) {
    stop("Source 00-fixtures.R before 04-export.R")
  }
  tileResults <- readRDS(file.path(fixture_dir, "tileResults.rds"))
  sampleTileMatrices <- readRDS(file.path(fixture_dir, "sampleTileMatrices.rds"))
  diffs <- readRDS(file.path(fixture_dir, "differentials.rds"))
}

outDir <- file.path(tempdir(), "mocha_export_tutorial")
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

# ---- pack-unpack ----
if (requireNamespace("zip", quietly = TRUE)) {
  zipPath <- packMOCHA(
    MOCHAObj = tileResults,
    zipfile = file.path(outDir, "mocha_results.zip")
  )
  unpacked <- unpackMOCHA(
    zipfile = zipPath,
    exdir = file.path(outDir, "unpacked_mocha")
  )
}

# ---- update-path ----
if (exists("unpacked")) {
  tileResults <- updateDirectoryPath(
    tileResults,
    directoryPath = file.path(outDir, "unpacked_mocha")
  )
}

# ---- export-coverage ----
if (requireNamespace("rtracklayer", quietly = TRUE)) {
  suppressWarnings(
    exportCoverage(
      SampleTileObject = tileResults,
      dir = file.path(outDir, "sample_specific_coverage"),
      cellPopulations = c("C2", "C5"),
      sampleSpecific = TRUE,
      saveFile = TRUE,
      numCores = 1
    )
  )
}

# ---- export-open-tiles ----
if (requireNamespace("rtracklayer", quietly = TRUE)) {
  exportOpenTiles(
    SampleTileObject = sampleTileMatrices,
    cellPopulation = "C2",
    outDir = file.path(outDir, "tiles_samplespecific"),
    verbose = FALSE
  )
}

# ---- export-diffs ----
if (requireNamespace("rtracklayer", quietly = TRUE) &&
    !is.null(diffs) && inherits(diffs, "GRanges") && length(diffs) > 0L) {
  exportDifferentials(
    SampleTileObject = sampleTileMatrices,
    DifferentialsGRList = list(C2 = diffs),
    outDir = file.path(outDir, "tiles_differential"),
    verbose = FALSE
  )
}

# ---- export-motifs ----
if (requireNamespace("chromVARmotifs", quietly = TRUE) &&
    requireNamespace("motifmatchr", quietly = TRUE) &&
    requireNamespace("rtracklayer", quietly = TRUE)) {
  motifsGRanges <- addMotifSet(
    SampleTileObj = sampleTileMatrices,
    motifPWMs = chromVARmotifs::human_pwms_v2,
    returnSTM = FALSE,
    motifSetName = "CISBP"
  )
  exportMotifs(
    SampleTileObject = tileResults,
    motifsGRanges = unlist(motifsGRanges),
    motifSetName = "CISBP",
    outDir = file.path(outDir, "motifs"),
    filterCellTypePeaks = TRUE,
    verbose = FALSE
  )
}

# ---- export-footprints ----
if (requireNamespace("rtracklayer", quietly = TRUE)) {
  exportLocalFootprints(
    SampleTileObj = sampleTileMatrices,
    cellPopulations = "C2",
    outDir = file.path(outDir, "footprints"),
    numCores = 1
  )
}

# ---- session-info ----
sessionInfo()
