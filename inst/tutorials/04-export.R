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
  tileResults <- unpacked
}

# ---- export-coverage ----
# Heavy: writes bigWig files from fragments. Run with MOCHA_HEAVY_TESTS=true.
if (tolower(Sys.getenv("MOCHA_HEAVY_TESTS", "false")) %in% c("true", "1", "yes") &&
    requireNamespace("rtracklayer", quietly = TRUE)) {
  tutorial_try_run(
    suppressWarnings(
      exportCoverage(
        SampleTileObject = tileResults,
        dir = file.path(outDir, "sample_specific_coverage"),
        cellPopulations = c("C2", "C5"),
        sampleSpecific = TRUE,
        saveFile = TRUE,
        numCores = 1
      )
    ),
    label = "export-coverage"
  )
}

# ---- export-open-tiles ----
if (requireNamespace("rtracklayer", quietly = TRUE)) {
  tutorial_try_run(
    exportOpenTiles(
      SampleTileObject = sampleTileMatrices,
      cellPopulation = "C2",
      outDir = file.path(outDir, "tiles_samplespecific"),
      verbose = FALSE
    ),
    label = "export-open-tiles"
  )
}

# ---- export-diffs ----
if (requireNamespace("rtracklayer", quietly = TRUE) &&
    !is.null(diffs) && inherits(diffs, "GRanges") && length(diffs) > 0L) {
  tutorial_try_run(
    exportDifferentials(
      SampleTileObject = sampleTileMatrices,
      DifferentialsGRList = list(C2 = diffs),
      outDir = file.path(outDir, "tiles_differential"),
      verbose = FALSE
    ),
    label = "export-diffs"
  )
}

# ---- export-motifs ----
# Full CIS-BP motif export is slow; enable with MOCHA_HEAVY_TESTS=true.
if (tolower(Sys.getenv("MOCHA_HEAVY_TESTS", "false")) %in% c("true", "1", "yes") &&
    requireNamespace("chromVARmotifs", quietly = TRUE) &&
    requireNamespace("motifmatchr", quietly = TRUE) &&
    requireNamespace("rtracklayer", quietly = TRUE)) {
  tutorial_try_run({
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
  }, label = "export-motifs")
}

# ---- export-footprints ----
if (requireNamespace("rtracklayer", quietly = TRUE)) {
  tutorial_try_run(
    exportLocalFootprints(
      SampleTileObj = sampleTileMatrices,
      cellPopulation = "C2",
      outDir = file.path(outDir, "footprints"),
      numCores = 1
    ),
    label = "export-footprints"
  )
}

# ---- session-info ----
sessionInfo()
