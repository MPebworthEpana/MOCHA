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
if (!exists("sampleTileMatrices")) {
  fixture_dir <- file.path(tempdir(), "mocha_tutorial_fixtures")
  if (!file.exists(file.path(fixture_dir, "sampleTileMatrices.rds"))) {
    stop("Source 00-fixtures.R before 05-advanced-modeling.R")
  }
  sampleTileMatrices <- readRDS(file.path(fixture_dir, "sampleTileMatrices.rds"))
  tutorial_stm_multisample <- readRDS(file.path(fixture_dir, "tutorial_stm_multisample.rds"))
}
modelFormula <- "exp ~ GroupA + PassQC + (1|Sample)"

# ---- pilot-lmem ----
if (requireNamespace("lmerTest", quietly = TRUE)) {
  tutorial_try_run({
    pilot <- pilotLMEM(
      ExperimentObj = tutorial_stm_multisample,
      modelFormula = modelFormula,
      assayName = "C2",
      pilotIndices = 1:5,
      verbose = FALSE
    )
  }, "pilotLMEM")
}

# ---- run-lmem ----
if (requireNamespace("lmerTest", quietly = TRUE)) {
  tutorial_try_run({
    modelList <- runLMEM(
      ExperimentObj = tutorial_stm_multisample,
      modelFormula = modelFormula,
      assayName = "C2",
      initialSampling = 5,
      numCores = 1,
      verbose = FALSE
    )
    coefs <- getModelValues(modelList, value = "Estimate")
    head(coefs)
  }, "runLMEM")
}

# ---- pilot-ziglmm ----
if (requireNamespace("glmmTMB", quietly = TRUE)) {
  tutorial_try_run({
    zigPilot <- pilotZIGLMM(
      TSAM_Object = tutorial_stm_multisample,
      cellPopulation = "C2",
      continuousFormula = exp ~ GroupA,
      ziformula = ~ GroupA,
      verbose = FALSE
    )
  }, "pilotZIGLMM")
}

# ---- run-ziglmm ----
if (requireNamespace("glmmTMB", quietly = TRUE)) {
  tutorial_try_run({
    zigModels <- runZIGLMM(
      TSAM_Object = tutorial_stm_multisample,
      cellPopulation = "C2",
      continuousFormula = exp ~ GroupA,
      ziformula = ~ GroupA,
      numCores = 1
    )
    variances <- varZIGLMM(zigModels)
  }, "runZIGLMM")
}

# ---- train-peak-model ----
# Slow; run with MOCHA_HEAVY_TESTS=true for a full training pass.
if (tolower(Sys.getenv("MOCHA_HEAVY_TESTS", "false")) %in% c("true", "1", "yes")) {
  tutorial_try_run({
    frags <- exampleFragments[[1]]
    peaks <- GenomicRanges::reduce(frags)
    peaks <- peaks[IRanges::width(peaks) >= 200L]
    customPeaks <- trainPeakModel(
      ATACFragments = frags,
      cellColData = exampleCellColData,
      blackList = exampleBlackList,
      groundTruthPeaks = peaks,
      tileSize = 250L,
      cellSubsetSizes = c(50L, 100L),
      replicatesFn = function(n) 2L,
      threshMethod = if (requireNamespace("cutpointr", quietly = TRUE)) "youden" else "f1",
      numCores = 1L,
      seed = 42L,
      verbose = FALSE
    )
  }, "trainPeakModel")
}

# ---- s4-classes ----
isMOCHAObject(sampleTileMatrices)

# ---- session-info ----
sessionInfo()
