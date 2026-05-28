# Advanced modeling tutorial — executable chunks.

# ---- libraries ----
library(MOCHA)

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

# ---- pilot-lmem ----
modelFormula <- "exp ~ GroupA + PassQC + (1|Sample)"
if (requireNamespace("lmerTest", quietly = TRUE) && exists("tutorial_stm_multisample")) {
  tutorial_try_run({
    pilot <- pilotLMEM(
      SampleTileObj = tutorial_stm_multisample,
      modelFormula = modelFormula,
      assayName = "C2",
      initialSampling = 5,
      numCores = 1,
      verbose = FALSE
    )
  }, "pilotLMEM")
}

# ---- run-lmem ----
if (requireNamespace("lmerTest", quietly = TRUE) && exists("tutorial_stm_multisample")) {
  tutorial_try_run({
    modelList <- runLMEM(
      SampleTileObj = tutorial_stm_multisample,
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
if (requireNamespace("glmmTMB", quietly = TRUE) && exists("tutorial_stm_multisample")) {
  tutorial_try_run({
    zigPilot <- pilotZIGLMM(
      SampleTileObj = tutorial_stm_multisample,
      modelFormula = "exp ~ GroupA + PassQC + (1|Sample)",
      assayName = "C2",
      initialSampling = 5,
      verbose = FALSE
    )
  }, "pilotZIGLMM")
}

# ---- run-ziglmm ----
if (requireNamespace("glmmTMB", quietly = TRUE) && exists("tutorial_stm_multisample")) {
  tutorial_try_run({
    zigModels <- runZIGLMM(
      SampleTileObj = tutorial_stm_multisample,
      modelFormula = "exp ~ GroupA + PassQC + (1|Sample)",
      assayName = "C2",
      numCores = 1
    )
    variances <- varZIGLMM(zigModels)
  }, "runZIGLMM")
}

# ---- train-peak-model ----
customPeaks <- trainPeakModel(
  ATACFragments = exampleFragments,
  cellColData = exampleCellColData,
  cellPopLabel = "Clusters",
  cellPopulations = "C2",
  genome = "BSgenome.Hsapiens.UCSC.hg19",
  outDir = tempdir(),
  numCores = 1,
  verbose = FALSE
)

# ---- s4-classes ----
isMOCHAObject(sampleTileMatrices)

# ---- session-info ----
sessionInfo()
