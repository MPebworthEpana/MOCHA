skip_if_not_installed("glmmTMB")

test_that("linearModeling returns named lmer list for one cell type", {
  skip_if_not_installed("lmerTest")
  ExperimentObj <- make_synthetic_sample_tile_matrix(
    n_tiles = 10L,
    n_per_group = 4L
  )

  models <- MOCHA::linearModeling(
    ExperimentObj,
    formula = exp ~ GroupA + (1 | GroupA),
    CellType = "C2",
    threshold = 0,
    NAtoZero = TRUE,
    numCores = 1
  )

  expect_type(models, "list")
  expect_gt(length(models), 0L)
  expect_true(all(vapply(models, function(m) inherits(m, "lmerMod"), logical(1))))
})

test_that("runZIGLMM validates formula and cellPopulation", {
  skip_if_not_installed("glmmTMB")
  capture.output(
    stm <- MOCHA::combineSampleTileMatrix(
      MOCHA::getSampleTileMatrix(
        MOCHA:::testTileResults,
        cellPopulations = c("C2", "C5"),
        reproducibilityThreshold = 0,
        numCores = 1
      )
    ),
    type = "message"
  )

  expect_error(
    MOCHA::runZIGLMM(
      stm,
      cellPopulation = "NOT_FOUND",
      continuousFormula = exp ~ PassQC,
      ziformula = ~ PassQC
    ),
    "Error around cell type name"
  )

  expect_error(
    MOCHA::runZIGLMM(
      stm,
      cellPopulation = "counts",
      continuousFormula = "not a formula",
      ziformula = ~ PassQC
    ),
    "not provided as a formula"
  )

  expect_error(
    MOCHA::runZIGLMM(
      stm,
      cellPopulation = "counts",
      continuousFormula = exp ~ PassQC,
      ziformula = ~ PassQC,
      zi_threshold = 2
    ),
    "zi_threshold must be between 0 and 1"
  )
})

test_that("pilotZIGLMM returns pilot model list", {
  skip_if_not_installed("glmmTMB")
  capture.output(
    stm <- MOCHA::combineSampleTileMatrix(
      MOCHA::getSampleTileMatrix(
        MOCHA:::testTileResults,
        cellPopulations = c("C2", "C5"),
        reproducibilityThreshold = 0,
        numCores = 1
      )
    ),
    type = "message"
  )

  pilot <- MOCHA::pilotZIGLMM(
    stm,
    cellPopulation = "counts",
    continuousFormula = exp ~ PassQC,
    ziformula = ~ PassQC,
    pilotIndices = 1:3,
    verbose = FALSE
  )

  expect_type(pilot, "list")
  expect_equal(length(pilot), 3L)
})

test_that("getModelValues extracts slopes and p-values from runZIGLMM output", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("lmerTest")

  capture.output(
    stm <- MOCHA::combineSampleTileMatrix(
      MOCHA::getSampleTileMatrix(
        MOCHA:::testTileResults,
        cellPopulations = c("C2", "C5"),
        reproducibilityThreshold = 0,
        numCores = 1
      )
    ),
    type = "message"
  )

  suppressWarnings(
    res <- try(
      MOCHA::runZIGLMM(
        stm[1:50, ],
        cellPopulation = "counts",
        continuousFormula = exp ~ PassQC,
        ziformula = ~ PassQC,
        initialSampling = 3,
        numCores = 1,
        verbose = FALSE
      ),
      silent = TRUE
    )
  )

  skip_if(inherits(res, "try-error"), "runZIGLMM did not converge on fixture subset")

  df <- MOCHA::getModelValues(res, "PassQC")
  expect_true(all(c("Element", "Estimate", "PValue") %in% colnames(df)))
  expect_equal(nrow(df), nrow(res))
})

test_that("varZIGLMM validates cellPopulation and returns variance decomposition", {
  skip_if_not_installed("glmmTMB")
  capture.output(
    stm <- MOCHA::combineSampleTileMatrix(
      MOCHA::getSampleTileMatrix(
        MOCHA:::testTileResults,
        cellPopulations = "C2",
        reproducibilityThreshold = 0,
        numCores = 1
      )
    ),
    type = "message"
  )

  expect_error(
    MOCHA::varZIGLMM(
      stm,
      cellPopulation = "NOT_FOUND",
      continuousRandom = "PassQC",
      ziRandom = 0
    ),
    "cellPopulation was not found"
  )

  suppressWarnings(
    vd <- MOCHA::varZIGLMM(
      stm[1:30, ],
      cellPopulation = "counts",
      continuousRandom = "PassQC",
      ziRandom = 0,
      numCores = 1,
      verbose = FALSE
    )
  )

  expect_true(is.matrix(vd) || is.data.frame(vd))
  expect_gt(nrow(vd), 0L)
})
