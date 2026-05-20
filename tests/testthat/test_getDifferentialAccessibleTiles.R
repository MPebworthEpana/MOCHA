test_that("getDifferentialAccessibleTiles validates inputs", {
  capture.output(
    SampleTileMatrix <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResultsMultisample,
      cellPopulations = c("C2", "C3"),
      reproducibilityThreshold = 0,
      numCores = 1
    ),
    type = "message"
  )

  expect_error(
    MOCHA::getDifferentialAccessibleTiles(
      SampleTileMatrix,
      cellPopulation = "NOT_A_POP",
      groupColumn = "Sample",
      foreground = "scATAC_BMMC_R1",
      background = "scATAC_CD34_BMMC_R1"
    ),
    "cellPopulation was not found"
  )

  expect_error(
    MOCHA::getDifferentialAccessibleTiles(
      SampleTileMatrix,
      cellPopulation = "C2",
      groupColumn = "NOT_A_COLUMN",
      foreground = "a",
      background = "b"
    ),
    "not found in the provided SampleTileObj"
  )

  expect_error(
    MOCHA::getDifferentialAccessibleTiles(
      SampleTileMatrix,
      cellPopulation = "C2",
      groupColumn = "Sample",
      foreground = "NOT_A_SAMPLE",
      background = "scATAC_CD34_BMMC_R1"
    ),
    "foreground value is not present"
  )

  expect_error(
    MOCHA::getDifferentialAccessibleTiles(
      SampleTileMatrix,
      cellPopulation = "C2",
      groupColumn = "Sample",
      foreground = "scATAC_BMMC_R1",
      background = "NOT_A_SAMPLE"
    ),
    "background value is not present"
  )

  expect_error(
    MOCHA::getDifferentialAccessibleTiles(
      SampleTileMatrix,
      cellPopulation = "C2",
      groupColumn = "Sample",
      foreground = "scATAC_BMMC_R1",
      background = "scATAC_CD34_BMMC_R1",
      qValueMethod = "invalid"
    ),
    "qValueMethod must either be set"
  )
})

test_that("getDifferentialAccessibleTiles returns NULL when fewer than three samples per group", {
  capture.output(
    SampleTileMatrix <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResultsMultisample,
      cellPopulations = "C3",
      reproducibilityThreshold = 0,
      numCores = 1
    ),
    type = "message"
  )

  expect_message(
    out <- MOCHA::getDifferentialAccessibleTiles(
      SampleTileMatrix,
      cellPopulation = "C3",
      groupColumn = "Sample",
      foreground = "scATAC_BMMC_R1",
      background = "scATAC_CD34_BMMC_R1",
      signalThreshold = 5,
      outputGRanges = FALSE,
      numCores = 1,
      verbose = FALSE
    ),
    "Less than three samples available per group"
  )
  expect_true(is.null(out) || (is.data.frame(out) && nrow(out) == 0))
})

test_that("getDifferentialAccessibleTiles runs with fixed signalThreshold on synthetic data", {
  SampleTileMatrix <- make_synthetic_sample_tile_matrix(n_tiles = 40L, n_per_group = 4L)

  suppressWarnings(
    differentials <- MOCHA::getDifferentialAccessibleTiles(
      SampleTileMatrix,
      cellPopulation = "C2",
      groupColumn = "GroupA",
      foreground = "A",
      background = "B",
      signalThreshold = 8,
      qValueMethod = "standard",
      outputGRanges = FALSE,
      numCores = 1,
      verbose = FALSE
    )
  )

  expect_true(is.data.frame(differentials))
  expect_true(all(c("Tile", "P_value", "FDR", "CellPopulation") %in% colnames(differentials)))
  expect_gt(nrow(differentials), 0L)
})

test_that("positional arguments after background match legacy order", {
  SampleTileMatrix <- make_synthetic_sample_tile_matrix(n_tiles = 40L, n_per_group = 4L)

  named <- suppressWarnings(
    MOCHA::getDifferentialAccessibleTiles(
      SampleTileMatrix,
      cellPopulation = "C2",
      groupColumn = "GroupA",
      foreground = "A",
      background = "B",
      minZeroDiff = 0.5,
      qValueMethod = "standard",
      signalThreshold = 8,
      qValueThreshold = 0.2,
      outputGRanges = FALSE,
      numCores = 1,
      verbose = FALSE
    )
  )

  positional <- suppressWarnings(
    MOCHA::getDifferentialAccessibleTiles(
      SampleTileMatrix,
      "C2",
      "GroupA",
      "A",
      "B",
      0.5,
      "standard",
      8,
      0.2,
      FALSE,
      1,
      FALSE
    )
  )

  expect_equal(positional$Tile, named$Tile)
  expect_equal(positional$P_value, named$P_value)
})

test_that("deprecated techThreshold and fdrToDisplay aliases work", {
  SampleTileMatrix <- make_synthetic_sample_tile_matrix(n_tiles = 40L, n_per_group = 4L)

  expect_warning(
    out_bio <- MOCHA::getDifferentialAccessibleTiles(
      SampleTileMatrix,
      cellPopulation = "C2",
      groupColumn = "GroupA",
      foreground = "A",
      background = "B",
      signalThreshold = 8,
      qValueMethod = "standard",
      techThreshold = 0.8,
      outputGRanges = FALSE,
      numCores = 1,
      verbose = FALSE
    ),
    "techThreshold"
  )
  expect_true(is.data.frame(out_bio))

  expect_warning(
    out_fdr <- MOCHA::getDifferentialAccessibleTiles(
      SampleTileMatrix,
      cellPopulation = "C2",
      groupColumn = "GroupA",
      foreground = "A",
      background = "B",
      signalThreshold = 8,
      qValueMethod = "standard",
      fdrToDisplay = 0.2,
      outputGRanges = FALSE,
      numCores = 1,
      verbose = FALSE
    ),
    "fdrToDisplay"
  )
  expect_true(is.data.frame(out_fdr))
})

test_that("dropoutAdjustment auto-runs assessDropout when assay missing", {
  stm <- make_synthetic_sample_tile_matrix(n_tiles = 40L, n_per_group = 4L)
  cd <- SummarizedExperiment::colData(stm)
  cd$FragNumber <- 5000L
  cd$CellCounts <- 100L
  SummarizedExperiment::colData(stm) <- cd
  expect_false("DropoutProb_C2" %in% SummarizedExperiment::assayNames(stm))

  expect_message(
    suppressWarnings(
      MOCHA::getDifferentialAccessibleTiles(
        stm,
        cellPopulation = "C2",
        groupColumn = "GroupA",
        foreground = "A",
        background = "B",
        signalThreshold = 8,
        qValueMethod = "standard",
        dropoutAdjustment = "biological_only",
        bioThreshold = 0.8,
        outputGRanges = FALSE,
        numCores = 1,
        verbose = FALSE
      )
    ),
    "assessDropout|dropout probability"
  )
})

test_that("dropoutAdjustment biological_only integrates with assessDropout", {
  stm <- make_synthetic_sample_tile_matrix(n_tiles = 40L, n_per_group = 4L)
  cd <- SummarizedExperiment::colData(stm)
  cd$FragNumber <- 5000L
  cd$CellCounts <- 100L
  SummarizedExperiment::colData(stm) <- cd
  stm <- MOCHA::assessDropout(stm)

  out <- suppressWarnings(
    MOCHA::getDifferentialAccessibleTiles(
      stm,
      cellPopulation = "C2",
      groupColumn = "GroupA",
      foreground = "A",
      background = "B",
      signalThreshold = 8,
      qValueMethod = "standard",
      dropoutAdjustment = "biological_only",
      bioThreshold = 0.8,
      outputGRanges = FALSE,
      numCores = 1,
      verbose = FALSE
    )
  )

  expect_true(is.data.frame(out))
  expect_gt(nrow(out), 0L)
})

test_that("dropoutAdjustment fails clearly when model cannot be fitted", {
  stm <- make_synthetic_dropout_tsam(n_bio_closed = 1L, n_tech_dropout = 1L, n_per_group = 2L)

  expect_error(
    MOCHA::getDifferentialAccessibleTiles(
      stm,
      cellPopulation = "C2",
      groupColumn = "GroupA",
      foreground = "A",
      background = "B",
      signalThreshold = 8,
      qValueMethod = "standard",
      dropoutAdjustment = "biological_only",
      outputGRanges = FALSE,
      numCores = 1,
      verbose = FALSE
    ),
    "DropoutProb|Dropout-adjusted|Insufficient data"
  )
})

test_that("getDifferentialAccessibleTiles can return GRanges", {
  SampleTileMatrix <- make_synthetic_sample_tile_matrix(n_tiles = 30L, n_per_group = 4L)

  suppressWarnings(
    gr <- MOCHA::getDifferentialAccessibleTiles(
      SampleTileMatrix,
      cellPopulation = "C2",
      groupColumn = "GroupA",
      foreground = "A",
      background = "B",
      signalThreshold = 8,
      qValueMethod = "standard",
      outputGRanges = TRUE,
      numCores = 1,
      verbose = FALSE
    )
  )

  expect_s4_class(gr, "GRanges")
})
