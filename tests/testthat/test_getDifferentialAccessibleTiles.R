test_that("getDifferentialAccessibleTiles validates inputs", {
  capture.output(
    SampleTileMatrix <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResultsMultisample,
      cellPopulations = c("C2", "C3"),
      threshold = 0,
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
      threshold = 0,
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
