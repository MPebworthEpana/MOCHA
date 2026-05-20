test_that("getCellPopMatrix returns sample by tile matrix for a cell population", {
  capture.output(
    stm <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResults,
      cellPopulations = "C2",
      reproducibilityThreshold = 0,
      numCores = 1
    ),
    type = "message"
  )

  mat <- MOCHA::getCellPopMatrix(stm, "C2", dropSamples = TRUE, NAtoZero = TRUE)
  expect_true(is.matrix(mat))
  expect_gt(nrow(mat), 0L)
  expect_gt(ncol(mat), 0L)
})

test_that("getCellPopMatrix errors when cell population is missing", {
  capture.output(
    stm <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResults,
      cellPopulations = "C2",
      reproducibilityThreshold = 0,
      numCores = 1
    ),
    type = "message"
  )

  expect_error(
    MOCHA::getCellPopMatrix(stm, "NOT_FOUND"),
    "Cell population not found"
  )
})
