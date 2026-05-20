test_that("getSampleTileMatrices works on a 1 sample test dataset", {
  cellPopulations <- c("C2", "C5")

  capture.output(
    tilemat <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResults,
      cellPopulations = cellPopulations,
      reproducibilityThreshold = 0,
      numCores = 1
    ),
    type = "message"
  )

  expect_s4_class(
    tilemat,
    "RangedSummarizedExperiment"
  )

  expect_snapshot_value(
    dim(tilemat),
    style = "json2"
  )

  expect_equal(
    names(SummarizedExperiment::assays(tilemat)),
    c("C2", "C5")
  )

  expect_true("summarizedData" %in% names(S4Vectors::metadata(tilemat)))
  expect_true(all(c("Sample", "PassQC") %in% colnames(SummarizedExperiment::colData(tilemat))))

  expect_snapshot_output(
    SummarizedExperiment::colData(tilemat)
  )
})
