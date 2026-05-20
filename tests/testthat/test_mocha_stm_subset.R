test_that("subset() on MochaSampleTileMatrix matches subsetMOCHAObject", {
  capture.output(
    stm <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResults,
      cellPopulations = c("C2", "C5"),
      reproducibilityThreshold = 0,
      numCores = 1,
      returnClass = "mocha"
    ),
    type = "message"
  )

  legacy <- MOCHA::getSampleTileMatrix(
    MOCHA:::testTileResults,
    cellPopulations = c("C2", "C5"),
    reproducibilityThreshold = 0,
    numCores = 1
  )

  by_fn <- MOCHA::subsetMOCHAObject(
    legacy,
    subsetBy = "celltype",
    groupList = "C2",
    subsetPeaks = TRUE
  )
  by_method <- subset(stm, cells = "C2", subsetPeaks = TRUE)

  expect_equal(SummarizedExperiment::assayNames(by_method), SummarizedExperiment::assayNames(by_fn))
  expect_equal(dim(by_method), dim(by_fn))
  expect_equal(
    rownames(SummarizedExperiment::assays(by_method)[[1]]),
    rownames(SummarizedExperiment::assays(by_fn)[[1]])
  )
  expect_s4_class(by_method, "MochaSampleTileMatrix")
})

test_that("stm[, j] preserves MochaSampleTileMatrix class", {
  capture.output(
    stm <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResults,
      cellPopulations = "C2",
      reproducibilityThreshold = 0,
      numCores = 1,
      returnClass = "mocha"
    ),
    type = "message"
  )
  j <- colnames(stm)[1]
  sub <- stm[, j]
  expect_s4_class(sub, "MochaSampleTileMatrix")
  expect_equal(ncol(sub), 1L)
})
