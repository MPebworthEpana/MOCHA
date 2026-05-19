test_that("Mocha coercion and validity work on bundled tileResults", {
  tr <- MOCHA:::testTileResults
  mtr <- MOCHA::asMochaTileResults(tr)
  expect_s4_class(mtr, "MochaTileResults")
  expect_s4_class(mtr, "MultiAssayExperiment")
  expect_true(MOCHA::isMOCHAObject(mtr))
  expect_equal(MOCHA::cellTypes(mtr), MOCHA::getCellTypes(mtr))
  expect_s4_class(MOCHA::openTiles(mtr, returnType = "GRangesList"), "GRangesList")

  capture.output(
    stm <- MOCHA::getSampleTileMatrix(
      tr,
      cellPopulations = c("C2", "C5"),
      threshold = 0,
      numCores = 1,
      returnClass = "mocha"
    ),
    type = "message"
  )
  expect_s4_class(stm, "MochaSampleTileMatrix")
  expect_s4_class(stm, "RangedSummarizedExperiment")
  expect_equal(MOCHA::cellTypes(stm), c("C2", "C5"))
  mstm <- MOCHA::asMochaSTM(stm)
  expect_s4_class(mstm, "MochaSampleTileMatrix")

  tmp <- tempfile(fileext = ".rds")
  saveRDS(mtr, tmp)
  loaded <- readRDS(tmp)
  expect_s4_class(loaded, "MochaTileResults")
  unlink(tmp)
})

test_that("returnClass legacy default is unchanged", {
  capture.output(
    stm <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResults,
      cellPopulations = c("C2", "C5"),
      threshold = 0,
      numCores = 1
    ),
    type = "message"
  )
  expect_s4_class(stm, "RangedSummarizedExperiment")
  expect_false(methods::is(stm, "MochaSampleTileMatrix"))
})
