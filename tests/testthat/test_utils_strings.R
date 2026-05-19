test_that("StringsToGRanges and GRangesToString round-trip", {
  regions <- c("chr1:1000-1499", "chr2:2000-2499")
  gr <- MOCHA::StringsToGRanges(regions)
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 2L)
  back <- MOCHA::GRangesToString(gr)
  expect_equal(back, regions)
})

test_that("StringsToGRanges errors on invalid strings", {
  expect_error(MOCHA::StringsToGRanges("not-a-region"))
})

test_that("isMOCHAObject identifies MOCHA object types", {
  expect_true(MOCHA::isMOCHAObject(MOCHA:::testTileResults))
  expect_equal(
    MOCHA::isMOCHAObject(MOCHA:::testTileResults, returnType = TRUE),
    "OpenTiles"
  )

  capture.output(
    stm <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResults,
      cellPopulations = c("C2", "C5"),
      threshold = 0,
      numCores = 1
    ),
    type = "message"
  )
  expect_true(MOCHA::isMOCHAObject(stm))
  expect_equal(MOCHA::isMOCHAObject(stm, returnType = TRUE), "SampleTileMatrix")
  expect_false(MOCHA::isMOCHAObject(data.frame(x = 1)))
})
