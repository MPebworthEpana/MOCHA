test_that("updateDirectoryPath updates Directory metadata on MOCHA objects", {
  capture.output(
    stm <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResults,
      cellPopulations = "C2",
      threshold = 0,
      numCores = 1
    ),
    type = "message"
  )

  new_dir <- normalizePath(tempdir(), winslash = "/")
  updated <- MOCHA::updateDirectoryPath(stm, directoryPath = new_dir)
  expect_equal(S4Vectors::metadata(updated)$Directory, new_dir)
})

test_that("updateDirectoryPath warns on non-MOCHA objects", {
  expect_warning(
    out <- MOCHA::updateDirectoryPath(data.frame(x = 1), directoryPath = tempdir()),
    "not a MOCHA"
  )
  expect_null(out)
})
