test_that("plotConsensus returns reproducibility data.frames", {
  dfs <- MOCHA::plotConsensus(
    MOCHA:::testTileResultsMultisample,
    cellPopulations = c("C2", "C3"),
    returnDFs = TRUE,
    numCores = 1
  )

  expect_type(dfs, "list")
  expect_true(all(c("C2", "C3") %in% names(dfs)))
  for (nm in names(dfs)) {
    expect_true(all(c("Reproducibility", "PeakNumber") %in% colnames(dfs[[nm]])))
    expect_gt(nrow(dfs[[nm]]), 0L)
  }
})

test_that("plotConsensus errors for unknown cell populations", {
  expect_error(
    MOCHA::plotConsensus(
      MOCHA:::testTileResultsMultisample,
      cellPopulations = "NOT_A_POP",
      returnDFs = TRUE
    ),
    "must present in tileResults"
  )
})

test_that("plotConsensus returns ggplot objects when requested", {
  skip_if_not_installed("ggplot2")
  plots <- MOCHA::plotConsensus(
    MOCHA:::testTileResultsMultisample,
    cellPopulations = "C2",
    returnPlotList = TRUE,
    numCores = 1
  )
  expect_type(plots, "list")
  expect_s3_class(plots[[1]], "ggplot")
})
