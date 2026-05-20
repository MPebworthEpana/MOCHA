test_that("filterCoAccessibleLinks filters by correlation and adds coordinates", {
  capture.output(
    SampleTileMatrix <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResults,
      cellPopulations = "C2",
      reproducibilityThreshold = 0
    )
  )
  links <- MOCHA::getCoAccessibleLinks(
    SampleTileMatrix,
    "C2",
    MOCHA::StringsToGRanges("chr1:101873000-101873499"),
    verbose = FALSE
  )
  expect_gt(nrow(links), 0L)
  skip_if(
    !any(abs(links$Correlation) > 0, na.rm = TRUE),
    "Fixture produced no non-zero correlations for filterCoAccessibleLinks"
  )

  filtered <- MOCHA::filterCoAccessibleLinks(links, threshold = 0)
  expect_true(all(abs(filtered$Correlation) > 0))
  expect_true(all(c("chr", "start", "end") %in% colnames(filtered)))
  expect_equal(nrow(filtered), nrow(links[abs(links$Correlation) > 0, , drop = FALSE]))
})

test_that("filterCoAccessibleLinks errors when no values pass threshold", {
  links <- data.frame(
    Tile1 = "chr1:1000-1499",
    Tile2 = "chr1:2000-2499",
    Correlation = c(0.1, -0.1),
    stringsAsFactors = FALSE
  )
  expect_error(
    MOCHA::filterCoAccessibleLinks(links, threshold = 0.5),
    "no values above the threshold"
  )
})
