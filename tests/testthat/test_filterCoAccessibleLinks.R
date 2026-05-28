test_that("filterCoAccessibleLinks filters by correlation and adds coordinates", {
  SampleTileMatrix <- make_synthetic_sample_tile_matrix(
    n_tiles = 40L,
    n_per_group = 4L
  )
  seed_tile <- SummarizedExperiment::rownames(SampleTileMatrix)[1L]
  regions <- MOCHA::StringsToGRanges(seed_tile)

  links <- MOCHA::getCoAccessibleLinks(
    SampleTileMatrix,
    "C2",
    regions,
    verbose = FALSE,
    ZI = FALSE
  )
  expect_gt(nrow(links), 0L)
  expect_true(any(abs(links$Correlation) > 0, na.rm = TRUE))

  links_zi <- MOCHA::getCoAccessibleLinks(
    SampleTileMatrix,
    "C2",
    regions,
    verbose = FALSE,
    ZI = TRUE
  )
  if (nrow(links_zi) > 0L) {
    expect_true(any(!is.na(links_zi$Correlation)))
  }

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
