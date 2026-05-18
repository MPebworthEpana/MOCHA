test_that(".resolveSubSamples returns all samples when ungrouped", {
  metaFile <- data.frame(
    Sample = c("S1", "S2", "S3"),
    Group = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )
  resolved <- MOCHA:::.resolveSubSamples(metaFile)

  expect_equal(resolved$subGroups, "All")
  expect_equal(unlist(resolved$subSamples), metaFile$Sample)
})

test_that(".resolveSubSamples respects groupColumn", {
  metaFile <- data.frame(
    Sample = c("S1", "S2", "S3"),
    Group = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )
  resolved <- MOCHA:::.resolveSubSamples(metaFile, groupColumn = "Group")

  expect_equal(sort(resolved$subGroups), c("A", "B"))
  expect_equal(unlist(resolved$subSamples[["A"]]), c("S1", "S2"))
  expect_equal(unlist(resolved$subSamples[["B"]]), "S3")
})

test_that(".resolveSubSamples respects subGroups within groupColumn", {
  metaFile <- data.frame(
    Sample = c("S1", "S2", "S3"),
    Group = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )
  resolved <- MOCHA:::.resolveSubSamples(
    metaFile,
    groupColumn = "Group",
    subGroups = "A"
  )

  expect_equal(resolved$subGroups, "A")
  expect_equal(unlist(resolved$subSamples[["A"]]), c("S1", "S2"))
})

test_that(".requireCellPopulationAssays errors on COUNTS-only assays", {
  expect_error(
    MOCHA:::.requireCellPopulationAssays("COUNTS"),
    "The only assay in the SummarizedExperiment is Counts"
  )
  expect_invisible(MOCHA:::.requireCellPopulationAssays(c("CD4", "CD8")))
})

test_that(".selectCoverageFromBundle handles named and legacy bundles", {
  legacy <- list(sample1 = "legacy")
  named <- list(Accessibility = list(sample1 = "acc"), Insertions = list(sample1 = "ins"))

  expect_equal(
    MOCHA:::.selectCoverageFromBundle(legacy, coverage = TRUE),
    legacy
  )
  expect_equal(
    MOCHA:::.selectCoverageFromBundle(named, coverage = TRUE),
    named$Accessibility
  )
  expect_equal(
    MOCHA:::.selectCoverageFromBundle(named, coverage = FALSE),
    named$Insertions
  )
})

test_that(".normalizeTilePairs accepts character and numeric indices", {
  accMat <- matrix(
    1:6,
    nrow = 3,
    dimnames = list(
      c("chr1:1-500", "chr1:501-1000", "chr1:1001-1500"),
      c("S1", "S2")
    )
  )
  fullObj <- SummarizedExperiment::SummarizedExperiment(list(counts = accMat))

  charPairs <- MOCHA:::.normalizeTilePairs(fullObj, c("chr1:1-500"), c("chr1:501-1000"))
  expect_equal(charPairs$tile1, "chr1:1-500")
  expect_equal(charPairs$tile2, "chr1:501-1000")

  numPairs <- MOCHA:::.normalizeTilePairs(fullObj, 1L, 2L)
  expect_equal(numPairs$tile1, "chr1:1-500")
  expect_equal(numPairs$tile2, "chr1:501-1000")
})

test_that(".computeCoAccessibilityPValues assigns directional p-values", {
  foreGround <- data.frame(
    Correlation = c(0.5, -0.5, NA),
    Tile1 = c("a", "b", "c"),
    Tile2 = c("d", "e", "f")
  )
  backGround <- data.frame(Correlation = c(0, 0.1, 0.2, -0.1, -0.2))

  out <- MOCHA:::.computeCoAccessibilityPValues(foreGround, backGround, verbose = FALSE)
  expect_true("pValues" %in% colnames(out))
  expect_equal(length(out$pValues), nrow(foreGround))
  expect_true(is.na(out$pValues[3]))
})
