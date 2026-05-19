test_that("getSampleCellTypeMetadata reads counts from summarizedData", {
  metaSE <- MOCHA::getSampleCellTypeMetadata(MOCHA:::testTileResults)

  expect_s4_class(metaSE, "SummarizedExperiment")
  expect_true(all(c("CellCounts", "FragmentCounts") %in% SummarizedExperiment::assayNames(metaSE)))
  expect_equal(rownames(SummarizedExperiment::colData(metaSE)), rownames(MOCHA:::testTileResults@colData))
  expect_equal(
    SummarizedExperiment::assay(metaSE, "CellCounts"),
    as.matrix(SummarizedExperiment::assays(MOCHA:::testTileResults@metadata$summarizedData)[["CellCounts"]])
  )
})

test_that("getSampleCellTypeMetadata works on getSampleTileMatrix output", {
  capture.output(
    STObj <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResults,
      cellPopulations = c("C2", "C5"),
      threshold = 0,
      numCores = 1
    ),
    type = "message"
  )

  metaSE <- MOCHA::getSampleCellTypeMetadata(STObj)
  expect_s4_class(metaSE, "SummarizedExperiment")
  expect_equal(nrow(SummarizedExperiment::colData(metaSE)), ncol(STObj))
})

test_that("getSampleCellTypeMetadata supports legacy top-level metadata", {
  STObj <- MOCHA::getSampleTileMatrix(
    MOCHA:::testTileResults,
    cellPopulations = c("C2", "C5"),
    threshold = 0,
    numCores = 1
  )

  summarizedData <- S4Vectors::metadata(STObj)$summarizedData
  legacyObj <- STObj
  legacyObj@metadata$summarizedData <- NULL
  legacyObj@metadata$CellCounts <- as.data.frame(
    SummarizedExperiment::assays(summarizedData)[["CellCounts"]]
  )
  legacyObj@metadata$FragmentCounts <- as.data.frame(
    SummarizedExperiment::assays(summarizedData)[["FragmentCounts"]]
  )

  metaSE <- MOCHA::getSampleCellTypeMetadata(legacyObj)
  expect_equal(
    SummarizedExperiment::assay(metaSE, "CellCounts"),
    SummarizedExperiment::assay(MOCHA::getSampleCellTypeMetadata(STObj), "CellCounts")
  )
})

test_that("getSampleCellTypeMetadata errors when counts are missing", {
  STObj <- MOCHA::getSampleTileMatrix(
    MOCHA:::testTileResults,
    cellPopulations = c("C2", "C5"),
    threshold = 0,
    numCores = 1
  )
  STObj@metadata$summarizedData <- NULL

  expect_error(
    MOCHA::getSampleCellTypeMetadata(STObj),
    regexp = "does not contain Sample-Celltype metadata"
  )
})
