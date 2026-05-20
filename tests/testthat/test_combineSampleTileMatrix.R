skip_if_not_installed("chromVAR")
skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg19")
if (requireNamespace("chromVAR", quietly = TRUE) &&
    requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
  test_that("combineSampleTileMatrix works on a 3-sample dataset", {
    cellPopulations <- c("C2", "C3")
    capture.output(
      STObj <- MOCHA::getSampleTileMatrix(
        MOCHA:::testTileResultsMultisample,
        cellPopulations = cellPopulations,
        reproducibilityThreshold = 0
      )
    )

    combinedObj <- MOCHA::combineSampleTileMatrix(STObj)

    coldata <- SummarizedExperiment::colData(combinedObj)
    expect_false(anyNA(coldata$CellCounts))
    expect_false(anyNA(coldata$FragmentCounts))

    summarizedData <- S4Vectors::metadata(STObj)$summarizedData
    cellCountsWide <- as.data.frame(
      SummarizedExperiment::assays(summarizedData)[["CellCounts"]]
    )
    fragCountsWide <- as.data.frame(
      SummarizedExperiment::assays(summarizedData)[["FragmentCounts"]]
    )
    for (i in seq_len(nrow(coldata))) {
      sampleKey <- coldata$Sample[i]
      cellType <- coldata$CellType[i]
      bioSample <- sub(paste0("^", cellType, "__"), "", sampleKey)
      expect_equal(coldata$CellCounts[i], cellCountsWide[cellType, bioSample])
      expect_equal(coldata$FragmentCounts[i], fragCountsWide[cellType, bioSample])
    }

    expect_snapshot_output(
      combinedObj,
      variant = "3sample"
    )
  })

  test_that("combineSampleTileMatrix works on a 1-sample dataset", {
    cellPopulations <- c("C2", "C5")
    capture.output(
      STObj <- MOCHA::getSampleTileMatrix(
        MOCHA:::testTileResults,
        cellPopulations = cellPopulations,
        reproducibilityThreshold = 0
      )
    )

    combinedObj <- MOCHA::combineSampleTileMatrix(STObj)

    coldata <- SummarizedExperiment::colData(combinedObj)
    expect_false(anyNA(coldata$CellCounts))
    expect_false(anyNA(coldata$FragmentCounts))

    summarizedData <- S4Vectors::metadata(STObj)$summarizedData
    cellCountsWide <- as.data.frame(
      SummarizedExperiment::assays(summarizedData)[["CellCounts"]]
    )
    for (i in seq_len(nrow(coldata))) {
      sampleKey <- coldata$Sample[i]
      cellType <- coldata$CellType[i]
      bioSample <- sub(paste0("^", cellType, "__"), "", sampleKey)
      expect_equal(coldata$CellCounts[i], cellCountsWide[cellType, bioSample])
    }

    expect_snapshot_output(
      combinedObj,
      variant = "1sample"
    )
  })
}
