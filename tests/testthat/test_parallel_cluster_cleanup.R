skip_if_not_installed("chromVAR")
skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg19")

if (requireNamespace("chromVAR", quietly = TRUE) &&
    requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
  test_that("getSampleTileMatrix stops parallel workers after errors", {
    expect_error(
      capture.output(
        MOCHA::getSampleTileMatrix(
          MOCHA:::testTileResultsMultisample,
          cellPopulations = c("C2", "C3"),
          threshold = 1.1,
          numCores = 1
        ),
        type = "message"
      ),
      regexp = "not a vector"
    )

    capture.output(
      STObj <- MOCHA::getSampleTileMatrix(
        MOCHA:::testTileResultsMultisample,
        cellPopulations = c("C2", "C3"),
        threshold = 0,
        numCores = 1
      ),
      type = "message"
    )

    expect_s4_class(STObj, "RangedSummarizedExperiment")
  })
}
