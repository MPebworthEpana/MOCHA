test_that("extractRegion errors when coverage files aren't saved locally", {
  capture.output(
    SampleTileMatrix <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResultsMultisample,
      cellPopulations = "all",
      threshold = 0,
      numCores = 1
    )
  )

  SampleTileMatrix@metadata$Directory <- file.path(tempdir(), "idontexist")
  expect_error(
    MOCHA::extractRegion(
      SampleTileObj = SampleTileMatrix,
      cellPopulations = "ALL",
      region = "chr1:18137866-38139912",
      numCores = 1,
      sampleSpecific = FALSE
    ),
    "does not exist"
  )
})

test_that("extractRegion works on multisample data when HemeTutorial coverage exists", {
  skip_on_cran()
  heme_dir <- skip_unless_heme_coverage()

  capture.output(
    SampleTileMatrix <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResultsMultisample,
      cellPopulations = "all",
      threshold = 0,
      numCores = 1
    )
  )

  SampleTileMatrix <- SampleTileMatrix[
    , SampleTileMatrix$Sample %in% c("scATAC_BMMC_R1", "scATAC_CD34_BMMC_R1")
  ]
  SampleTileMatrix@metadata$Directory <- heme_dir

  capture.output(
    countSE <- MOCHA::extractRegion(
      SampleTileObj = SampleTileMatrix,
      cellPopulations = "C3",
      region = "chr1:18137866-38139912",
      numCores = 1,
      sampleSpecific = FALSE
    )
  )

  expect_s4_class(countSE, "RangedSummarizedExperiment")
  expect_gt(SummarizedExperiment::nrow(countSE), 0L)
})

test_that("extractRegion errors when there is no fragment coverage for a cell population", {
  skip_on_cran()
  heme_dir <- skip_unless_heme_coverage()

  capture.output(
    SampleTileMatrix <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResultsMultisample,
      cellPopulations = "all",
      threshold = 0,
      numCores = 1
    )
  )

  SampleTileMatrix@metadata$Directory <- heme_dir
  expect_error(
    MOCHA::extractRegion(
      SampleTileObj = SampleTileMatrix,
      cellPopulations = "C3",
      region = "chr1:18137866-38139912",
      numCores = 1,
      sampleSpecific = FALSE
    ),
    "There is no fragment coverage for cell population"
  )
})
