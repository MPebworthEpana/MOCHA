test_that("extractRegion errors when coverage files aren't saved locally", {
  capture.output(
    SampleTileMatrix <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResultsMultisample,
      cellPopulations = "all",
      reproducibilityThreshold = 0,
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

test_that("extractRegion works on multisample data when coverage exists", {
  skip_on_cran()
  # Region under approxLimit (100kb) so extractRegion uses base-pair coverage path.
  region <- "chr1:18137866-18237865"
  sample_ids <- c("scATAC_BMMC_R1", "scATAC_CD34_BMMC_R1")

  capture.output(
    SampleTileMatrix <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResultsMultisample,
      cellPopulations = "all",
      reproducibilityThreshold = 0,
      numCores = 1
    )
  )

  SampleTileMatrix <- SampleTileMatrix[
    , SampleTileMatrix$Sample %in% sample_ids
  ]
  SampleTileMatrix@metadata$Directory <- make_synthetic_coverage_dir(
    cell_populations = "C3",
    sample_ids = sample_ids,
    region = region
  )

  capture.output(
    countSE <- MOCHA::extractRegion(
      SampleTileObj = SampleTileMatrix,
      cellPopulations = "C3",
      region = region,
      numCores = 1,
      sampleSpecific = FALSE
    )
  )

  expect_s4_class(countSE, "SummarizedExperiment")
  expect_gt(SummarizedExperiment::nrow(countSE), 0L)
})

test_that("extractRegion errors when there is no fragment coverage for a cell population", {
  skip_on_cran()
  region <- "chr1:18137866-18237865"

  capture.output(
    SampleTileMatrix <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResultsMultisample,
      cellPopulations = "all",
      reproducibilityThreshold = 0,
      numCores = 1
    )
  )

  SampleTileMatrix@metadata$Directory <- make_synthetic_coverage_dir(
    cell_populations = "C3",
    sample_ids = c("dummy_no_coverage"),
    region = region
  )
  expect_error(
    MOCHA::extractRegion(
      SampleTileObj = SampleTileMatrix,
      cellPopulations = "C3",
      region = region,
      numCores = 1,
      sampleSpecific = FALSE,
      skipEmpty = FALSE
    ),
    "There is no fragment coverage for cell population"
  )
})

test_that("extractRegion works with HemeTutorial coverage (heavy)", {
  skip_unless_mocha_heavy()
  skip_on_cran()
  heme_dir <- skip_unless_heme_coverage()

  region <- "chr1:18137866-18237865"
  sample_ids <- c("scATAC_BMMC_R1", "scATAC_CD34_BMMC_R1")

  capture.output(
    SampleTileMatrix <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResultsMultisample,
      cellPopulations = "all",
      reproducibilityThreshold = 0,
      numCores = 1
    )
  )

  SampleTileMatrix <- SampleTileMatrix[
    , SampleTileMatrix$Sample %in% sample_ids
  ]
  SampleTileMatrix@metadata$Directory <- heme_dir

  capture.output(
    countSE <- MOCHA::extractRegion(
      SampleTileObj = SampleTileMatrix,
      cellPopulations = "C3",
      region = region,
      numCores = 1,
      sampleSpecific = FALSE
    ),
    type = "message"
  )

  expect_s4_class(countSE, "SummarizedExperiment")
  expect_gt(SummarizedExperiment::nrow(countSE), 0L)
})
