skip_on_cran()

test_that("MOCHAPeakModel validation rejects invalid objects", {
  expect_error(
    MOCHA:::.validatePeakModel(list()),
    "MOCHAPeakModel"
  )
  default <- MOCHA:::.defaultPeakModel()
  expect_s3_class(default, "MOCHAPeakModel")
  expect_invisible(MOCHA:::.validatePeakModel(default))
})

test_that("default peak model matches bundled coefficients", {
  default <- MOCHA:::.defaultPeakModel()
  expect_equal(default$tileSize, MOCHA:::.MOCHA_DEFAULT_TILE_SIZE)
  expect_equal(default$trainingMedian, 3668)
  expect_identical(default$finalModelObject, MOCHA::finalModelObject)
  expect_identical(default$youden_threshold, MOCHA::youden_threshold)
})

test_that("trainPeakModel returns MOCHAPeakModel with expected structure", {
  frags <- MOCHA::exampleFragments[[1]]
  peaks <- GenomicRanges::reduce(frags)
  peaks <- peaks[IRanges::width(peaks) >= 200L]

  model <- MOCHA::trainPeakModel(
    ATACFragments = frags,
    cellColData = MOCHA::exampleCellColData,
    blackList = MOCHA::exampleBlackList,
    groundTruthPeaks = peaks,
    tileSize = 250L,
    cellSubsetSizes = c(50L, 100L),
    replicatesFn = function(n) 2L,
    threshMethod = mocha_thresh_method_for_tests(),
    numCores = 1L,
    seed = 42L,
    verbose = FALSE
  )

  expect_s3_class(model, "MOCHAPeakModel")
  expect_equal(model$tileSize, 250L)
  expect_true(model$trainingMedian > 0)
  expect_true(inherits(model$youden_threshold, c("loess", "lm")))
  expect_true(all(c("Loess", "Linear") %in% names(model$finalModelObject)))
  expect_output(print(model), "MOCHA peak-calling model")
})

test_that("trainPeakModel youden threshold uses cutpointr", {
  skip_unless_cutpointr()

  frags <- MOCHA::exampleFragments[[1]]
  peaks <- GenomicRanges::reduce(frags)
  peaks <- peaks[IRanges::width(peaks) >= 200L]

  model <- MOCHA::trainPeakModel(
    ATACFragments = frags,
    cellColData = MOCHA::exampleCellColData,
    blackList = MOCHA::exampleBlackList,
    groundTruthPeaks = peaks,
    tileSize = 250L,
    cellSubsetSizes = c(50L, 100L),
    replicatesFn = function(n) 2L,
    threshMethod = "youden",
    numCores = 1L,
    seed = 7L,
    verbose = FALSE
  )

  expect_s3_class(model, "MOCHAPeakModel")
  expect_equal(model$method$thresh, "youden")
})

test_that("callOpenTiles accepts custom peakModel and uses tileSize", {
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  skip_if_not_installed("org.Hs.eg.db")

  frags <- MOCHA::exampleFragments[[1]]
  peaks <- GenomicRanges::reduce(frags)
  peaks <- peaks[IRanges::width(peaks) >= 200L]

  model <- MOCHA::trainPeakModel(
    ATACFragments = frags,
    cellColData = MOCHA::exampleCellColData,
    blackList = MOCHA::exampleBlackList,
    groundTruthPeaks = peaks,
    tileSize = 250L,
    cellSubsetSizes = c(50L, 100L),
    replicatesFn = function(n) 2L,
    threshMethod = "f1",
    numCores = 1L,
    seed = 1L
  )

  capture.output(
    tiles <- MOCHA::callOpenTiles(
      ATACFragments = MOCHA::exampleFragments,
      cellColData = MOCHA::exampleCellColData,
      blackList = MOCHA::exampleBlackList,
      genome = "hg19",
      TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
      OrgDb = "org.Hs.eg.db",
      outDir = tempdir(),
      cellPopLabel = "Clusters",
      cellPopulations = "C2",
      peakModel = model,
      numCores = 1L
    ),
    type = "message"
  )

  expect_true(methods::is(tiles, "MultiAssayExperiment"))
  exp <- tiles[["C2"]]
  rr <- SummarizedExperiment::rowRanges(exp)
  tileWidths <- IRanges::width(rr[IRanges::width(rr) == 250L])
  expect_gt(length(tileWidths), 0L)
  expect_true(all(tileWidths == 250L))
})

test_that("callOpenTiles without peakModel matches default peak model path", {
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  skip_if_not_installed("org.Hs.eg.db")

  capture.output(
    tiles_default <- MOCHA::callOpenTiles(
      ATACFragments = MOCHA::exampleFragments,
      cellColData = MOCHA::exampleCellColData,
      blackList = MOCHA::exampleBlackList,
      genome = "hg19",
      TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
      OrgDb = "org.Hs.eg.db",
      outDir = tempdir(),
      cellPopLabel = "Clusters",
      cellPopulations = "C2",
      peakModel = NULL,
      numCores = 1L
    ),
    type = "message"
  )

  capture.output(
    tiles_explicit <- MOCHA::callOpenTiles(
      ATACFragments = MOCHA::exampleFragments,
      cellColData = MOCHA::exampleCellColData,
      blackList = MOCHA::exampleBlackList,
      genome = "hg19",
      TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
      OrgDb = "org.Hs.eg.db",
      outDir = tempdir(),
      cellPopLabel = "Clusters",
      cellPopulations = "C2",
      peakModel = MOCHA:::.defaultPeakModel(),
      numCores = 1L
    ),
    type = "message"
  )

  exp1 <- tiles_default[["C2"]]
  exp2 <- tiles_explicit[["C2"]]
  expect_equal(
    SummarizedExperiment::rowRanges(exp1),
    SummarizedExperiment::rowRanges(exp2)
  )
})

test_that("MACS2 wrapper is gated when binary absent", {
  skip_if_not_installed("rtracklayer")
  if (Sys.which("macs2") != "") {
    skip("macs2 is installed; skipping absence test")
  }
  expect_error(
    MOCHA::trainPeakModel(
      ATACFragments = MOCHA::exampleFragments[[1]],
      cellColData = MOCHA::exampleCellColData,
      blackList = MOCHA::exampleBlackList,
      callPeaks = TRUE,
      threshMethod = "f1"
    ),
    "MACS2"
  )
})
