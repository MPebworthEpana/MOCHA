test_that("core APIs accept a trivial RangedSummarizedExperiment subclass", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("GenomicRanges")

  stm <- make_synthetic_sample_tile_matrix(n_tiles = 10L, n_per_group = 2L)
  stm@metadata <- list(
    summarizedData = SummarizedExperiment::SummarizedExperiment(
      assays = list(
        CellCounts = matrix(1L, nrow = 1L, ncol = ncol(stm)),
        FragmentCounts = matrix(2L, nrow = 1L, ncol = ncol(stm))
      )
    ),
    Genome = "hg19",
    TxDb = list(pkgname = "TxDb.Hsapiens.UCSC.hg38.knownGene"),
    OrgDb = list(pkgname = "org.Hs.eg.db"),
    Directory = tempdir(),
    History = list("getSampleTileMatrix 1.2.0")
  )
  rownames(stm@metadata$summarizedData) <- "C2"

  cls_env <- new.env(parent = emptyenv())
  setClass("TestMOCHARSE", contains = "RangedSummarizedExperiment", where = cls_env)
  on.exit(methods::removeClass("TestMOCHARSE", where = cls_env), add = TRUE)

  sub_stm <- methods::new("TestMOCHARSE", stm)

  expect_s4_class(sub_stm, "RangedSummarizedExperiment")
  expect_equal(MOCHA::getCellTypes(sub_stm), "C2")
  expect_s4_class(MOCHA::getCellTypeTiles(sub_stm, "C2"), "GRanges")

  sub2 <- MOCHA::subsetMOCHAObject(
    sub_stm,
    subsetBy = "celltype",
    groupList = "C2",
    subsetPeaks = FALSE
  )
  expect_s4_class(sub2, "TestMOCHARSE")
})
