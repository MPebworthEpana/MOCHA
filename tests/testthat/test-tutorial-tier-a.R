skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg19")
skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
skip_if_not_installed("org.Hs.eg.db")

test_that("Tier A tutorial scripts run", {
  tutorial_dir <- system.file("tutorials", package = "MOCHA")
  if (!nzchar(tutorial_dir)) {
    tutorial_dir <- normalizePath(
      file.path(testthat::test_path(), "..", "..", "inst", "tutorials"),
      mustWork = TRUE
    )
  }
  shared <- file.path(tutorial_dir, "_shared.R")
  if (file.exists(shared)) {
    source(shared, local = FALSE)
  }
  for (f in c("01-workflow.R", "02-import-bundled.R")) {
    env <- new.env(parent = globalenv())
    expect_error(
      source(file.path(tutorial_dir, f), local = env),
      NA
    )
  }
})
