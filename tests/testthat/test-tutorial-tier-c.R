skip_on_cran()
skip_unless_mocha_heavy()

test_that("Tier C reference tutorial scripts run when fixtures exist", {
  tutorial_dir <- system.file("tutorials", package = "MOCHA")
  if (!nzchar(tutorial_dir)) {
    tutorial_dir <- normalizePath(
      file.path(testthat::test_path(), "..", "..", "inst", "tutorials"),
      mustWork = TRUE
    )
  }
  archr_script <- file.path(tutorial_dir, "reference", "archr-call-open-tiles.R")
  if (requireNamespace("ArchR", quietly = TRUE)) {
    archr_dir <- mocha_archr_project_dir("PBMCSmall")
    if (!is.na(archr_dir)) {
      env <- new.env(parent = globalenv())
      expect_error(source(archr_script, local = env), NA)
    } else {
      skip("PBMCSmall ArchR project not found")
    }
  } else {
    skip("ArchR not installed")
  }
})
