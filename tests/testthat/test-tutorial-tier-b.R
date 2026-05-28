skip_on_cran()

test_that("Tier B tutorial scripts run after fixtures", {
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
  env <- new.env(parent = globalenv())
  source(file.path(tutorial_dir, "00-fixtures.R"), local = env)
  for (f in c(
    "03-downstream.R",
    "04-export.R",
    "05-advanced-modeling.R",
    "06-alt-tss-motifs.R"
  )) {
    child <- new.env(parent = env)
    expect_error(source(file.path(tutorial_dir, f), local = child), NA)
  }
})
