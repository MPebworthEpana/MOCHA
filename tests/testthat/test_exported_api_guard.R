# Guardrail: every exported symbol should be referenced in tests or allowlisted.

test_that("exported APIs are referenced in tests or explicitly allowlisted", {
  ns <- readLines(system.file("NAMESPACE", package = "MOCHA"))
  if (length(ns) == 0) {
    ns <- readLines(file.path(mocha_repo_root(), "NAMESPACE"))
  }
  exports <- grep("^export\\(", ns, value = TRUE)
  exports <- gsub("^export\\((.+)\\)$", "\\1", exports)

  allowlist <- c(
    # Optional / heavy integration
    "bulkUMAP",
    "motifFootprint",
    "plotMotifs",
    "plotIntensityDistribution",
    "exportLocalFootprints",
    "testCoAccessibilityChromVar",
    "testCoAccessibilityRandom",
    # Thin wrappers or vignette-only utilities
    "addAccessibilityShift",
    "addInsertionBias",
    "correctGenome",
    "exportDifferentials",
    "getAltTSS",
    "getPromoterGenes",
    "getCellTypeMotifs",
    "getCellTypeTiles",
    "getCellTypes",
    "mergeTileResults",
    "renameCellTypes"
  )

  test_files <- list.files(
    mocha_testthat_dir(),
    pattern = "^test_.*\\.R$",
    full.names = TRUE
  )
  test_text <- paste(vapply(test_files, readLines, character()), collapse = "\n")

  missing <- character()
  for (fn in exports) {
    if (fn %in% allowlist) {
      next
    }
    if (!grepl(paste0("\\b", fn, "\\b"), test_text, fixed = FALSE)) {
      missing <- c(missing, fn)
    }
  }

  expect_equal(
    missing,
    character(),
    info = paste(
      "Add tests or allowlist entries for:",
      paste(missing, collapse = ", ")
    )
  )
})
