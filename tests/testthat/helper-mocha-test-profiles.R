#' Centralized MOCHA test profile helpers.
#'
#' Default CI runs deterministic unit tests only. Set MOCHA_HEAVY_TESTS=true
#' (or MOCHA_HEAVY_TESTS=1) to enable integration tests that need large local
#' ArchR projects or long runtimes.

mocha_heavy_tests_enabled <- function() {
  val <- Sys.getenv("MOCHA_HEAVY_TESTS", unset = "false")
  isTRUE(tolower(val) %in% c("1", "true", "yes"))
}

skip_unless_mocha_heavy <- function() {
  if (!mocha_heavy_tests_enabled()) {
    testthat::skip("Set MOCHA_HEAVY_TESTS=true to run heavy integration tests")
  }
}

skip_unless_cutpointr <- function() {
  if (!requireNamespace("cutpointr", quietly = TRUE)) {
    testthat::skip(
      "cutpointr not installed (Suggests); needed for threshMethod = 'youden'"
    )
  }
}

mocha_thresh_method_for_tests <- function() {
  if (requireNamespace("cutpointr", quietly = TRUE)) {
    "youden"
  } else {
    "f1"
  }
}

mocha_testthat_dir <- function() {
  normalizePath(testthat::test_path(), winslash = "/", mustWork = TRUE)
}

mocha_repo_root <- function() {
  normalizePath(file.path(mocha_testthat_dir(), "..", ".."), winslash = "/", mustWork = TRUE)
}

mocha_fixture_path <- function(...) {
  file.path(mocha_testthat_dir(), "fixtures", ...)
}

mocha_external_path <- function(...) {
  file.path(mocha_repo_root(), ...)
}

mocha_archr_project_dir <- function(name = c("PBMCSmall", "FullCovid", "HemeTutorial")) {
  name <- match.arg(name)
  candidates <- c(
    file.path(mocha_repo_root(), name),
    file.path(mocha_repo_root(), "..", name),
    file.path(mocha_testthat_dir(), "..", "..", "..", name)
  )
  found <- candidates[dir.exists(candidates)]
  if (length(found) == 0) {
    return(NA_character_)
  }
  normalizePath(found[[1]], winslash = "/", mustWork = TRUE)
}

mocha_heme_mocha_dir <- function() {
  candidates <- c(
    mocha_external_path("HemeTutorial", "MOCHA"),
    file.path(mocha_testthat_dir(), "..", "..", "..", "HemeTutorial", "MOCHA")
  )
  found <- candidates[dir.exists(candidates)]
  if (length(found) == 0) {
    return(NA_character_)
  }
  normalizePath(found[[1]], winslash = "/", mustWork = TRUE)
}

skip_unless_heme_coverage <- function() {
  dir <- mocha_heme_mocha_dir()
  if (is.na(dir)) {
    testthat::skip("HemeTutorial/MOCHA coverage directory not found (optional heavy fixture)")
  }
  invisible(dir)
}

skip_unless_archr_project <- function(name = "PBMCSmall") {
  dir <- mocha_archr_project_dir(name)
  if (is.na(dir)) {
    testthat::skip(sprintf("ArchR project '%s' not found (optional heavy fixture)", name))
  }
  if (!requireNamespace("ArchR", quietly = TRUE)) {
    testthat::skip("ArchR not installed")
  }
  invisible(dir)
}

#' Minimal SampleTileMatrix for differential-accessibility unit tests.
make_synthetic_sample_tile_matrix <- function(n_tiles = 40L, n_per_group = 4L) {
  n_samples <- n_per_group * 2L
  tile_names <- paste0(
    "chr1:",
    seq(10000, by = 500, length.out = n_tiles),
    "-",
    seq(10499, by = 500, length.out = n_tiles)
  )
  set.seed(42)
  counts <- matrix(
    stats::runif(n_tiles * n_samples, min = 8, max = 16),
    nrow = n_tiles,
    dimnames = list(tile_names, paste0("sample_", seq_len(n_samples)))
  )
  gr <- MOCHA::StringsToGRanges(tile_names)
  S4Vectors::mcols(gr)$C2 <- TRUE
  cd <- data.frame(
    Sample = colnames(counts),
    PassQC = TRUE,
    GroupA = rep(c("A", "B"), each = n_per_group),
    stringsAsFactors = FALSE
  )
  SummarizedExperiment::SummarizedExperiment(
    assays = list(C2 = counts),
    rowRanges = gr,
    colData = cd
  )
}
