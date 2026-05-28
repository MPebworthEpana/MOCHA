#!/usr/bin/env Rscript
# Run the full testthat suite serially (avoids parallel makeCluster port clashes).
Sys.setenv(TESTTHAT_PARALLEL = "false")
if (Sys.getenv("NOT_CRAN", unset = "") == "") {
  Sys.setenv(NOT_CRAN = "true")
}

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("pkgload is required: install.packages('pkgload')")
}
if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("testthat is required: install.packages('testthat')")
}

pkgload::load_all()
testthat::set_max_fails(Inf)
invisible(testthat::test_local())
