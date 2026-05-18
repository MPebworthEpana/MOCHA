# Optional ArchR test project cleanup (heavy local setup only).
if (mocha_heavy_tests_enabled() &&
    requireNamespace("ArchR", quietly = TRUE) &&
    dir.exists("PBMCSmall")) {
  withr::local_options(list(timeout = 600))
  capture.output(proj <- ArchR::getTestProject(), type = "message")
  outdir <- dirname(ArchR::getOutputDirectory(proj))
  withr::defer(unlink(ArchR::getOutputDirectory(proj), recursive = TRUE), teardown_env())
  withr::defer(unlink(file.path(outdir, "__MACOSX"), recursive = TRUE), teardown_env())
  withr::defer(unlink(file.path(outdir, "ArchRLogs"), recursive = TRUE), teardown_env())
}

# testthat sets NOT_CRAN during R CMD check; local devtools::test() does not.
if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
  Sys.setenv(NOT_CRAN = "true")
}

# Prefer tabix/bgzip from active conda env (mocha-test) when running tests in WSL.
conda_prefix <- Sys.getenv("CONDA_PREFIX", unset = "")
if (nzchar(conda_prefix)) {
  conda_bin <- file.path(conda_prefix, "bin")
  if (dir.exists(conda_bin)) {
    Sys.setenv(PATH = paste(conda_bin, Sys.getenv("PATH"), sep = .Platform$path.sep))
  }
}
