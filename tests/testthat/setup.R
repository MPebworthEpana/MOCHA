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
