#!/usr/bin/env Rscript
# Run subclass / mocha subset tests (used from Bioconductor Docker in WSL).
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages(c("remotes", "BiocManager"), repos = "https://cloud.r-project.org")
}
repo_root <- Sys.getenv("MOCHA_ROOT", unset = "/MOCHA")
remotes::install_deps(
  repo_root,
  dependencies = TRUE,
  suggests = FALSE,
  upgrade = FALSE
)

if (!dir.exists(repo_root)) {
  stop("Set MOCHA_ROOT to the package source directory (e.g. /MOCHA in Docker).")
}
if (!requireNamespace("pkgload", quietly = TRUE)) {
  install.packages("pkgload", repos = "https://cloud.r-project.org")
}
message("Loading MOCHA from: ", repo_root)
pkgload::load_all(repo_root, export_all = FALSE)

res <- testthat::test_local(
  filter = "subclass|mocha_",
  stop_on_failure = FALSE
)
df <- as.data.frame(res)
message("Tests: ", nrow(df), " | failed: ", sum(df$failed), " | error: ", sum(df$error))
if (sum(df$failed) > 0 || sum(df$error) > 0) {
  print(df[df$failed > 0 | df$error > 0, , drop = FALSE])
  quit(status = 1)
}
quit(status = 0)
