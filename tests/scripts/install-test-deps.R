# Install MOCHA test dependencies not available via conda (optional suggests).
pkgs_cran <- c(
  "wCorr", "mixtools", "zip", "uwot", "lifecycle", "tidyselect", "rlang"
)
for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# chromVAR is optional for several tests
optional_bioc <- c("chromVAR")
for (p in optional_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) {
    tryCatch(
      BiocManager::install(p, ask = FALSE, update = FALSE),
      error = function(e) message("Optional package ", p, " not installed: ", conditionMessage(e))
    )
  }
}

message("Dependency install script finished.")
