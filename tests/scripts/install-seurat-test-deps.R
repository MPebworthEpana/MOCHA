# Install Seurat/Signac test dependencies (not all available via conda).
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_pkgs <- c("Signac", "sparseMatrixStats")
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing ", p, " ...")
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

if (!requireNamespace("Seurat", quietly = TRUE)) {
  message("Installing Seurat from CRAN ...")
  install.packages("Seurat", repos = "https://cloud.r-project.org")
}

message("Done. Verify with:")
message('  Rscript -e \'print(sapply(c("Seurat","Signac"), requireNamespace, quietly=TRUE))\'')
