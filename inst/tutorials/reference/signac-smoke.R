# Tier C: Signac ingest smoke test wrapper.
repo_root <- normalizePath(
  file.path(Sys.getenv("MOCHA_REPO_ROOT", getwd()), "."),
  mustWork = FALSE
)
smoke <- file.path(repo_root, "tests", "scripts", "smoke_seurat_ingest.R")
if (!file.exists(smoke)) {
  smoke <- normalizePath(
    file.path(dirname(getwd()), "tests", "scripts", "smoke_seurat_ingest.R"),
    mustWork = FALSE
  )
}
if (!file.exists(smoke)) {
  stop("Could not locate tests/scripts/smoke_seurat_ingest.R")
}
source(smoke, local = FALSE)
