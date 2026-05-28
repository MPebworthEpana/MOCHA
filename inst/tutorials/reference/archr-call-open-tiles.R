# Tier C: ArchR import path (requires local PBMCSmall + ArchR).
# Run with MOCHA_HEAVY_TESTS=true.

if (!requireNamespace("ArchR", quietly = TRUE)) {
  stop("ArchR is required for this reference tutorial script")
}

archr_dir <- Sys.getenv("MOCHA_ARCHR_PBMC_SMALL", unset = "")
if (!nzchar(archr_dir)) {
  candidates <- c(
    file.path("..", "PBMCSmall"),
    file.path("..", "..", "PBMCSmall")
  )
  found <- candidates[dir.exists(candidates)]
  if (length(found) == 0) {
    stop("PBMCSmall ArchR project not found; set MOCHA_ARCHR_PBMC_SMALL")
  }
  archr_dir <- normalizePath(found[1], mustWork = TRUE)
}

library(MOCHA)
library(ArchR)

# ---- archr-setup ----
testProj <- ArchR::loadArchRProject(archr_dir)

# ---- archr-call-open-tiles ----
tileResults <- callOpenTiles(
  ATACFragments = testProj,
  cellPopLabel = "Clusters",
  cellPopulations = c("C2", "C5"),
  TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
  OrgDb = "org.Hs.eg.db",
  outDir = tempdir(),
  numCores = 1,
  verbose = FALSE
)
tileResults
