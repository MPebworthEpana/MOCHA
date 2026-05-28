# Data Import Tutorial — bundled GRangesList section (Bioconductor-safe).

# ---- check-deps ----
has_deps <- all(
  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE),
  requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE),
  requireNamespace("org.Hs.eg.db", quietly = TRUE)
)
if (!has_deps) {
  stop(
    "Install suggested packages: BSgenome.Hsapiens.UCSC.hg19, ",
    "TxDb.Hsapiens.UCSC.hg38.knownGene, org.Hs.eg.db"
  )
}

# ---- bundled-libraries ----
library(MOCHA)

# ---- bundled-call-open-tiles ----
tileResults <- callOpenTiles(
  ATACFragments = exampleFragments,
  cellColData = exampleCellColData,
  blackList = exampleBlackList,
  genome = "BSgenome.Hsapiens.UCSC.hg19",
  TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
  OrgDb = "org.Hs.eg.db",
  outDir = tempdir(),
  cellPopLabel = "Clusters",
  cellPopulations = c("C2", "C5"),
  numCores = 1,
  verbose = TRUE
)
tileResults

# ---- bundled-get-open-tiles ----
openC2 <- getOpenTiles(tileResults, cellPopulations = "C2")
length(openC2[[1]])

# ---- import-utilities ----
# Per-sample fragment GRanges from an ArchR project (requires ArchR):
# frags <- getPopFrags(ArchRProj, cellPopLabel = "Clusters")

# Quick access to called peaks from tileResults
openC2_util <- getOpenTiles(tileResults, cellPopulations = "C2")
# equivalent: getOpenTiles(tileResults, cellPopulations = "C2")

# ---- session-info ----
sessionInfo()
