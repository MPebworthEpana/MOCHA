# MOCHA workflow tutorial — executable chunks (source of truth for vignette code).

# ---- libraries ----
library(MOCHA)
library(MultiAssayExperiment)
library(SummarizedExperiment)

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

# ---- call-open-tiles ----
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

# ---- get-open-tiles ----
openC2 <- getOpenTiles(tileResults, cellPopulations = "C2")
length(openC2[[1]])

# ---- plot-consensus ----
consensusDF <- plotConsensus(
  tileResults,
  cellPopulations = c("C2", "C5"),
  returnDFs = TRUE,
  numCores = 1
)
head(consensusDF$C2)

# ---- suggest-threshold ----
suggested <- suggestConsensusThreshold(
  tileResults,
  cellPopulations = c("C2", "C5"),
  method = "kneedle",
  numCores = 1
)
suggested
if (nrow(suggested) == 0L) {
  message(
    "Bundled PBMCSmall has one sample per population; ",
    "threshold suggestion needs multi-sample curves. Use reproducibilityThreshold = 0 below."
  )
}

# ---- sample-tile-matrix ----
sampleTileMatrices <- getSampleTileMatrix(
  tileResults,
  cellPopulations = c("C2", "C5"),
  reproducibilityThreshold = 0
)
sampleTileMatrices

# ---- sample-tile-matrix-tuned ----
# Example using a non-zero reproducibility threshold from suggestConsensusThreshold():
# thr <- suggested$Reproducibility[suggested$CellPopulation == "C2"]
# sampleTileMatrices <- getSampleTileMatrix(
#   tileResults,
#   cellPopulations = c("C2", "C5"),
#   reproducibilityThreshold = thr,
#   numCores = 1
# )

# ---- subset ----
c2Matrix <- subsetMOCHAObject(
  sampleTileMatrices,
  subsetBy = "celltype",
  groupList = "C2"
)
dim(c2Matrix)

# ---- plot-intensity ----
plotIntensityDistribution(c2Matrix, cellPopulation = "C2", returnDF = TRUE)

# ---- annotate ----
annotated <- annotateTiles(c2Matrix)
head(SummarizedExperiment::rowData(annotated))

# ---- add-motif-set ----
# annotatedWithMotifs <- addMotifSet(
#   SampleTileObj = annotated,
#   motifPWMs = chromVARmotifs::human_pwms_v2,
#   motifSetName = "CISBP"
# )

# ---- differentials ----
samples <- unique(SummarizedExperiment::colData(c2Matrix)$Sample)
if (length(samples) >= 2L) {
  diffs <- getDifferentialAccessibleTiles(
    SampleTileObj = c2Matrix,
    cellPopulations = "C2",
    groupColumn = "Sample",
    foreground = samples[1],
    background = samples[2],
    numCores = 1,
    verbose = TRUE
  )
  head(diffs)
}

# ---- differentials-advanced ----
# diffs_paired <- getDifferentialAccessibleTiles(
#   SampleTileObj = c2Matrix,
#   cellPopulations = "C2",
#   groupColumn = "Sample",
#   foreground = samples[1],
#   background = samples[2],
#   method = "wilcoxon_paired",
#   pairColumn = "Subject",
#   dropoutAdjustment = "biological_only",
#   bioThreshold = 0.8,
#   numCores = 4
# )

# ---- session-info ----
sessionInfo()
