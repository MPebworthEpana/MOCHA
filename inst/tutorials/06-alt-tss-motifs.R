# Alternative TSS and TF regulation tutorial — executable chunks.

# ---- libraries ----
library(MOCHA)

if (!exists("differentials")) {
  fixture_dir <- file.path(tempdir(), "mocha_tutorial_fixtures")
  if (!file.exists(file.path(fixture_dir, "differentials.rds"))) {
    stop("Source 00-fixtures.R before 06-alt-tss-motifs.R")
  }
  sampleTileMatrices <- readRDS(file.path(fixture_dir, "sampleTileMatrices.rds"))
  differentials <- readRDS(file.path(fixture_dir, "differentials.rds"))
}

# ---- add-motif-set ----
if (requireNamespace("chromVARmotifs", quietly = TRUE) &&
    requireNamespace("motifmatchr", quietly = TRUE)) {
  sampleTileMatrices <- addMotifSet(
    SampleTileObj = sampleTileMatrices,
    motifPWMs = chromVARmotifs::human_pwms_v2,
    motifSetName = "CISBP"
  )
}

# ---- motif-enrichment ----
sig <- differentials[!is.na(differentials$FDR) &
                     differentials$FDR <= 0.1 &
                     abs(differentials$Log2FC_C) >= 1]
bg <- differentials[!(differentials %in% sig)]
if (length(sig) > 0L && length(bg) > 0L &&
    "CISBP" %in% names(S4Vectors::metadata(sampleTileMatrices))) {
  motifPosList <- S4Vectors::metadata(sampleTileMatrices)$CISBP
  enr <- MotifEnrichment(
    Group1 = sig,
    Group2 = bg,
    motifPosList = motifPosList
  )
  head(enr[order(enr$adjp_val), ])
}

# ---- motifset-enrichment ----
# upstream <- MotifSetEnrichmentAnalysis(
#   ligandTFMatrix = ligandTFMatrix,
#   motifEnrichmentDF = enr
# )

# ---- get-alt-tss ----
altTSS <- getAltTSS(
  completeDAPs = differentials,
  threshold = 0.2,
  TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
  OrgDb = "org.Hs.eg.db"
)
table(altTSS$type)
head(altTSS[altTSS$type == "ii", ])

# ---- plot-alt-tss ----
# candidate <- "MYD88"
# regionGR <- GenomicRanges::GRanges("chr3:38179000-38186000")
# countSE <- getCoverage(
#   sampleTileMatrices,
#   cellPopulations = "C2",
#   regions = regionGR,
#   groupColumn = "Sample"
# )
# plotRegion(countSE = countSE, whichGene = candidate)

# ---- motif-footprint ----
# fp <- motifFootprint(
#   SampleTileObj = sampleTileMatrices,
#   motifName = "CISBP",
#   specMotif = "MA0080.4_SPI1",
#   cellPopulations = "C2",
#   windowSize = 500,
#   normTn5 = TRUE,
#   smoothTn5 = 10,
#   groupColumn = "Sample"
# )

# ---- plot-motifs ----
# fp_stats <- plotMotifs(
#   fp,
#   footprint = "MA0080.4_SPI1",
#   groupColumn = "Sample",
#   returnDF = TRUE,
#   plotIndividualRegions = FALSE
# )

# ---- session-info ----
sessionInfo()
