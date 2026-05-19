
#' @title Merge the TSAM from multiple cell populations into a single matrix
#'
#' @description \code{combineSampleTileMatrix} combines all celltypes in a
#'   SampleTileMatrix object into a SummarizedExperiment with one single matrix
#'   across all cell types and samples,
#'
#' @param SampleTileObj The SummarizedExperiment object output from
#'   getSampleTileMatrix containing your sample-tile matrices
#' @param NAtoZero Set NA values in the sample-tile matrix to zero
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param returnClass Return type: \code{"legacy"} (default) or \code{"mocha"}
#'   (\code{MochaSampleTileMatrix} when the input is a MOCHA sample-tile matrix).
#' @return A \code{RangedSummarizedExperiment} with one matrix across cell types.
#'
#' @examples
#' \donttest{
#' if (
#'   requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE) &&
#'     requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE) &&
#'     requireNamespace("org.Hs.eg.db", quietly = TRUE)
#' ) {
#'   tiles <- MOCHA::callOpenTiles(
#'     ATACFragments = MOCHA::exampleFragments,
#'     cellColData = MOCHA::exampleCellColData,
#'     blackList = MOCHA::exampleBlackList,
#'     genome = "BSgenome.Hsapiens.UCSC.hg19",
#'     TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
#'     OrgDb = "org.Hs.eg.db",
#'     outDir = tempdir(),
#'     cellPopLabel = "Clusters",
#'     cellPopulations = c("C2", "C5"),
#'     numCores = 1
#'   )
#'   stm <- MOCHA::getSampleTileMatrix(
#'     tiles,
#'     cellPopulations = c("C2", "C5"),
#'     threshold = 0
#'   )
#'   combined <- MOCHA::combineSampleTileMatrix(stm)
#' }
#' }
#'
#' @export
#' @keywords utils
combineSampleTileMatrix <- function(SampleTileObj,
                                    NAtoZero = TRUE,
                                    verbose = FALSE,
                                    returnClass = c("legacy", "mocha")) {
  returnClass <- match.arg(returnClass)
  CellTypes <- FragNumber <- NULL

  genome <- S4Vectors::metadata(SampleTileObj)$Genome
  genome <- BSgenome::getBSgenome(genome)

  Sample <- Freq <- . <- NULL
  # Extract all the Sample-Tile Matrices for each cell type
  assays <- SummarizedExperiment::assays(SampleTileObj)

  coldata <- SampleTileObj@colData

  # Let's generate a new assay, that will contain the
  # the intensity for a given cell, as well as the
  # median intensity per sample-tile for all other cell types (i.e. the background)
    
  newAssays <- list(do.call("cbind", methods::as(assays, "list")))
  newSamplesNames <- unlist(lapply(names(assays), function(x) {
    gsub(" ", "_", paste(x, colnames(SampleTileObj), sep = "__"))
  }))

  if (length(newSamplesNames) != dim(newAssays[[1]])[2]) {
    stop('columns corrupted in sample tile object')
  }

  names(newAssays) <- "counts"
  colnames(newAssays[[1]]) <- newSamplesNames

  if (NAtoZero) {
    newAssays[[1]][is.na(newAssays[[1]])] <- 0
  }

  # Combine Sample and Cell type for the new columns.
  # This takes the colData and repeats it across cell types for a given sample
  # Rows are now CellType__Sample. So for example 'CD16 Mono' and 'Sample1' becomes "CD16_Mono__Sample1"
  allSampleData <- as.data.frame(do.call("rbind", lapply(names(assays), function(x) {
    tmp_meta <- coldata
    tmp_meta$Sample <- gsub(" ", "_", paste(x, tmp_meta$Sample, sep = "__"))
    tmp_meta$CellType <- rep(x, dim(tmp_meta)[1])
    rownames(tmp_meta) <- tmp_meta$Sample
    tmp_meta
  })))

  # This is where cell counts and fragments counts are pivoted into a long format and merged (left-Joined) into the new allSampleData
  cellTypeLabelList <- Var1 <- NULL

  summarizedData <- S4Vectors::metadata(SampleTileObj)$summarizedData
  .pivot_counts_to_sample_keys <- function(countsWide) {
    countsWide$CellPop <- rownames(countsWide)
    rownames(countsWide) <- NULL
    countsLong <- tidyr::pivot_longer(
      countsWide,
      cols = -CellPop,
      names_to = "bioSample",
      values_to = "countValue"
    )
    countsLong$Sample <- gsub(" ", "_", paste(countsLong$CellPop, countsLong$bioSample, sep = "__"))
    countsLong
  }

  cellCountsWide <- as.data.frame(
    SummarizedExperiment::assays(summarizedData)[["CellCounts"]]
  )
  cellCountsLong <- .pivot_counts_to_sample_keys(cellCountsWide)
  cellCountsLong <- dplyr::rename(cellCountsLong, CellCounts = countValue)

  fragCountsWide <- as.data.frame(
    SummarizedExperiment::assays(summarizedData)[["FragmentCounts"]]
  )
  fragCountsLong <- .pivot_counts_to_sample_keys(fragCountsWide)
  fragCountsLong <- dplyr::rename(fragCountsLong, FragmentCounts = countValue)

  if (anyDuplicated(cellCountsLong$Sample) > 0 || anyDuplicated(fragCountsLong$Sample) > 0) {
    stop("Duplicate CellType__Sample keys found while joining count metadata.")
  }

  missingKeys <- setdiff(allSampleData$Sample, cellCountsLong$Sample)
  if (length(missingKeys) > 0) {
    stop(
      "Could not map all combined sample keys to CellCounts metadata. Missing: ",
      paste(missingKeys, collapse = ", ")
    )
  }

  allSampleData <- dplyr::left_join(
    allSampleData,
    dplyr::select(cellCountsLong, Sample, CellCounts),
    by = "Sample"
  )

  allSampleData <- dplyr::left_join(
    allSampleData,
    dplyr::select(fragCountsLong, Sample, FragmentCounts),
    by = "Sample"
  )

  if (nrow(allSampleData) != length(newSamplesNames)) {
    stop("Count metadata join changed the number of combined samples.")
  }

  # Artificially set all the cell type columns in rowRanges to TRUE, incase of later subsetting.
  allRanges <- SummarizedExperiment::rowRanges(SampleTileObj)

  newMetadata <- S4Vectors::metadata(SampleTileObj)
  newMetadata$History <- append(newMetadata$History, paste("combineSampleTileMatrix", utils::packageVersion("MOCHA")))

  newObj <- SummarizedExperiment::SummarizedExperiment(
    assays = newAssays,
    colData = allSampleData,
    rowRanges = allRanges,
    metadata = newMetadata
  )
  .mocha_promote_return(newObj, returnClass)
}
