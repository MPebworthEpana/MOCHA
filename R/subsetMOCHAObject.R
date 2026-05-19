#' @title Subset a tileResults object by metadata
#'
#' @description \code{subsetMOCHAObject} subsets a tileResults-type object (from
#'   callOpenTiles), or a SummarizedExperiment-type object (from
#'   getSampleTileMatrix), either by cell type or sample metadata.
#'
#' @param Object A MultiAssayExperiment or RangedSummarizedExperiment,
#' @param subsetBy The variable to subset by. Can either be 'celltype', or a
#'   column from the sample metadata (see `colData(Object)`).
#' @param groupList the list of cell type names or sample-associated data that
#'   should be used to subset the Object
#' @param removeNA If TRUE, removes groups in groupList that are NA. If FALSE,
#'   keep groups that are NA.
#' @param subsetPeaks If `subsetBy` = 'celltype', subset the tile set to tiles
#'   only called in those cell types. Default is TRUE.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#'
#' @return Object the input Object, filtered down to either the cell type or
#'   samples desired.
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
#'   subsetTiles <- MOCHA::subsetMOCHAObject(
#'     tiles,
#'     subsetBy = "celltype",
#'     groupList = "C2"
#'   )
#' }
#' }
#'
#' @export
#' @keywords utils
subsetMOCHAObject <- function(Object,
                              subsetBy,
                              groupList,
                              removeNA = TRUE,
                              subsetPeaks = TRUE,
                              verbose = FALSE) {
  .subset_mocha_object(
    Object,
    subsetBy = subsetBy,
    groupList = groupList,
    removeNA = removeNA,
    subsetPeaks = subsetPeaks,
    verbose = verbose
  )
}

#' @title Modify the cell population names in a Sample-Tile Object from
#'   \code{getSampleTileMatrix()}
#'
#' @description \code{renameCellTypes} Allows you to modify the cell type names
#'   for a MOCHA SampleTileObject, from the assay names, GRanges column names,
#'   and summarizedData (within the metadata), all at once.
#'
#' @param MOCHAObject A  RangedSummarizedExperiment,
#' @param oldNames A list of cell type names that you want to change.
#' @param newNames A list of new cell type names to replace the old names with.
#' @return A MOCHA SampleTile object with new cell types.
#'
#' @export
#' @keywords utils
renameCellTypes <- function(MOCHAObject,
                            oldNames,
                            newNames) {
  if (methods::is(MOCHAObject, "SummarizedExperiment")) {
    if (!any(grepl("getSampleTileMatrix", unlist(MOCHAObject@metadata$History)))) {
      stop("MOCHAObject is not an SampleTile object from MOCHA.")
    }

    if (!all(oldNames %in% names(SummarizedExperiment::assays(MOCHAObject)))) {
      stop("Not all of the provided oldNames exist in the current MOCHAObject")
    }

    if (length(oldNames) != length(newNames)) {
      stop("oldNames and newNames are different lengths.")
    }

    # assay names edits
    assayNames <- names(SummarizedExperiment::assays(MOCHAObject))
    assayNames[match(oldNames, assayNames)] <- newNames
    names(SummarizedExperiment::assays(MOCHAObject)) <- assayNames

    # rowRanges edits
    mColData <- GenomicRanges::mcols(SummarizedExperiment::rowRanges(MOCHAObject))
    colnames(mColData)[match(oldNames, colnames(mColData))] <- newNames
    GenomicRanges::mcols(SummarizedExperiment::rowRanges(MOCHAObject)) <- mColData

    # summarized cell type metadata edits
    oldSumData <- rownames(MOCHAObject@metadata$summarizedData)
    oldSumData[match(oldNames, oldSumData)] <- newNames

    rownames(MOCHAObject@metadata$summarizedData) <- oldSumData

    return(MOCHAObject)
  } else {
    stop("MOCHAObject is not an SampleTile object from MOCHA.")
  }
}
