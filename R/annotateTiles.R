#' @title Annotate tiles with gene annotations
#'
#' @description \code{annotateTiles} annotates a set of sample-tile matrices
#'   given with gene annotations. Details on TxDb and Org annotation packages
#'   and available annotations can be found at Bioconductor:
#'   https://bioconductor.org/packages/3.15/data/annotation/
#'
#' @param Obj A RangedSummarizedExperiment generated from getSampleTileMatrix,
#'   containing TxDb and Org in the metadata. This may also be a GRanges object.
#' @param TxDb The annotation package for TxDb object for your genome.
#'   Optional, only required if Obj is a GRanges.
#' @param Org The genome-wide annotation for your organism.
#'   Optional, only required if Obj is a GRanges.
#' @param promoterRegion Optional list containing the window size in basepairs
#' defining the promoter region. The format is (upstream, downstream).
#' Default is (2000, 100).
#'
#' @return Obj, the input data structure with added gene annotations (whether GRanges or SampleTileObj)
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
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
#'     cellPopulations = "C2",
#'     numCores = 1
#'   )
#'   stm <- MOCHA::getSampleTileMatrix(tiles, cellPopulations = "C2", threshold = 0)
#'   annotated <- MOCHA::annotateTiles(stm)
#' }
#' }
#'
#' @export
#' @keywords utils
annotateTiles <- function(Obj,
                          TxDb = NULL,
                          Org = NULL,
                          promoterRegion = c(2000, 100)) {
  . <- Type <- NULL
  objMeta <- S4Vectors::metadata(Obj)
  if (methods::is(Obj, "RangedSummarizedExperiment") & is.null(TxDb) & is.null(Org)) {
    if (!all(c("TxDb", "OrgDb") %in% names(objMeta))) {
      stop("Error: SampleTileObj as a RangedSummarizedExperiment does not contain a TxDb and/or OrgDb in the metadata. Please provide these as input.")
    }
    tileGRanges <- SummarizedExperiment::rowRanges(Obj)
    TxDb <- getAnnotationDbFromInstalledPkgname(objMeta$TxDb$pkgname, "TxDb")
    Org <- getAnnotationDbFromInstalledPkgname(objMeta$OrgDb$pkgname, "OrgDb")
  } else if (methods::is(Obj, "GRanges") & !is.null(TxDb) & !is.null(Org)) {
    tileGRanges <- Obj
  } else if (methods::is(Obj, "RangedSummarizedExperiment")) {
    tileGRanges <- SummarizedExperiment::rowRanges(Obj)
  } else {
    stop("Error: Invalid inputs. Verify Obj is a RangedSummarizedExperiment and tiles were called from an ArchR project. If Obj is a GRanges, TxDb and Org must be provided.")
  }

  txList <- suppressWarnings(GenomicFeatures::transcriptsBy(TxDb, by = ("gene")))
  names(txList) <- .map_gene_ids(Org, names(txList), column = "SYMBOL", keytype = "ENTREZID")

  txs <- IRanges::stack(txList) %>%
    GenomicRanges::trim() %>%
    S4Vectors::unique(.)

  promoterSet <- IRanges::stack(txList) %>%
    GenomicRanges::trim(.) %>%
    S4Vectors::unique(.) %>%
    GenomicRanges::promoters(., upstream = promoterRegion[1], downstream = promoterRegion[2])


  getOverlapNameList <- function(rowTiles, annotGR) {
    overlapGroup <- IRanges::findOverlaps(rowTiles, annotGR) %>% as.data.frame()
    overlapGroup$Genes <- as.character(annotGR$name[overlapGroup$subjectHits])
    last <- overlapGroup %>%
      dplyr::group_by(.data$queryHits) %>%
      dplyr::summarize(Genes = paste(unique(.data$Genes), collapse = ", "))
    return(last)
  }

  txs_overlaps <- getOverlapNameList(tileGRanges, txs)
  promo_overlaps <- getOverlapNameList(tileGRanges, promoterSet)

  tileType <- as.data.frame(GenomicRanges::mcols(tileGRanges)) %>%
    dplyr::mutate(Index = 1:nrow(.)) %>%
    dplyr::mutate(Type = dplyr::case_when(
      Index %in% promo_overlaps$queryHits ~ "Promoter",
      Index %in% txs_overlaps$queryHits ~ "Intragenic",
      TRUE ~ "Distal"
    )) %>%
    dplyr::left_join(promo_overlaps, by = c("Index" = "queryHits")) %>%
    dplyr::rename("Promo" = .data$Genes) %>%
    dplyr::left_join(txs_overlaps, by = c("Index" = "queryHits")) %>%
    dplyr::rename("Txs" = .data$Genes) %>%
    dplyr::mutate(Genes = ifelse(Type == "Promoter", .data$Promo, NA)) %>%
    dplyr::mutate(Genes = ifelse(Type == "Intragenic", .data$Txs, .data$Genes))

  tileGRanges$tileType <- tileType$Type
  tileGRanges$Gene <- tileType$Genes

  # If input was as Ranged SE, then edit the rowRanges for the SE and return it.
  # Else, return the annotated tile GRanges SampleTileObject.
  if (methods::is(Obj, "RangedSummarizedExperiment")) {
    SummarizedExperiment::rowRanges(Obj) <- tileGRanges
    return(Obj)
  } else {
    return(tileGRanges)
  }
}

#' @title Extract the list of promoter genes from a GRanges annotated with
#'   \code{annotateTiles()}
#'
#' @description \code{getPromoterGenes} Takes rowRanges from annotateTiles and
#'   extracts a unique list of genes.
#'
#' @param GRangesObj a GRanges object with a metadata column for tileType and
#'   Gene.
#' @return vector of strings with gene names.
#'
#' @export
#' @keywords utils
getPromoterGenes <- function(GRangesObj) {
  tileType <- NULL
  if (!methods::is(GRangesObj, "GRanges")) {
    stop("Object provided is not a GRanges object.")
  }

  if (!all(c("tileType", "Gene") %in% colnames(GenomicRanges::mcols(GRangesObj)))) {
    stop("GRanges object does not contain tileType or Gene column. Run annotateTiles on this GRanges object and try again.")
  }

  promoterTiles <- dplyr::filter(GRangesObj, tileType == "Promoter")
  geneString <- paste0(unlist(GenomicRanges::mcols(promoterTiles)$Gene), collapse = ", ")
  genes <- unique(unlist(stringr::str_split(geneString, pattern = ", "), recursive = TRUE))

  return(genes)
}
