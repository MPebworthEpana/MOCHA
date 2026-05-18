#' @title Perform peak-calling on a set of fragments or an 
#'   ArchR Project
#'
#' @description \code{callOpenTiles} is the main peak-calling function in MOCHA
#'   that serves as a wrapper function to call peaks provided a set of fragment
#'   files and an ArchR Project for meta-data purposes
#'
#'
#' @param ATACFragments an ArchR Project, a Seurat object with a Signac
#'   ChromatinAssay (fragment paths in the assay), or a GRangesList of fragments.
#'   Two input formats are supported for non-ArchR GRangesList input:
#'   \itemize{
#'     \item **Legacy:** each list element is named `CellPopulation#Sample`
#'       (optional `__normalizationFactor` suffix), with mixed cell types already
#'       split per sample.
#'     \item **Sample-level:** each list element is named by sample ID (matching
#'       `cellColData[[sampleColumn]]`), containing fragments from all cell types
#'       in that sample. Cell populations are derived from `cellColData` using
#'       `cellPopLabel` and `cellCol`.
#'   }
#'   Each GRanges must contain unique cell IDs in the column given by `cellCol`.
#' @param cellPopLabel string indicating which column in the ArchRProject
#'   metadata contains the cell population label.
#' @param cellPopulations vector of strings. Cell subsets for which to call
#'   peaks. This list of group names must be identical to names that appear in
#'   the ArchRProject metadata.  Optional, if cellPopulations='ALL', then peak
#'   calling is done on all cell populations in the ArchR project metadata.
#'   Default is 'ALL'.
#' @param cellColData A DataFrame containing cell-level metadata. This must
#'   contain both a column 'Sample' with unique sample IDs and the column
#'   specified by 'cellPopLabel'.
#' @param sampleColumn The name of the metadata column with biological sample names. Standard is 'Sample', as fixed by ArchR
#' @param blackList A GRanges of blacklisted regions
#' @param genome A BSgenome object, or the full name of an installed BSgenome
#'   data package, or a short string specifying the name of an NCBI assembly
#'   (e.g. "GRCh38", "TAIR10.1", etc...) or UCSC genome (e.g. "hg38", "bosTau9",
#'   "galGal6", "ce11", etc...). The supplied short string must refer
#'   unambiguously to an installed BSgenome data package. See
#'   \link[BSgenome]{getBSgenome}.
#' @param studySignal The median signal (number of fragments) in your study. If
#'   not set, this will be calculated using the input ArchR project but relies
#'   on the assumption that the ArchR project encompasses your whole study (i.e.
#'   is not a subset).
#' @param generalizeStudySignal If `studySignal` is not provided, calculate the
#'   signal as the mean of the mean & median number of fragments for of
#'   individual samples within each cell population. This may improve MOCHA's
#'   ability to generalize to datasets with XXXXXX #TODO. Default is FALSE, use
#'   the median number of fragments.
#' @param cellCol The column in cellColData specifying unique cell ids or
#'   barcodes. Default is "RG", the unique cell identifier used by ArchR.
#' @param TxDb The exact package name of a TxDb-class transcript annotation
#'   package for your organism (e.g. "TxDb.Hsapiens.UCSC.hg38.refGene"). This
#'   must be installed. See
#'   \href{https://bioconductor.org/packages/release/data/annotation/}{
#'   Bioconductor AnnotationData Packages}.
#' @param OrgDb The exact package name of a OrgDb-class genome wide annotation
#'   package for your organism (e.g. "org.Hs.eg.db"). This must be installed.
#'   See \href{https://bioconductor.org/packages/release/data/annotation/}{
#'   Bioconductor AnnotationData Packages}
#' @param outDir is a string describing the output directory for coverage files.
#'   Must be a complete directory string. With ArchR input, set outDir to NULL
#'   to create a directory within the input ArchR project directory named MOCHA
#'   for saving files.
#' @param numCores integer. Number of cores to parallelize peak-calling across
#'   multiple cell populations.
#' @param force Optional, whether to force creation of coverage files if they
#'   already exist. Default is FALSE.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param peakModel Optional \code{MOCHAPeakModel} from \code{\link{trainPeakModel}}.
#'   When \code{NULL}, the bundled 500 bp model is used. Custom models set the
#'   tile width and study-signal calibration (\code{trainingMedian}).
#'
#' @return tileResults A MultiAssayExperiment object containing ranged data for
#'   each tile
#' @examples
#' \dontrun{
#' # Starting from an ArchR Project:
#' tileResults <- MOCHA::callOpenTiles(
#'   ArchRProj = myArchRProj,
#'   cellPopLabel = "celltype_labeling",
#'   cellPopulations = "CD4",
#'   TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene",
#'   OrgDb = "org.Hs.eg.db",
#'   numCores = 1
#' )
#' }
#' \donttest{
#' # Starting from GRangesList:
#' if (
#'   requireNamespace("BSgenome.Hsapiens.UCSC.hg19") &&
#'   requireNamespace("TxDb.Hsapiens.UCSC.hg38.refGene") &&
#'   requireNamespace("org.Hs.eg.db")
#' ) {
#'   tiles <- MOCHA::callOpenTiles(
#'     ATACFragments = MOCHA::exampleFragments,
#'     cellColData = MOCHA::exampleCellColData,
#'     blackList = MOCHA::exampleBlackList,
#'     genome = "BSgenome.Hsapiens.UCSC.hg19",
#'     TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene",
#'     OrgDb = "org.Hs.eg.db",
#'     outDir = tempdir(),
#'     cellPopLabel = "Clusters",
#'     cellPopulations = c("C2", "C5"),
#'     numCores = 1
#'   )
#' }
#' # Sample-level GRangesList (one element per sample):
#' # sampleFragments <- GRangesList(Sample1 = frags1, Sample2 = frags2)
#' # callOpenTiles(sampleFragments, cellColData = meta, cellPopLabel = "Clusters", ...)
#' }
#'
#' @export
#' @docType methods
#' @rdname callOpenTiles-methods
#' @keywords core
setGeneric(
  "callOpenTiles",
  function(ATACFragments,
           cellColData,
           blackList,
           genome,
           cellPopLabel,
           cellPopulations = "ALL",
           sampleColumn = 'Sample',
           studySignal = NULL,
           generalizeStudySignal = FALSE,
           cellCol = "RG",
           TxDb,
           OrgDb,
           outDir,
           numCores = 30,
           verbose = FALSE,
           force = FALSE,
           peakModel = NULL) {
    standardGeneric("callOpenTiles")
  },
  signature = "ATACFragments"
)

#' Parse CellPopulation and Sample from legacy ATACFragments element names.
#' @noRd
.parseLegacyFragmentNames <- function(fragNames) {
  cellPopList <- vapply(fragNames, function(x) {
    unlist(stringr::str_split(x, "#"))[1]
  }, character(1))
  sampleList <- vapply(fragNames, function(x) {
    parts <- unlist(stringr::str_split(x, "#"))
    if (length(parts) < 2) {
      return(NA_character_)
    }
    unlist(stringr::str_split(parts[2], "__"))[1]
  }, character(1))
  list(cellPopulations = cellPopList, samples = sampleList)
}

#' Detect whether ATACFragments use legacy or sample-level naming.
#' @noRd
.detectATACFragmentsInputMode <- function(
    ATACFragments,
    cellColData,
    cellPopLabel,
    sampleColumn) {
  fragNames <- names(ATACFragments)
  if (is.null(fragNames) || any(!nzchar(fragNames))) {
    stop("ATACFragments must be a named list or GRangesList.")
  }

  parsed <- .parseLegacyFragmentNames(fragNames)
  metaSamples <- unique(as.character(cellColData[[sampleColumn]]))
  metaCellPops <- unique(as.character(cellColData[[cellPopLabel]]))

  legacyOk <- all(grepl("#", fragNames, fixed = TRUE)) &&
    all(parsed$cellPopulations %in% metaCellPops, na.rm = TRUE) &&
    all(parsed$samples %in% metaSamples, na.rm = TRUE)

  sampleLevelOk <- all(fragNames %in% metaSamples) &&
    !any(grepl("#", fragNames, fixed = TRUE))

  if (legacyOk) {
    return("legacy")
  }
  if (sampleLevelOk) {
    return("sample_level")
  }

  if (all(fragNames %in% metaSamples)) {
    stop(
      "ATACFragments names match sample IDs in cellColData but also contain '#'. ",
      "Use sample-level names without '#', or legacy names in format ",
      "`CellPopulation#Sample`."
    )
  }

  unknownSamples <- setdiff(fragNames, metaSamples)
  if (length(unknownSamples) > 0 && !any(grepl("#", fragNames, fixed = TRUE))) {
    stop(
      "Sample names in ATACFragments not found in cellColData[[", sampleColumn, "]]: ",
      paste(unknownSamples, collapse = ", ")
    )
  }

  stop(
    "Could not interpret ATACFragments names. Provide either sample-level names ",
    "(matching cellColData[[sampleColumn]]) or legacy names in format ",
    "`CellPopulation#Sample`."
  )
}

#' Cell identifiers from cellColData (rownames or cellCol column).
#' @noRd
.getCellColDataCellIds <- function(cellColData, cellCol) {
  meta <- as.data.frame(cellColData)
  cellIds <- rownames(meta)
  if (is.null(cellIds) || !any(nzchar(cellIds))) {
    if (cellCol %in% colnames(meta)) {
      cellIds <- as.character(meta[[cellCol]])
    } else {
      stop(
        "cellColData must have rownames matching cell barcodes, or include ",
        "column '", cellCol, "'."
      )
    }
  }
  cellIds
}

#' Split sample-level fragments into legacy CellPopulation#Sample GRangesList.
#' @noRd
.normalizeSampleLevelFragments <- function(
    ATACFragments,
    cellColData,
    cellPopLabel,
    sampleColumn,
    cellCol,
    cellPopulations,
    verbose = FALSE) {
  meta <- as.data.frame(cellColData)
  cellIds <- .getCellColDataCellIds(cellColData, cellCol)
  rownames(meta) <- cellIds

  if (all(tolower(cellPopulations) == "all")) {
    cellPops <- unique(as.character(meta[[cellPopLabel]]))
    cellPops <- cellPops[!is.na(cellPops)]
  } else {
    cellPops <- cellPopulations
    missingPops <- cellPops[!cellPops %in% meta[[cellPopLabel]]]
    if (length(missingPops) > 0) {
      stop(
        "cellPopulations not found in cellColData column '", cellPopLabel, "': ",
        paste(missingPops, collapse = ", ")
      )
    }
  }

  metaSamples <- unique(as.character(meta[[sampleColumn]]))
  unknownSamples <- setdiff(names(ATACFragments), metaSamples)
  if (length(unknownSamples) > 0) {
    stop(
      "Sample names in ATACFragments not found in cellColData[[", sampleColumn, "]]: ",
      paste(unknownSamples, collapse = ", ")
    )
  }

  out <- list()
  for (sampleName in names(ATACFragments)) {
    frags <- ATACFragments[[sampleName]]
    if (length(frags) == 0) {
      next
    }

    fragCells <- unique(as.character(GenomicRanges::mcols(frags)[[cellCol]]))
    sampleCells <- cellIds[meta[[sampleColumn]] == sampleName]
    unknownCells <- setdiff(fragCells, sampleCells)
    if (length(unknownCells) > 0) {
      stop(
        "Fragments in sample '", sampleName, "' contain cell IDs not found in ",
        "cellColData: ",
        paste(head(unknownCells, 5), collapse = ", "),
        if (length(unknownCells) > 5) paste0(" (and ", length(unknownCells) - 5, " more)") else ""
      )
    }

    for (cp in cellPops) {
      cellsInPop <- cellIds[
        meta[[sampleColumn]] == sampleName & meta[[cellPopLabel]] == cp
      ]
      if (length(cellsInPop) == 0) {
        next
      }
      idx <- as.character(GenomicRanges::mcols(frags)[[cellCol]]) %in% cellsInPop
      subFrags <- frags[idx]
      if (length(subFrags) == 0) {
        next
      }
      key <- paste(cp, sampleName, sep = "#")
      out[[key]] <- subFrags
    }
  }

  if (length(out) == 0) {
    stop(
      "No fragments remained after splitting sample-level ATACFragments by ",
      "cell population. Check cellPopLabel, cellPopulations, and cellColData."
    )
  }

  if (verbose) {
    message(
      "Converted sample-level ATACFragments to ",
      length(out),
      " CellPopulation#Sample fragment sets."
    )
  }

  GenomicRanges::GRangesList(out)
}

#' @rdname callOpenTiles-methods
#' @aliases callOpenTiles, GRangesList-method
.callOpenTiles_default <- function(ATACFragments,
                                   cellColData,
                                   blackList,
                                   genome,
                                   cellPopLabel,
                                   cellPopulations = "ALL",
                                   sampleColumn = 'Sample',
                                   studySignal = NULL,
                                   generalizeStudySignal = FALSE,
                                   cellCol = "RG",
                                   TxDb,
                                   OrgDb,
                                   outDir,
                                   numCores = 30,
                                   verbose = FALSE,
                                   force = FALSE,
                                   peakModel = NULL) {
  Sample <- seqnames <- NULL

  genome <- BSgenome::getBSgenome(genome)

  validGRanges <- sapply(ATACFragments, function(obj) {
    all(
      methods::is(obj, "GRanges"),
      cellCol %in% colnames(GenomicRanges::mcols(obj))
    )
  })
  if (!all(validGRanges)) {
    stop(
      "Invalid ATACFragments. ATACFragments must be a list of GRanges or ",
      "GRangesList, each containing a column given by 'cellCol' containing ",
      "unique cell IDs."
    )
  }

  if (!methods::is(blackList, "GRanges")) {
    stop(
      "Invalid blackList. blackList must be a GRanges"
    )
  }
  if (!(sampleColumn %in% colnames(cellColData))) {
    stop("Invalid cellColData. sampleColumn must be found within cellColData.")
  }
  if (!(cellPopLabel %in% colnames(cellColData))) {
    stop(
      stringr::str_interp("cellPopLabel {cellPopLabel} not found in "),
      "cellColData. cellColData must contain column 'cellPopLabel'."
    )
  }

  if (!file.exists(outDir)) {
    if (verbose) {
      message(stringr::str_interp("Creating directory for MOCHA at ${outDir}"))
    }
    dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
  }

  # Filter out fragments that are not aligned to the Genome.
  allnames <- names(ATACFragments)
  beforeLengths <- lengths(ATACFragments)
  ATACFragments <- lapply(ATACFragments, function(x) {
    dplyr::filter(x, seqnames %in% GenomeInfoDb::seqnames(genome))
  })
  names(ATACFragments) <- allnames

  if (verbose & sum(beforeLengths) - sum(lengths(ATACFragments)) > 0) {
    warning(
      "Some fragments are not aligned to BS.Genome provided. ",
      sum(beforeLengths) - sum(lengths(ATACFragments)),
      " fragments were removed."
    )
  }

  if (verbose & any(lengths(ATACFragments) == 0)) {
    warning(
      "Some cell type sample combinations did not have any aligned fragments, including ",
      paste(names(ATACFragments)[lengths(ATACFragments) == 0], collapse = ", ")
    )
  }

  inputMode <- .detectATACFragmentsInputMode(
    ATACFragments,
    cellColData,
    cellPopLabel,
    sampleColumn
  )
  if (inputMode == "sample_level") {
    ATACFragments <- .normalizeSampleLevelFragments(
      ATACFragments,
      cellColData,
      cellPopLabel,
      sampleColumn,
      cellCol,
      cellPopulations,
      verbose = verbose
    )
  }

  parsedNames <- .parseLegacyFragmentNames(names(ATACFragments))
  cellPopList <- parsedNames$cellPopulations
  sampleList <- parsedNames$samples

  # Verify that cell populations and samples in ATACFragments names match those
  # in cellColData
  if (!all(cellPopList %in% cellColData[[cellPopLabel]])) {
    stop(
      "Cell populations in names of ATACFragments do not match those in cellColData. ",
      "Names must be in format `CellPopulation#Sample`, or use sample-level ",
      "names matching cellColData[[", sampleColumn, "]]."
    )
  }
  if (!all(sampleList %in% cellColData[[sampleColumn]])) {
    stop(
      "Sample names in names of ATACFragments do not match those in cellColData. ",
      "Names must be in format `CellPopulation#Sample`, or use sample-level ",
      "names matching cellColData[[", sampleColumn, "]]."
    )
  }

  .callOpenTiles(
    ATACFragments,
    cellColData,
    blackList,
    genome,
    cellPopLabel,
    cellPopulations,
    sampleColumn,
    studySignal,
    generalizeStudySignal,
    cellCol,
    TxDb,
    OrgDb,
    outDir,
    numCores,
    verbose,
    force,
    useArchR = FALSE,
    peakModel = peakModel
  )
}
#' @rdname callOpenTiles-methods
#' @aliases callOpenTiles,GRangesList-method
setMethod(
  "callOpenTiles",
  signature(ATACFragments = "GRangesList"),
  .callOpenTiles_default
)
#' @rdname callOpenTiles-methods
#' @aliases callOpenTiles,list-method
setMethod(
  "callOpenTiles",
  signature(ATACFragments = "list"),
  .callOpenTiles_default
)


#' @rdname callOpenTiles-methods
#' @aliases callOpenTiles,ArchRProject-method
.callOpenTiles_ArchR <- function(ATACFragments,
                                 cellPopLabel,
                                 cellPopulations = "ALL",
                                 sampleColumn = 'Sample',
                                 studySignal = NULL,
                                 generalizeStudySignal = FALSE,
                                 TxDb,
                                 OrgDb,
                                 outDir = NULL,
                                 numCores = 30,
                                 verbose = FALSE,
                                 force = FALSE,
                                 peakModel = NULL) {
  Sample <- nFrags <- NULL
  # Load Genome
  genome <- ArchR::validBSgenome(ArchR::getGenome(ATACFragments))

  if (is.null(outDir)) {
    outDir <- paste(
      ArchR::getOutputDirectory(ATACFragments),
      "/MOCHA",
      sep = ""
    )
  }

  if (!file.exists(outDir)) {
    if (verbose) {
      message(stringr::str_interp("Creating directory for MOCHA at ${outDir}"))
    }
    dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
  }

  # Get cell metadata and blacklisted regions from ArchR Project
  cellColData <- ArchR::getCellColData(ATACFragments)
  blackList <- ArchR::getBlacklist(ATACFragments)

  if (!(cellPopLabel %in% colnames(cellColData))) {
    stop(
      stringr::str_interp("cellPopLabel ${cellPopLabel} not found in "),
      "cellColData. cellColData must contain column cellPopLabel."
    )
  }

  .callOpenTiles(
    ATACFragments,
    cellColData,
    blackList,
    genome,
    cellPopLabel,
    cellPopulations,
    sampleColumn,
    studySignal,
    generalizeStudySignal,
    cellCol = "RG",
    TxDb,
    OrgDb,
    outDir,
    numCores,
    verbose,
    force,
    useArchR = TRUE,
    peakModel = peakModel
  )
}
setMethod(
  "callOpenTiles",
  signature(ATACFragments = "ArchRProject"),
  .callOpenTiles_ArchR
)


.callOpenTiles <- function(ATACFragments,
                           cellColData,
                           blackList,
                           genome,
                           cellPopLabel,
                           cellPopulations,
                           sampleColumn,
                           studySignal,
                           generalizeStudySignal,
                           cellCol,
                           TxDb,
                           OrgDb,
                           outDir,
                           numCores,
                           verbose,
                           force,
                           useArchR,
                           peakModel) {
  Sample <- meanValues <- NULL
  peakModel <- if (is.null(peakModel)) {
    .defaultPeakModel()
  } else {
    .validatePeakModel(peakModel)
  }
  trainingMedian <- peakModel$trainingMedian
  # Load databases and save names for use in metadata
  TxDbName <- TxDb
  OrgDbName <- OrgDb
  TxDb <- getAnnotationDbFromInstalledPkgname(dbName = TxDb, type = "TxDb")
  OrgDb <- getAnnotationDbFromInstalledPkgname(dbName = OrgDb, type = "OrgDb")
  
  if (any(is.na(cellColData[[cellPopLabel]]))){
    
    warning(
      stringr::str_interp("Some cells within the column ${cellPopLabel} are labeled as NA. Those cells will be ignored.")
    )
    
  }

  # Get cell populations
  cellTypeLabelList <- cellColData[, cellPopLabel]
  
  if (!all(tolower(cellPopulations) == "all")) {
    if (!all(cellPopulations %in%  unique(cellTypeLabelList))) {
      missingCellPopulations <- cellPopulations[
        !cellPopulations %in%  unique(cellTypeLabelList)
      ]
      stop(
        stringr::str_interp(paste0(
          "Some or all of the cell populations provided were not found in the ",
          "cellColData column '${cellPopLabel}'. Missing cell populations: "
        )),
        paste0(missingCellPopulations, collapse=", "), "."
      )
    }
  }

  #################
  # Begin constructing the additional metadata SummarizedExperiment
  # Done prior to the expensive computation so there's no heartbreak if this
  # step fails.

  cellCounts <- as.data.frame(table(cellColData[, sampleColumn], cellTypeLabelList))
  names(cellCounts) <- c(sampleColumn, "CellPop", "CellCount")
  cellCounts <- tidyr::pivot_wider(
    cellCounts,
    id_cols = "CellPop",
    names_from = sampleColumn,
    values_from = "CellCount"
  )
  allCellCounts <- as.data.frame(cellCounts[, -1])
  rownames(allCellCounts) <- cellCounts$CellPop

  # Create a dummy data.frame for storing fragment information
  # as it is calculated later (when iterating over cell populations)
  allFragmentCounts <- allCellCounts
  allFragmentCounts[!is.na(allFragmentCounts)] <- NA

  # Filter allCellCounts and allFragmentCounts to just cell populations (rows)
  # of interest
  if (all(tolower(cellPopulations) == "all")) {
    cellPopulations <- rownames(allCellCounts)
  } else if(any(rownames(allCellCounts) %in% cellPopulations)){
    allCellCounts <- allCellCounts[rownames(allCellCounts) %in% cellPopulations, , drop = FALSE]
    allFragmentCounts <- allFragmentCounts[rownames(allFragmentCounts) %in% cellPopulations, , drop = FALSE]
  } else{
  
    stop(
      stringr::str_interp("Some or all of the cell populations provided were not found in the cellColData column ${cellPopLabel}. 
        These cell populations include ${cellPopulations}.")
    )
      
  }

  # For additional metadata:
  # Some numeric columns may be stored as character - convert these to numeric
  # Make a copy to preserve original columns.

  # This filtering introduced a bug where a cell population was only present
  # in one sample - resulted in entire samples being dropped from
  # summarizedData.
  # cellColDataCopy <- as.data.frame(dplyr::filter(data.frame(cellColData),
  # !!as.name(cellPopLabel) %in% cellPopulations))
  cellColDataCopy <- copy(cellColData)

  cellColDataCopy[] <- lapply(cellColDataCopy, function(x) {
    utils::type.convert(as.character(x), as.is = TRUE)
  })

  # Assume all numeric columns are to be saved as additionalCellData, unless it's unique to a sample. 
  isNumericCol <- unlist(lapply(cellColDataCopy, function(x) is.numeric(x)))
  additionalCellData <- colnames(cellColDataCopy)[isNumericCol]
                                
  # Group by Sample (rows) and cellPop (columns)
  if (!is.null(additionalCellData)) {
    if (verbose) {
      message(
        "Summarizing additional metadata: ",
        paste(additionalCellData, collapse = ", ")
      )
    }
    additionalMetaData <- lapply(additionalCellData, function(x) {
      suppressMessages(
        summarizedData <- as.data.frame(cellColDataCopy) %>%
          dplyr::group_by(!!as.name(cellPopLabel), !!as.name(sampleColumn)) %>%
          dplyr::summarize(meanValues = mean(!!as.name(x), .groups = "drop")) %>%
          tidyr::pivot_wider(
            id_cols = tidyselect::all_of(cellPopLabel),
            names_from = !!as.name(sampleColumn),
            values_from = meanValues
          )
      )

      summarizedData <- as.data.frame(summarizedData)
        
      #Remove any columns with missing labels
      cellTypeList <- summarizedData[[cellPopLabel]]
      summarizedData <- summarizedData[!is.na(cellTypeList),,drop=FALSE]
      rownames(summarizedData) <- summarizedData[[cellPopLabel]]
      summarizedData <- summarizedData[, -1, drop = FALSE]

      # Filter to specific cellPopulations
      summarizedData <- summarizedData[
        rownames(summarizedData) %in% cellPopulations, ,
        drop = FALSE
      ]

      summarizedData
    })
    names(additionalMetaData) <- additionalCellData
  } else if (is.null(additionalCellData)) {
    additionalMetaData <- NULL
  }
  remove(cellColDataCopy)
  ################# End metadata construction

  # Add prefactor multiplier across datasets
  if (is.null(studySignal)) {
    if (generalizeStudySignal) {
      if (verbose) {
        message(
          "Parameter `studySignal` was not provided. ",
          "Calculating study signal on cellColData as the mean of the mean ",
          "and median nFrags of individual samples within each cell population."
        )
      }
      study_prefactor <- NULL
    } else {
      # Calculate study prefactor with study-wide median nfrags
      if (verbose) {
        message(
          "Parameter `studySignal` was not provided. ",
          "Calculating study signal on cellColData as the median ",
          "nFrags with the assumption that all cell populations are ",
          "present in cellColData."
        )
      }
      if (!("nFrags" %in% colnames(cellColData))) {
        stop(
          "cellColData is missing fragment count information. ",
          "To calculate study signal, cellColData must contain a column",
          " 'nFrags' representing the number of fragments per cell. ",
          "Alternatively, provide a value for the parameter 'studySignal'."
        )
      }
      studySignal <- stats::median(cellColData$nFrags)
      study_prefactor <- trainingMedian / studySignal
    }
  } else {
    if (generalizeStudySignal) {
      if (verbose) {
        message(
          "Ignoring provided studySignal since `generalizeStudySignal` = TRUE. ",
          "Calculating study signal on cellColData as the mean of the mean ",
          "and median nFrags of individual samples within each cell population."
        )
      }
      study_prefactor <- NULL
    } else {
      # Use user-provided studySignal
      study_prefactor <- trainingMedian / studySignal
    }
  }


  # Main loop over all cell populations
  experimentList <- list()
  for (cellPop in cellPopulations) {
    if (verbose) {
      message(stringr::str_interp(
        "Calling open tiles for cell population ${cellPop}"
      ))
    }

    cl <- parallel::makeCluster(numCores)
    tryCatch({
      if (useArchR) {
        parallel::clusterEvalQ(cl, {
          library(ArchR)
          library(rhdf5)
        })
        # Get our fragments for this cellPop
        frags <- MOCHA::getPopFrags(
          ArchRProj = ATACFragments,
          cellPopLabel = cellPopLabel,
          cellSubsets = cellPop,
          numCores = cl,
          returnGRangesList = FALSE,
          verbose = verbose
        )
      } else {
        cellPopList <- sapply(names(ATACFragments), function(x) {
          unlist(stringr::str_split(x, "#"))[1]
        })
        frags <- ATACFragments[which(cellPop == cellPopList)]
        if (length(frags) == 0) {
          stop(
            "Provided ATACFragments does not contain any sample fragments ",
            stringr::str_interp("belonging to cellPopulations: ${cellPop}. "),
            "ATACFragments must have at least one sample GRanges for each ",
            "cell population."
          )
        }
      }

      # Simplify sample names to remove celltype and normalization factor
      # Removes everything before the first "#" and after the first "__"
      # This is a convention from MOCHA::getPopFrags which is ArchR-dependent
      # E.g. C1#PBMCSmall__0.527239 becomes PBMCSmall
      sampleNames <- gsub("__.*", "", gsub(".*#", "", names(frags)))
      names(frags) <- sampleNames

      ##Check whether fragments lists are hashed or not.
      ##If fragments are hashed, and the user is using ArchR, then dehash them by true sample here.
      ## If ArchR is not being used, the user is expected to do this themselves.
      if (all(!sampleNames %in% unique(cellColData[, sampleColumn])) & useArchR) {
        frags <- dehashArchR(frags, cellColData, sampleColumn, cl = cl)

        ## Check if any samples are now empty of fragments for a given cell type.
        ## Remove that index from the list, and then update the sample names for this cell type.
        if (any(lengths(frags) == 0)) {
          frags <- frags[lengths(frags) != 0]
        }

        sampleNames <- names(frags)
      }

      # Calculate normalization factors as the number of fragments for
      # each celltype_sample
      normalization_factors <- as.integer(lengths(frags))

      # Assign the number of fragments into the fragment count
      allFragmentCounts[cellPop, sampleNames] <- normalization_factors

      # save coverage files to folder.
      if (!file.exists(
        paste(outDir, "/", cellPop, "_CoverageFiles.RDS", sep = "")
      ) || force) {
        if (verbose) {
          message(stringr::str_interp(
            "Saving coverage files for cell population ${cellPop}"
          ))
        }
        covFiles <- getCoverage(
          popFrags = frags,
          normFactor = normalization_factors / 10^6,
          filterEmpty = TRUE,
          cl = cl, TxDb = TxDb
        )
        saveRDS(
          covFiles,
          paste(outDir, "/", cellPop, "_CoverageFiles.RDS", sep = "")
        )
        rm(covFiles)
      }

      if (is.null(study_prefactor)) {
        if (verbose) {
          message("Calculating generalized study signal...")
        }
        allmeans <- list()
        allmedians <- list()
        for (sample in names(frags)) {
          fragsdf <- as.data.frame(frags[[sample]])
          cellFragsTable <- table(fragsdf[[cellCol]])
          allmeans <- append(allmeans, mean(cellFragsTable))
          allmedians <- append(allmedians, stats::median(cellFragsTable))
        }
        mean_nfrags <- mean(unlist(allmeans))
        median_nfrags <- mean(unlist(allmedians))
        combinedSignal <- mean(c(median_nfrags, mean_nfrags))
        study_prefactor <- trainingMedian / combinedSignal

        if (verbose) {
          message(
            "Mean fragments per cell: ", mean_nfrags,
            "\nMedian fragments per cell: ", median_nfrags,
            "\nCombined signal: ", combinedSignal
          )
        }
        rm(fragsdf)
        rm(cellFragsTable)
      }

      iterList <- lapply(seq_along(frags), function(x) {
        list(blackList, frags[[x]], cellCol, verbose, study_prefactor, peakModel)
      })

      tilesGRangesList <- pbapply::pblapply(
        cl = cl,
        X = iterList,
        FUN = simplifiedTilesBySample
      )
    }, finally = {
      try(parallel::stopCluster(cl), silent = TRUE)
    })
    gc()
    names(tilesGRangesList) <- names(frags)

    # Where samples have no cells, add an empty GRanges placeholder
    if (!all(colnames(allCellCounts) %in% names(tilesGRangesList))) {
      emptySamples <- colnames(allCellCounts)[
        !colnames(allCellCounts) %in% names(tilesGRangesList)
      ]
      emptyGRanges <- lapply(emptySamples, function(x) {
        NULL
      })
      names(emptyGRanges) <- emptySamples
      tilesGRangesList <- append(tilesGRangesList, emptyGRanges)
    }

    tilesGRangesList <- tilesGRangesList[base::sort(names(tilesGRangesList))]

    # Cannot make peak calls with < 5 cells (see make_prediction.R)
    # so NULL will occur for those samples. We need to fill in
    # dummy data so that we preserve the existence of the sample, while
    # also not including any information from it.
    emptyGroups <- which(unlist(lapply(tilesGRangesList, is.null)))

    if (length(emptyGroups) > 0) {
      if (verbose) {
        warning(
          "The following samples have too few cells (<5) of this celltype (",
          cellPop,
          ") and will be ignored: ",
          paste(names(tilesGRangesList)[emptyGroups], collapse = ", ")
        )
      }
    }

    for (i in emptyGroups) {
      # This is an empty region placeholder that represents an empty sample
      # And is functionally ignored downstream but required for
      # the RaggedExperiment structure
      tilesGRangesList[[i]] <- data.table(
        tileID = "chr1:1-499", seqnames = "chr1", start = 1,
        end = 499, strand = "*", TotalIntensity = 0, maxIntensity = 0,
        numCells = 0, Prediction = 0, PredictionStrength = 0, peak = FALSE
      )
    }

    # Package rangeList into a RaggedExperiment
    ragExp <- RaggedExperiment::RaggedExperiment(
      tilesGRangesList
    )

    # And add it to the experimentList for this cell population
    experimentList <- append(experimentList, ragExp)
  }

  allFragmentCounts[is.na(allFragmentCounts)] <- 0
                                
  # Create sample metadata from cellColData using util function
  # "Sample" is the enforced col in ArchR containing the
  # sample IDs which correspond to the arrow files.
  # We are assuming samples are synonymous to wells
  # (If hashed, samples were un-hashed during ArchR project generation)
  sampleData <- suppressWarnings(
    sampleDataFromCellColData(cellColData, sampleLabel = sampleColumn)
  )

 # Match cell populations in allCellCounts/allFragmentCounts to those
  # in additionalMetaData, in case of NA cell populations
  if (length(additionalMetaData) >= 1) {
    allCellCounts <- allCellCounts[
      match(rownames(additionalMetaData[[1]]), rownames(allCellCounts)),
      ,
      drop = FALSE
    ]
    
    allFragmentCounts <- allFragmentCounts[
      match(rownames(additionalMetaData[[1]]), rownames(allFragmentCounts)),
      ,
      drop = FALSE
    ]
  }

 ### Let's put a full list together, and then order all of the additional metadata the same way.
  sumDataAssayList <- append(
    list(
      "CellCounts" = allCellCounts,
      "FragmentCounts" = allFragmentCounts
    ),
    additionalMetaData
  )

  # Enforce Row and Column orders in sumDataAssayList and sampleData
  colOrder <- colnames(allCellCounts)
  rowOrder <- rownames(allCellCounts)

  for (i in seq_along(sumDataAssayList)) {
    assayName <- names(sumDataAssayList[i])
    assay <- sumDataAssayList[[i]]
    sumDataAssayList[assayName] <- list(assay[rowOrder, colOrder, drop = FALSE])
  }
                                
  if(sampleColumn != 'Sample'){
  
    sampleData$Sample = sampleData[,sampleColumn]
      
  }

  sampleData <- dplyr::arrange(
    sampleData, factor(Sample, levels = colOrder)
  )
  rownames(sampleData) <- sampleData[, sampleColumn]

  # Validate Row and Column orders
  if (verbose) {
    if (!all(rownames(sampleData) == colnames(allCellCounts))) {
      warning(
        "SampleData and allCellCounts samples mismatch:",
        "sampleData rownames:", rownames(sampleData),
        "allCellCounts colnames:", rownames(allCellCounts)
      )
    }
    if (length(unique(lapply(sumDataAssayList, dim))) > 1) {
      warning(
        "Assays in sumDataAssayList have different dimensions: ",
        paste(unique(lapply(sumDataAssayList, dim)), collapse = "\n")
      )
    }
    if (length(unique(lapply(sumDataAssayList, rownames))) > 1) {
      warning(
        "Assays in sumDataAssayList have different rownames: ",
        paste(unique(lapply(sumDataAssayList, rownames)), collapse = "\n")
      )
    }
    if (length(unique(lapply(sumDataAssayList, colnames))) > 1) {
      warning(
        "Assays in sumDataAssayList have different colnames: ",
        paste(unique(lapply(sumDataAssayList, colnames)), collapse = "\n")
      )
    }
  }
  
  summarizedData <- SummarizedExperiment::SummarizedExperiment(
   sumDataAssayList,
    colData = sampleData
  )

  # Add experimentList to MultiAssayExperiment
  names(experimentList) <- cellPopulations
  tileResults <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = experimentList,
    colData = sampleData,
    metadata = list(
      "summarizedData" = summarizedData,
      "Genome" = genome@pkgname,
      "TxDb" = list(pkgname = TxDbName, metadata = S4Vectors::metadata(TxDb)),
      "OrgDb" = list(pkgname = OrgDbName, metadata = S4Vectors::metadata(OrgDb)),
      "Directory" = outDir,
      "History" = list(paste("callOpenTiles", utils::packageVersion("MOCHA")))
    )
  )
  return(tileResults)
}
