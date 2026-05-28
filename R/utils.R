# Internal helpers for SampleTileObject coverage workflows

# Map gene identifiers using either a Bioconductor OrgDb (via AnnotationDbi)
# or a biomaRt::Mart connection. Routes by class so callers can transparently
# pass either, supporting offline (AnnotationDbi) and online (biomaRt) workflows.
.map_gene_ids <- function(db, keys, column = "SYMBOL", keytype = "ENTREZID") {
  if (methods::is(db, "Mart")) {
    if (!requireNamespace("biomaRt", quietly = TRUE)) {
      stop("Package 'biomaRt' is required when passing a Mart connection.")
    }
    attr_col <- .biomart_attribute(column)
    attr_key <- .biomart_attribute(keytype)
    res <- biomaRt::getBM(
      attributes = unique(c(attr_key, attr_col)),
      filters = attr_key,
      values = as.character(keys),
      mart = db
    )
    out <- res[[attr_col]][match(as.character(keys), as.character(res[[attr_key]]))]
    out[!nzchar(out)] <- NA_character_
    return(out)
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop(
      "Gene-id mapping requires either AnnotationDbi (with an OrgDb) ",
      "or a biomaRt::useEnsembl() Mart passed in place of the OrgDb. ",
      "Install AnnotationDbi or provide a Mart."
    )
  }
  suppressWarnings(AnnotationDbi::mapIds(db, keys, column, keytype))
}

# biomaRt attribute names for the OrgDb-style labels used across MOCHA.
.biomart_attribute <- function(label) {
  switch(
    toupper(label),
    "ENTREZID" = "entrezgene_id",
    "ENSEMBL" = "ensembl_gene_id",
    "ENSEMBLTRANS" = "ensembl_transcript_id",
    "TXNAME" = "ensembl_transcript_id",
    "SYMBOL" = "external_gene_name",
    "GENENAME" = "external_gene_name",
    "REFSEQ" = "refseq_mrna",
    stop(sprintf(
      "No biomaRt attribute is configured for OrgDb label '%s'. Pass an OrgDb (AnnotationDbi) or extend .biomart_attribute().",
      label
    ))
  )
}

.resolveSubSamples <- function(metaFile, groupColumn = NULL, subGroups = NULL) {
  if (!is.null(subGroups) & !is.null(groupColumn)) {
    subSamples <- lapply(subGroups, function(x) {
      metaFile[metaFile[, groupColumn] %in% x, "Sample"]
    })
    names(subSamples) <- subGroups
    subGroupsOut <- subGroups
  } else if (!is.null(groupColumn)) {
    subGroupsOut <- unique(metaFile[, groupColumn])
    subSamples <- lapply(subGroupsOut, function(x) {
      metaFile[metaFile[, groupColumn] %in% x, "Sample"]
    })
    names(subSamples) <- subGroupsOut
  } else {
    subGroupsOut <- "All"
    subSamples <- list("All" = metaFile[, "Sample"])
  }

  list(subGroups = subGroupsOut, subSamples = subSamples)
}

.validateCoverageDirectory <- function(object, objectName = "SampleTileObj") {
  outDir <- object@metadata$Directory
  if (is.na(outDir)) {
    stop(
      "Missing coverage file directory. ",
      objectName,
      "$metadata must contain 'Directory'."
    )
  }
  if (!file.exists(outDir)) {
    stop(
      "Directory given by ",
      objectName,
      "@metadata$Directory does not exist."
    )
  }
  outDir
}

.requireCellPopulationAssays <- function(cellNames) {
  if (all(toupper(cellNames) == "COUNTS")) {
    stop(
      "The only assay in the SummarizedExperiment is Counts. The names of assays must reflect cell types,",
      " such as those in the Summarized Experiment output of getSampleTileMatrix."
    )
  }
}

.selectCoverageFromBundle <- function(bundle, coverage = TRUE) {
  if (coverage) {
    if ("Accessibility" %in% names(bundle)) {
      bundle[["Accessibility"]]
    } else {
      bundle
    }
  } else if ("Insertions" %in% names(bundle)) {
    bundle[["Insertions"]]
  } else {
    stop("Error around reading coverage files. Check that coverage files are not corrupted.")
  }
}

.selectAccessibilityCoverage <- function(bundle) {
  if ("Accessibility" %in% names(bundle)) {
    bundle$Accessibility
  } else {
    bundle
  }
}

.prepSampleTileGrouping <- function(metaFile, groupColumn = NULL, subGroups = NULL) {
  .resolveSubSamples(metaFile, groupColumn, subGroups)
}

.loadInsertionBias <- function(SampleTileObj, normTn5 = TRUE) {
  if (!normTn5) {
    return(list(insertBias = NULL, genome = NULL, genome_db = NULL))
  }

  if (any(grepl("InsertionBias", names(SampleTileObj@metadata)))) {
    genome_db <- SampleTileObj@metadata$Genome
    genome <- getAnnotationDbFromInstalledPkgname(dbName = genome_db, type = "BSgenome")
    insertBias <- SampleTileObj@metadata$InsertionBias
    insertBias <- insertBias[!is.na(insertBias[, "Norm"]), ]
    list(insertBias = insertBias, genome = genome, genome_db = genome_db)
  } else {
    stop("Attempting to normalize by Tn5 insertion bias, but no bias calculated. Please run addInsertionBias.")
  }
}

.callWithPackedArgs <- function(fn, packed) {
  do.call(fn, packed)
}

# Co-accessibility internal helpers
.normalizeTilePairs <- function(fullObj, tile1, tile2) {
  if (length(tile1) != length(tile2)) {
    stop("tile1 and tile2 must be the same length.")
  }

  if (is.character(tile1) && is.character(tile2)) {
    nTile1 <- match(tile1, rownames(fullObj))
    nTile2 <- match(tile2, rownames(fullObj))
  } else if (is.numeric(tile1) && is.numeric(tile2)) {
    nTile1 <- tile1
    nTile2 <- tile2
    tile1 <- rownames(fullObj)[nTile1]
    tile2 <- rownames(fullObj)[nTile2]
  } else {
    stop("tile1 and tile 2 must both be either numbers (indices) or strings")
  }

  list(tile1 = tile1, tile2 = tile2, nTile1 = nTile1, nTile2 = nTile2)
}

.validateBackNumber <- function(backNumber, fullObj, tile1, tile2, verbose = TRUE) {
  if (!is.null(dim(backNumber))) {
    return(backNumber)
  }

  if (backNumber >= length(rownames(fullObj)) - length(unique(c(tile1, tile2)))) {
    backNumber <- length(rownames(fullObj)) - length(unique(c(tile1, tile2)))
    if (verbose) {
      warning("backNumber too high. Reset to all background combinations.")
    }
  } else if (backNumber <= 10) {
    stop("backNumber too low (<=10). We recommend 1000.")
  }

  backNumber
}

.generateRandomBackgroundPairs <- function(accMat, tile1, tile2, backNumber, verbose = TRUE) {
  if (verbose) {
    message("Finding background peak pairs")
  }

  backGroundTiles <- rownames(accMat)[!rownames(accMat) %in% c(tile1, tile2)]
  backgroundCombos <- data.frame(
    Tile1 = sample(backGroundTiles, backNumber),
    Tile2 = sample(backGroundTiles, backNumber)
  )
  backgroundCombos[backgroundCombos[, 1] != backgroundCombos[, 2], , drop = FALSE]
}

.parseUserBackgroundPairs <- function(backNumber, fullObj, verbose = TRUE) {
  if (verbose) {
    message("Using user-defined background pairs")
  }

  backNumber <- as.data.frame(backNumber)
  if (!all(c("Tile1", "Tile2") %in% colnames(backNumber))) {
    stop("User-defined background pairs requires a column for Tile1 and Tile2")
  } else if (!all(grepl(":", c(backNumber[, "Tile1"], backNumber[, "Tile2"])) &
    grepl("-", backNumber[, "Tile1"], backNumber[, "Tile2"]))) {
    stop("User-defined background pairs must be in the form ChrX:100-2000")
  } else if (!all(c(backNumber[, "Tile1"], backNumber[, "Tile2"]) %in% rownames(fullObj))) {
    stop("User-defined background pairs includes regions not found within the sample tile accessibility matrix.")
  }

  backgroundCombos <- as.data.frame(backNumber)[, c("Tile1", "Tile2")]

  if (sum(backgroundCombos[, "Tile1"] != backgroundCombos[, "Tile2"]) < 10) {
    stop("User-defined background pairs are fewer than 10. Please provide a larger background.")
  }

  backgroundCombos[backgroundCombos[, "Tile1"] != backgroundCombos[, "Tile2"], , drop = FALSE]
}

.resolveBackgroundPairs <- function(backNumber, fullObj, tile1, tile2, verbose = TRUE) {
  if (is.null(dim(backNumber))) {
    .generateRandomBackgroundPairs(
      accMat = SummarizedExperiment::assays(fullObj)[[1]],
      tile1 = tile1,
      tile2 = tile2,
      backNumber = backNumber,
      verbose = verbose
    )
  } else if (dim(backNumber)[2] > 1) {
    .parseUserBackgroundPairs(backNumber, fullObj, verbose = verbose)
  } else {
    stop(paste(
      "Incorrect backNumber provided. Please provider either a number, or a data.frame",
      "with columns entitled Tile1 and Tile2, describing pairs to test.",
      "The tile names should be in the format ChrX:100-2000."
    ))
  }
}

.computeCoAccessibilityPValues <- function(foreGround, backGround, verbose = TRUE) {
  if (verbose) {
    message("Generating p-values.")
  }

  greatList <- unlist(pbapply::pblapply(
    foreGround$Correlation[which(foreGround$Correlation > 0)],
    function(x) {
      sum(x > backGround$Correlation)
    },
    cl = 1
  )) / length(backGround$Correlation)

  lesserList <- unlist(pbapply::pblapply(
    foreGround$Correlation[which(foreGround$Correlation < 0)],
    function(x) {
      sum(x < backGround$Correlation)
    },
    cl = 1
  )) / length(backGround$Correlation)

  foreGround$pValues <- rep(NA, length(foreGround$Correlation))
  foreGround$pValues[which(foreGround$Correlation > 0)] <- 1 - greatList
  foreGround$pValues[which(foreGround$Correlation < 0)] <- 1 - lesserList
  foreGround
}


# Function to get sample-level metadata,
# from an ArchR project's colData
sampleDataFromCellColData <- function(cellColData, sampleLabel) {
  if (!(sampleLabel %in% colnames(cellColData))) {
    stop(paste(
      "`sampleLabel` must be present in cellColData",
      "Check colnames(cellColData) for possible sample columns."
    ))
  }

  # Drop columns where all values are NA
  cellColDataNoNA <- BiocGenerics::Filter(function(x) {
    !all(is.na(x))
  }, cellColData)

  # Convert to data.table
  cellColDT <- data.table::as.data.table(cellColDataNoNA)
  BoolDT <- cellColDT[, lapply(.SD, function(x) {
      unique(x)
    length(unique(x)) == 1
  }), by = c(sampleLabel)]

  trueCols <- apply(BoolDT, 2, all)
  trueCols[[sampleLabel]] <- TRUE
  cellColDF <- as.data.frame(cellColDT)

  sampleData <- dplyr::distinct(cellColDF[, names(which(trueCols)), drop = F])

  # Set sampleIDs as rownames
  rownames(sampleData) <- sampleData[[sampleLabel]]
  return(sampleData)
}


# Function to de-hash a list of fragments files by the true cell type. 
# The function will identify the true list of samples, 
# And then iterate over each arrow file and select all the fragments from cells from that sample.
# And then combine the fragments into a new GRanges list. 
dehashArchR <- function(frags, cellColData, sampleColumn, cl = cl) {

  ##Identify cells for each true sample
  trueSamples <- unique(cellColData[,sampleColumn])
    
  trueCells <- lapply(trueSamples, function(XX){
                      trueCells = rownames(cellColData[cellColData[,sampleColumn] == XX,])
                  })
    
  trueFrags <- pbapply::pblapply(cl = cl, X = trueCells, dehashIter, oldfrags = frags)
  rm(frags)
    
  names(trueFrags) = trueSamples
  return(trueFrags)
}


# To effectively iterate over fragments without without memory leaks,
# This function is necessary. 
dehashIter <- function(cellIDs, oldfrags){
    RG <- NULL
    newFrags <- lapply(oldfrags, function(ZZ){
          
                  dplyr::filter(ZZ, RG %in% cellIDs)
          
          })
    
    newFrags <- unlist(methods::as(newFrags, 'GRangesList'))

    return(newFrags)
}


# To fix wierd cell type names within metadata, so that it doesn't break rownames and assay names
fixCellTypeNames <- function(cellNames){

    ##fix spaces and dashes
    cellNames = gsub(" |-","_", cellNames)
    
    #fix + and - signs
    cellNames = gsub("\\+","pos", cellNames)
    cellNames = gsub("\\-","neg", cellNames)
    cellNames = gsub("\\\\", "_", cellNames)

    return(cellNames)
}



# Function to split the output of getPopFrags into a list
# of lists of GRanges, one named for each celltype.
# Resulting in a list containing each celltype, and each
# celltype has a list of GRanges name for each sample.
splitFragsByCellPop <- function(frags) {

  # Rename frags by cell population
  renamedFrags <- lapply(
    1:length(frags),
    function(y) {
      # Split out celltype and sample from the name
      x <- frags[y]
      celltype_sample <- names(x)
      splits <- unlist(stringr::str_split(celltype_sample, "#"))
      celltype <- splits[1]
      sample <- unlist(stringr::str_split(splits[2], "__"))[1]
      # Rename the fragments with just the sample
      names(x) <- sample
      # Return as a list named for celltype
      output <- list(x)
      names(output) <- celltype
      output
    }
  )

  # Group frags by cell population
  renamedFrags <- unlist(renamedFrags, recursive = FALSE)
  splitFrags <- split(renamedFrags, f = names(renamedFrags))
  return(splitFrags)
}


# Tests if a string is a in the correct format to convert to GRanges
validRegionString <- function(regionString) {
  if (!is.character(regionString)) {
    return(FALSE)
  }

  pattern <- "([0-9]{1,2}|chr[0-9]{1,2}|chr[X-Y]{1,1}):[0-9]*-[0-9]*"
  matchedPattern <- stringr::str_extract(regionString, pattern)

  if (any(is.na(matchedPattern))) {
    return(FALSE)
  } else if (any(!matchedPattern == regionString)) {
    return(FALSE)
  }

  splits <- stringr::str_split(regionString, "[:-]")[[1]]
  start <- splits[2]
  end <- splits[3]
  if (any(start > end)) {
    return(FALSE)
  }

  # All conditions satisfied
  return(TRUE)
}

#' @title Convert a list of strings in the format "chr1:100-200" into a GRanges
#'
#' @description \code{StringsToGRanges} Turns a list of strings defining genomic 
#'  regions in the format chr1:100-200 into a GRanges object
#'
#' @param regionString A string or list of strings each in the format chr1:100-200
#' @return a GRanges object with ranges representing the input string(s)
#'
#' @export
#' @keywords utils
StringsToGRanges <- function(regionString) {
  # if (length(regionString)>1){
  #   boolList <- lapply(regionString, function(x){validRegionString(x)})
  #   if (!all(boolList)){
  #     stop("Some region strings are invalid. Given regions must all be strings matching format 'seqname:start-end', where start<end e.g. chr1:123000-123500")
  #   }
  # } else if(!validRegionString(regionString)) {
  #   stop("Region must be a string matching format 'seqname:start-end', where start<end e.g. chr1:123000-123500")
  # }
  . <- NULL
  regionString <- gsub(",", "", regionString)
  chrom <- gsub(":.*", "", regionString)
  startSite <- gsub(".*:", "", regionString) %>%
    gsub("-.*|-.*", "", .) %>%
    as.numeric()
  endSite <- gsub(".*-|.*-", "", regionString) %>% as.numeric()

  if (any(is.na(startSite) | is.na(endSite))) {
    stop(
      "Region must be a string matching format 'seqname:start-end', ",
      "where start < end (e.g. chr1:123000-123500)."
    )
  }

  if (any(startSite >= endSite)) {
    stop("Error in region string: Make sure the start of the genomic range occurs before the end")
  }
  regionGRanges <- GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::IRanges(start = startSite, end = endSite), strand = "*")
  return(regionGRanges)
}

#' @title Convert a GRanges object to a string in the format 'chr1:100-200'
#'
#' @description \code{GRangesToString} Turns a GRanges Object into
#'  a list of strings in the format chr1:100-200
#'
#' @param GR_obj the GRanges object to convert to a string
#' @return A string or list of strings in the format 'chr1:100-200' representing
#'  ranges in the input GRanges
#'
#' @export
#' @keywords utils
GRangesToString <- function(GR_obj) {
  paste(GenomicRanges::seqnames(GR_obj), ":", GenomicRanges::start(GR_obj), "-", GenomicRanges::end(GR_obj), sep = "")
}

#' @title Convert a data.frame or matrix to a GRanges
#'   
#' @param differentials a matrix/data.frame with a column tileColumn containing
#'   region strings in the format "chr:start-end"
#' @param tileColumn name of column containing region strings. Default is "Tile".
#'
#' @return a GRanges containing all original information
#' @export
#' @keywords utils
differentialsToGRanges <- function(differentials, tileColumn = "Tile") {
  regions <- MOCHA::StringsToGRanges(differentials[[tileColumn]])
  GenomicRanges::mcols(regions) <- differentials
  regions
}

# Helpers for getAnnotationDb
.TXDB_PREFIX <- "TxDb."
.has_TxDb_prefix <- function(x) {
  substr(x, 1L, nchar(.TXDB_PREFIX)) == .TXDB_PREFIX
}

.ORGDB_PREFIX <- "Org."
.has_OrgDb_prefix <- function(x) {
  substr(x, 1L, nchar(.ORGDB_PREFIX)) == .ORGDB_PREFIX
}

#' @title isMOCHAObject
#'
#' @description \code{isMOCHAObject} Checks with the object is a MOCHA object, and/or returns the type of object
#'
#' @param Object The output of callOpenTiles or getSampleTileMatrix. Both have a path to saved files. 
#' @param returnType A Boolean. If True, returns object type (OpenTiles, SampleTileMatrix, or NULL). If FALSE, returns True/FALSE 
#' @return A string or list of strings in the format 'chr1:100-200' representing
#'  ranges in the input GRanges
#'
#' @export
#' @keywords utils
isMOCHAObject <- function(Object, returnType = FALSE) {
    required_names <- c("summarizedData", "Genome", "TxDb", "OrgDb", "Directory", "History")
    if (inherits(Object, "MultiAssayExperiment")) {
        type = 'OpenTiles'
        metadata_names <- names(S4Vectors::metadata(Object))
    } else if (inherits(Object, "SummarizedExperiment")) {
        type = 'SampleTileMatrix'
        metadata_names <- names(S4Vectors::metadata(Object))
    }else if(!returnType){
        return(FALSE)
    }else{
        return(NULL)
    }

    if(all(required_names %in% metadata_names) & returnType){
        return(type)
    }else if(all(required_names %in% metadata_names)){
        return(TRUE)
    }
    return(FALSE)
}


#' @title updateDirectoryPath
#'
#' @description \code{updateDirectoryPath} Updated the path to save files, in case the directory moved.
#'
#' @param Object The output of callOpenTiles or getSampleTileMatrix. Both have a path to saved files. 
#' @param directoryPath A string, containing the absolute path to the directory. 
#' @return A string or list of strings in the format 'chr1:100-200' representing
#'  ranges in the input GRanges
#'
#' @export
#' @keywords utils
updateDirectoryPath <- function(Object, directoryPath = NULL) {

    if(!isMOCHAObject(Object)){
        warning("Object is not a MOCHA-related object. Returning NULL")
        return(NULL)
    }else if(!is.null(directoryPath) & dir.exists(directoryPath)){
        Object@metadata$Directory = directoryPath
        return(Object)
    }else if(is.null(directoryPath)) {
        warning("Directory path was not provided. No changes made")
        return(Object)
    }else {
         warning("Directory path not detected. No changes made")
        return(Object)
    }
}


#' @title Loads and attaches an installed TxDb or OrgDb-class Annotation database package.
#'
#' @description See \link[BSgenome]{getBSgenome}
#'
#' @param dbName Exact name of installed annotation data package.
#' @param type Expected class of the annotation data package, must be
#'   either "OrgDb" or "TxDb".
#' @return the loaded Annotation database object.#' @noRd
#' @keywords internal
getAnnotationDbFromInstalledPkgname <- function(dbName, type) {
  if (!methods::is(dbName, "character")) {
    stop(
      "dbName must be a character string. ",
      "Please provide TxDb, OrgDb, or BSgenome as a string."
    )
  }

  if (!type %in% c("OrgDb", "TxDb", 'BSgenome')) {
    stop('Invalid type. Type must be either "OrgDb", "TxDb", or "BSgenome".')
  }

  ok <- suppressWarnings(require(dbName,
    quietly = TRUE,
    character.only = TRUE
  ))

  if (!ok) {
    stop(
      "Package ", dbName, " is not available. Please install and provide ",
      "the name of a Bioconductor AnnotationData Package for your organism."
    )
  }

  pkgenvir <- as.environment(paste("package", dbName, sep = ":"))

  db <- try(get(dbName, envir = pkgenvir, inherits = FALSE), silent = TRUE)

  if (!methods::is(db, type)) {
    stop(dbName, " doesn't look like a valid ", type, " data package")
  }
  return(db)
}


#' @title Extract cell population names from a Tile Results or Sample Tile object.
#'
#' @description \code{getCellTypes} Returns a vector of cell names from a Tile Results or Sample Tile object.
#'
#' @param object tileResults object from callOpenTiles or SummarizedExperiment from getSampleTileMatrix
#' @return a vector of cell type names.
#' 
#' @export
#' @keywords utils
getCellTypes <- function(object) {
  if (methods::is(object, "MultiAssayExperiment")) {
    return(names(object))
  } else if (methods::is(object, "RangedSummarizedExperiment")) {
    return(names(SummarizedExperiment::assays(object)))
  } else {
    stop("Object not recognized. Please provide an object from callOpenTiles or getSampleTileMatrix.")
  }
}


#' @title Extract the GRanges for a particular cell population
#'
#' @description \code{getCellTypeTiles} Returns a GRanges object of all tiles called for a certain cell type
#'
#' @param object A SampleTileObject.
#' @param cellType A string describing one cell type.
#' @return a vector of cell type names.
#'
#' @export
#' @keywords utils
getCellTypeTiles <- function(object, cellType) {
  if (methods::is(object, "MultiAssayExperiment")) {
    stop("This is a MultiAssayExperiment, and thus like a tileResults object. Please provide a SampleTileMatrix object.")
  } else if (methods::is(object, "RangedSummarizedExperiment")) {
    all_ranges <- SummarizedExperiment::rowRanges(object)

    if (!all(cellType %in% SummarizedExperiment::assayNames(object))) {
      stop("Please provide cell types that are present in the SampleTileObject.")
    } else {

      peak_subSet = GenomicRanges::mcols(all_ranges)[, cellType, drop = FALSE]
      subRange <- all_ranges[apply(peak_subSet,1, any)]
    }

    return(subRange)
  } else {
    stop("Object not recognized. Please provide a SampleTileMatrix object (SummarizedExperiment).")
  }
}


#' @keywords internal
#' @noRd
.get_sample_celltype_count_tables <- function(object) {
  if (!is.null(object@metadata$summarizedData)) {
    summarizedData <- object@metadata$summarizedData
    assayNames <- SummarizedExperiment::assayNames(summarizedData)
    if (all(c("CellCounts", "FragmentCounts") %in% assayNames)) {
      return(list(
        CellCounts = as.data.frame(SummarizedExperiment::assays(summarizedData)[["CellCounts"]]),
        FragmentCounts = as.data.frame(SummarizedExperiment::assays(summarizedData)[["FragmentCounts"]])
      ))
    }
  }

  if (all(c("FragmentCounts", "CellCounts") %in% names(object@metadata))) {
    return(list(
      CellCounts = object@metadata$CellCounts,
      FragmentCounts = object@metadata$FragmentCounts
    ))
  }

  NULL
}

#' @title Extract Sample-celltype specific metadata
#' @description \code{getSampleCellTypeMetadata} Extract Sample-celltype
#'   specific metadata like fragment and cell counts
#'
#' @param object tileResults object from callOpenTiles or SummarizedExperiment
#'   from getSampleTileMatrix
#' @return a SummarizedExperiment where each assay is a different type of
#'   metadata.
#'
#' @export
#' @keywords utils
getSampleCellTypeMetadata <- function(object) {
  countTables <- .get_sample_celltype_count_tables(object)
  if (is.null(countTables)) {
    stop(
      "Object does not contain Sample-Celltype metadata. ",
      "Expected CellCounts and FragmentCounts in metadata$summarizedData ",
      "or legacy top-level metadata slots."
    )
  }

  sampleData <- SummarizedExperiment::colData(object)
  bioSamples <- rownames(sampleData)
  cellCounts <- countTables$CellCounts
  fragCounts <- countTables$FragmentCounts

  if (!all(bioSamples %in% colnames(cellCounts)) ||
      !all(bioSamples %in% colnames(fragCounts))) {
    stop(
      "Sample IDs in colData do not match columns in CellCounts/FragmentCounts. ",
      "Regenerate the MOCHA object or align sample metadata."
    )
  }

  # Count tables are cell population (row) by biological sample (column).
  cellMat <- as.matrix(cellCounts[, bioSamples, drop = FALSE])
  fragMat <- as.matrix(fragCounts[, bioSamples, drop = FALSE])

  metaList <- list(FragmentCounts = fragMat, CellCounts = cellMat)
  SE <- SummarizedExperiment::SummarizedExperiment(metaList, colData = sampleData)
  return(SE)
}


#' @title Plots the distribution of sample-tile intensities for a given cell
#'   population
#' @description \code{plotIntensityDistribution}  Plots the distribution of
#'   sample-tile intensities for a give cell population.
#'
#' @param TSAM_object  SummarizedExperiment from getSampleTileMatrix
#' @param cellPopulation Cell type names (assay name) within the TSAM_object
#' @param density Boolean to determine whether to plot density or histogram.
#'   Default is TRUE (plots density).
#' @param returnDF If TRUE, return the data frame without plotting. Default is
#'   FALSE.
#' @return data.frame or ggplot histogram.
#'
#' @export
#' @keywords plotting
plotIntensityDistribution <- function(TSAM_object, cellPopulation, returnDF = FALSE, density = TRUE) {
  Values <- NULL
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Package 'Biobase' is not installed. 'Biobase' is required for `plotIntensityDistribution`.")
  }
  mat <- unlist(Biobase::rowMedians(log2(getCellPopMatrix(TSAM_object, cellPopulation) + 1)))
  plotMat <- data.frame(Values = mat)

  if (returnDF) {
    return(plotMat)
  }
  if (density) {
    p1 <- ggplot2::ggplot(plotMat, ggplot2::aes(x = Values)) +
      ggplot2::geom_density() +
      ggplot2::theme_bw()
  } else {
    p1 <- ggplot2::ggplot(plotMat, ggplot2::aes(x = Values)) +
      ggplot2::geom_histogram() +
      ggplot2::theme_bw()
  }

  return(p1)
}


#' @title Add a column to the sample-level colData of a MOCHA object
#'
#' @description \code{addCellColData} adds a new column to the sample-level
#'   colData of a MOCHA tileResults (\code{MultiAssayExperiment} from
#'   \code{callOpenTiles}) or SampleTileMatrix
#'   (\code{RangedSummarizedExperiment} from \code{getSampleTileMatrix}).
#'   MOCHA pseudobulks by sample x cell-population, so colData rows are
#'   biological samples; the name mirrors \code{ArchR::addCellColData} for
#'   API familiarity.
#'
#' @param object A MOCHA tileResults or SampleTileMatrix object.
#' @param name Character scalar. Name of the column to add.
#' @param value Vector of values to add. Either length \code{nrow(colData)} (in
#'   which case it is assumed aligned with \code{samples}) or a named vector
#'   with names matching sample identifiers.
#' @param samples Optional character vector of sample identifiers that
#'   \code{value} corresponds to. Defaults to \code{rownames(colData(object))}.
#'   Missing samples will receive \code{NA}.
#' @param force Logical. If \code{TRUE}, an existing column with the same
#'   \code{name} is overwritten. Default \code{FALSE}.
#'
#' @return The input object with the new colData column attached.
#'
#' @export
#' @keywords utils
addCellColData <- function(object, name, value, samples = NULL, force = FALSE) {
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    stop("`name` must be a single non-empty character string.")
  }

  isMAE <- methods::is(object, "MultiAssayExperiment")
  isSE <- methods::is(object, "SummarizedExperiment")
  if (!isMAE && !isSE) {
    stop("`object` must be a MOCHA tileResults (MultiAssayExperiment) or SampleTileMatrix (SummarizedExperiment).")
  }

  cd <- if (isMAE) {
    MultiAssayExperiment::colData(object)
  } else {
    SummarizedExperiment::colData(object)
  }

  if (name %in% colnames(cd) && !force) {
    stop(sprintf("Column '%s' already exists in colData. Use force = TRUE to overwrite.", name))
  }

  rn <- rownames(cd)
  if (is.null(samples)) {
    if (length(value) != nrow(cd)) {
      stop(sprintf(
        "`value` has length %d but colData has %d rows. Provide `samples` to specify alignment or match the colData row count.",
        length(value), nrow(cd)
      ))
    }
    aligned <- value
  } else {
    if (length(samples) != length(value)) {
      stop("`samples` and `value` must have the same length.")
    }
    aligned <- rep(NA, nrow(cd))
    storage.mode(aligned) <- storage.mode(value)
    idx <- match(samples, rn)
    if (any(is.na(idx))) {
      stop(sprintf(
        "These samples were not found in colData: %s",
        paste(samples[is.na(idx)], collapse = ", ")
      ))
    }
    aligned[idx] <- value
  }

  cd[[name]] <- aligned
  if (isMAE) {
    MultiAssayExperiment::colData(object) <- cd
  } else {
    SummarizedExperiment::colData(object) <- cd
  }
  return(object)
}


#' @title Get per-cell-population open tiles from a MOCHA tileResults object
#'
#' @description \code{getOpenTiles} extracts the called open tiles (peaks) for
#'   one or more cell populations from a \code{MultiAssayExperiment} returned
#'   by \code{callOpenTiles}. By default tiles are returned as a
#'   \code{GRangesList} keyed by cell population; with \code{returnType =
#'   "data.frame"} they are flattened into a single data frame with a
#'   \code{CellPopulation} column.
#'
#' @param tileResults A \code{MultiAssayExperiment} from \code{callOpenTiles}.
#' @param cellPopulations Character vector of cell population names, or
#'   \code{"all"} (default) to return all populations.
#' @param returnType One of \code{"GRangesList"} (default) or
#'   \code{"data.frame"}.
#'
#' @return A \code{GRangesList} or \code{data.frame} of open tiles per cell
#'   population.
#'
#' @export
#' @keywords utils
getOpenTiles <- function(tileResults,
                         cellPopulations = "all",
                         returnType = c("GRangesList", "data.frame")) {
  if (!methods::is(tileResults, "MultiAssayExperiment")) {
    stop("`tileResults` must be a MultiAssayExperiment from callOpenTiles().")
  }
  returnType <- match.arg(returnType)

  available <- names(tileResults)
  if (length(cellPopulations) == 1L && tolower(cellPopulations) == "all") {
    cellPopulations <- available
  } else if (!all(cellPopulations %in% available)) {
    missing <- setdiff(cellPopulations, available)
    stop(sprintf(
      "These cell populations were not found in tileResults: %s",
      paste(missing, collapse = ", ")
    ))
  }

  grList <- lapply(cellPopulations, function(pop) {
    re <- tileResults[[pop]]
    peakMat <- RaggedExperiment::compactAssay(re, i = "peak")
    isPeak <- rowSums(peakMat == TRUE, na.rm = TRUE) > 0
    tiles <- SummarizedExperiment::rowRanges(re)
    if (length(tiles) != length(isPeak)) {
      keep <- seq_len(min(length(tiles), length(isPeak)))
      tiles <- tiles[keep]
      isPeak <- isPeak[keep]
    }
    tiles[isPeak]
  })
  names(grList) <- cellPopulations
  grList <- GenomicRanges::GRangesList(grList)

  if (returnType == "GRangesList") {
    return(grList)
  }

  dfs <- lapply(cellPopulations, function(pop) {
    gr <- grList[[pop]]
    if (length(gr) == 0L) {
      return(NULL)
    }
    data.frame(
      CellPopulation = pop,
      as.data.frame(gr),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, dfs)
}
