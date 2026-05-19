#' Require Seurat/Signac (Suggests) with an actionable error.
#' @noRd
.requireSeuratSignac <- function() {
  missing <- character(0)
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    missing <- c(missing, "Seurat")
  }
  if (!requireNamespace("Signac", quietly = TRUE)) {
    missing <- c(missing, "Signac")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    missing <- c(missing, "data.table")
  }
  if (length(missing) > 0) {
    stop(
      "The following packages are required for Seurat/Signac input: ",
      paste(missing, collapse = ", "),
      ". Install them with BiocManager::install() or remotes::install_github()."
    )
  }
  invisible(TRUE)
}

#' Resolve ChromatinAssay name on a Seurat object.
#' @noRd
.resolveSeuratAssay <- function(seuratObj, assay = NULL) {
  if (is.null(assay) || !nzchar(assay)) {
    assay <- Seurat::DefaultAssay(seuratObj)
  }
  available <- Seurat::Assays(seuratObj)
  if (!assay %in% available) {
    stop(
      "Assay '", assay, "' not found. Available assays: ",
      paste(available, collapse = ", ")
    )
  }
  chromAssays <- vapply(
    available,
    function(a) {
      methods::is(seuratObj[[a]], "ChromatinAssay")
    },
    logical(1)
  )
  if (!methods::is(seuratObj[[assay]], "ChromatinAssay")) {
    chromNames <- available[chromAssays]
    stop(
      "Assay '", assay, "' is not a ChromatinAssay. ",
      if (length(chromNames) > 0) {
        paste0("ChromatinAssays in this object: ", paste(chromNames, collapse = ", "))
      } else {
        "No ChromatinAssay found in this Seurat object."
      }
    )
  }
  assay
}

#' Map Fragment objects to sample IDs.
#' @noRd
.mapFragmentsToSamples <- function(fragObjs, meta, sampleColumn, cellIds) {
  samples <- unique(as.character(meta[[sampleColumn]]))
  nSamples <- length(samples)
  nFrags <- length(fragObjs)

  if (nFrags == 0) {
    stop("No Fragment objects found in the Signac assay.")
  }

  fragMap <- stats::setNames(vector("list", nSamples), samples)

  if (nFrags == 1L && nSamples >= 1L) {
    for (s in samples) {
      fragMap[[s]] <- fragObjs[[1L]]
    }
    return(fragMap)
  }

  if (nFrags == nSamples) {
    overlapMat <- matrix(
      0L,
      nrow = nFrags,
      ncol = nSamples,
      dimnames = list(names(fragObjs), samples)
    )
    if (is.null(rownames(overlapMat))) {
      rownames(overlapMat) <- paste0("Fragment", seq_len(nFrags))
    }

    for (i in seq_len(nFrags)) {
      fo <- fragObjs[[i]]
      cellsMap <- methods::slot(fo, "cells")
      if (length(cellsMap) == 0) {
        next
      }
      seuratIdsInFrag <- names(cellsMap)
      for (s in samples) {
        sampleCells <- cellIds[meta[[sampleColumn]] == s]
        overlapMat[i, s] <- length(intersect(seuratIdsInFrag, sampleCells))
      }
    }

    assigned <- character(nFrags)
    for (i in seq_len(nFrags)) {
      best <- which.max(overlapMat[i, , drop = TRUE])
      if (overlapMat[i, best] == 0L) {
        stop(
          "Could not map Fragment object ", rownames(overlapMat)[i],
          " to any sample. Check Fragment@cells and metadata[[", sampleColumn, "]]."
        )
      }
      assigned[i] <- samples[best]
    }

    dupSamples <- assigned[duplicated(assigned)]
    if (length(dupSamples) > 0) {
      stop(
        "Multiple Fragment objects map to the same sample(s): ",
        paste(unique(dupSamples), collapse = ", "),
        ". Provide fragmentPathColumn or one Fragment per sample."
      )
    }

    for (i in seq_len(nFrags)) {
      fragMap[[assigned[i]]] <- fragObjs[[i]]
    }
    return(fragMap)
  }

  stop(
    "Found ", nFrags, " Fragment object(s) for ", nSamples,
    " sample(s). Use one Fragment per sample, a single combined Fragment file, ",
    "or provide fragmentPathColumn in metadata."
  )
}

#' Read fragments for one sample from a Signac Fragment object.
#' @noRd
.readSampleFragmentsFromFragment <- function(
    fragObj,
    sampleCells,
    cellCol,
    verbose = FALSE) {
  cellsMap <- methods::slot(fragObj, "cells")
  if (length(cellsMap) == 0) {
    stop(
      "Fragment object has empty @cells slot. ",
      "Re-create the ChromatinAssay with CreateChromatinAssay() and valid cells."
    )
  }

  seuratIdsInFrag <- names(cellsMap)
  if (length(intersect(seuratIdsInFrag, sampleCells)) == 0) {
    stop(
      "Barcode prefix mismatch between Seurat cells and Signac Fragments. ",
      "No Fragment@cells entries match sample cell IDs. ",
      "Check merge prefixes and sample assignments."
    )
  }

  fileBcs <- unname(cellsMap)[names(cellsMap) %in% sampleCells]
  if (length(fileBcs) == 0) {
    return(.emptyFragmentGRanges(cellCol))
  }

  path <- Signac::GetFragmentData(fragObj, slot = "path")
  if (!file.exists(path)) {
    stop("Fragment file not found: ", path)
  }

  if (verbose) {
    message("Reading fragments from ", path)
  }

  df <- .readFragmentFile(path, barcodes = fileBcs)
  if (nrow(df) == 0) {
    return(.emptyFragmentGRanges(cellCol))
  }

  df$barcode <- names(cellsMap)[match(df$barcode, cellsMap)]  # file barcode -> Seurat cell ID

  gr <- GenomicRanges::GRanges(
    seqnames = df$chr,
    ranges = IRanges::IRanges(start = df$start, end = df$end),
    strand = "*"
  )
  GenomicRanges::mcols(gr)[[cellCol]] <- df$barcode
  gr
}

#' Resolve file barcodes to read for a sample (Fragment@cells aware).
#' @noRd
.file_barcodes_for_sample <- function(sampleCells, cellsMap = NULL) {
  if (!is.null(cellsMap) && length(cellsMap) > 0) {
    # cellsMap: names = Seurat cell IDs, values = file barcodes
    file_bcs <- unname(cellsMap)[names(cellsMap) %in% sampleCells]
    if (length(file_bcs) > 0) {
      return(list(file_barcodes = file_bcs, cells_map = cellsMap))
    }
  }
  list(file_barcodes = sampleCells, cells_map = NULL)
}

#' Read a fragment file and filter to requested barcodes.
#' @noRd
.readFragmentFile <- function(path, barcodes) {
  barcodes <- unique(as.character(barcodes))
  if (length(barcodes) == 0) {
    return(data.table::data.table(
      chr = character(),
      start = integer(),
      end = integer(),
      barcode = character(),
      count = integer()
    ))
  }

  colNames <- c("chr", "start", "end", "barcode", "count")
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    df <- data.table::fread(
      cmd = paste("zcat", shQuote(path)),
      header = FALSE,
      col.names = colNames,
      showProgress = FALSE
    )
  } else {
    df <- data.table::fread(
      file = path,
      header = FALSE,
      col.names = colNames,
      showProgress = FALSE
    )
  }

  df <- df[barcode %in% barcodes]
  df[, `:=`(start = as.integer(start), end = as.integer(end))]
  df
}

#' Empty GRanges with the required cell column.
#' @noRd
.emptyFragmentGRanges <- function(cellCol) {
  gr <- GenomicRanges::GRanges()
  mcols(gr) <- S4Vectors::DataFrame(
    stats::setNames(list(character()), cellCol)
  )
  gr
}

#' Convert Seurat/Signac object to sample-level MOCHA fragment inputs.
#'
#' @description Converts a Seurat object with a Signac \code{ChromatinAssay}
#'   into a sample-level \code{GRangesList} and matching \code{cellColData} for
#'   use with \code{\link{callOpenTiles}}.
#'
#' @param seuratObj A Seurat object with a ChromatinAssay.
#' @param assay Name of the ChromatinAssay. Default uses
#'   \code{Seurat::DefaultAssay(seuratObj)}.
#' @param sampleColumn Metadata column with biological sample IDs.
#' @param cellPopLabel Metadata column with cell population labels.
#' @param cellCol Column name for cell barcodes in fragment \code{mcols}.
#'   Default is \code{"RG"} (MOCHA/ArchR convention).
#' @param fragmentPathColumn Optional metadata column with per-cell fragment
#'   file paths. If provided, takes precedence over Signac \code{Fragments()}.
#'   File barcodes must match \code{rownames(meta)} unless \code{Fragment@cells}
#'   is available on the assay (then file barcodes are translated automatically).
#' @param dropEmptySamples If \code{TRUE}, drop samples with no fragments and
#'   warn.
#' @param verbose Print progress messages.
#'
#' @return A list with \code{ATACFragments} (\code{GRangesList}) and
#'   \code{cellColData} (\code{data.frame}).
#' @export
#' @keywords utils
seuratToMOCHAInputs <- function(
    seuratObj,
    assay = NULL,
    sampleColumn = "Sample",
    cellPopLabel,
    cellCol = "RG",
    fragmentPathColumn = NULL,
    dropEmptySamples = TRUE,
    verbose = FALSE) {
  .seuratToSampleLevelFragments(
    seuratObj = seuratObj,
    assay = assay,
    sampleColumn = sampleColumn,
    cellPopLabel = cellPopLabel,
    cellCol = cellCol,
    fragmentPathColumn = fragmentPathColumn,
    dropEmptySamples = dropEmptySamples,
    verbose = verbose
  )
}

#' Internal: Seurat/Signac -> sample-level GRangesList + metadata.
#' @noRd
.seuratToSampleLevelFragments <- function(
    seuratObj,
    assay = NULL,
    sampleColumn = "Sample",
    cellPopLabel,
    cellCol = "RG",
    fragmentPathColumn = NULL,
    dropEmptySamples = TRUE,
    verbose = FALSE) {
  .requireSeuratSignac()

  if (!methods::is(seuratObj, "Seurat")) {
    stop("seuratObj must be a Seurat object.")
  }

  assay <- .resolveSeuratAssay(seuratObj, assay = assay)
  meta <- as.data.frame(seuratObj@meta.data)

  if (!(sampleColumn %in% colnames(meta))) {
    stop(
      "sampleColumn '", sampleColumn, "' not found in Seurat metadata. ",
      "Available columns: ", paste(colnames(meta), collapse = ", ")
    )
  }
  if (!(cellPopLabel %in% colnames(meta))) {
    stop(
      "cellPopLabel '", cellPopLabel, "' not found in Seurat metadata. ",
      "Available columns: ", paste(colnames(meta), collapse = ", ")
    )
  }

  cellIds <- rownames(meta)
  if (is.null(cellIds) || !any(nzchar(cellIds))) {
    if (cellCol %in% colnames(meta)) {
      cellIds <- as.character(meta[[cellCol]])
    } else {
      stop(
        "Seurat metadata must have rownames matching cell barcodes, or include ",
        "column '", cellCol, "'."
      )
    }
  }
  rownames(meta) <- cellIds

  naSample <- is.na(meta[[sampleColumn]])
  naPop <- is.na(meta[[cellPopLabel]])
  if (any(naSample | naPop)) {
    nDrop <- sum(naSample | naPop)
    if (verbose) {
      warning(
        "Dropping ", nDrop, " cell(s) with NA in ", sampleColumn,
        " or ", cellPopLabel, "."
      )
    }
    keep <- !(naSample | naPop)
    meta <- meta[keep, , drop = FALSE]
    cellIds <- cellIds[keep]
  }

  samples <- unique(as.character(meta[[sampleColumn]]))
  out <- stats::setNames(vector("list", length(samples)), samples)

  if (!is.null(fragmentPathColumn)) {
    if (!(fragmentPathColumn %in% colnames(meta))) {
      stop(
        "fragmentPathColumn '", fragmentPathColumn, "' not found in metadata."
      )
    }
    pathMap <- lapply(samples, function(s) {
      unique(as.character(meta[[fragmentPathColumn]][meta[[sampleColumn]] == s]))
    })
    bad <- samples[vapply(pathMap, length, integer(1)) != 1L]
    if (length(bad) > 0) {
      stop(
        "Each sample must have exactly one unique path in ",
        fragmentPathColumn, ". Check sample(s): ",
        paste(bad, collapse = ", ")
      )
    }
    names(pathMap) <- samples

    frag_objs <- tryCatch(
      Signac::Fragments(seuratObj[[assay]]),
      error = function(e) NULL
    )
    default_cells_map <- NULL
    if (!is.null(frag_objs) && length(frag_objs) > 0) {
      default_cells_map <- methods::slot(frag_objs[[1L]], "cells")
    }

    for (s in samples) {
      sampleCells <- cellIds[meta[[sampleColumn]] == s]
      path <- pathMap[[s]][[1L]]
      if (!file.exists(path)) {
        stop("Fragment file not found for sample '", s, "': ", path)
      }
      if (verbose) {
        message("Reading fragments for sample ", s, " from ", path)
      }
      bc_info <- .file_barcodes_for_sample(sampleCells, default_cells_map)
      df <- .readFragmentFile(path, barcodes = bc_info$file_barcodes)
      if (nrow(df) == 0) {
        out[[s]] <- .emptyFragmentGRanges(cellCol)
        next
      }
      if (!is.null(bc_info$cells_map)) {
        df$barcode <- names(bc_info$cells_map)[match(df$barcode, bc_info$cells_map)]
      }
      gr <- GenomicRanges::GRanges(
        seqnames = df$chr,
        ranges = IRanges::IRanges(start = df$start, end = df$end),
        strand = "*"
      )
      GenomicRanges::mcols(gr)[[cellCol]] <- df$barcode
      out[[s]] <- gr
    }
  } else {
    fragObjs <- Signac::Fragments(seuratObj[[assay]])
    if (is.null(fragObjs) || length(fragObjs) == 0) {
      stop(
        "No Fragment objects in assay '", assay, "'. ",
        "Add fragments with CreateChromatinAssay() or provide fragmentPathColumn."
      )
    }

    fragMap <- .mapFragmentsToSamples(
      fragObjs,
      meta = meta,
      sampleColumn = sampleColumn,
      cellIds = cellIds
    )

    for (s in samples) {
      fo <- fragMap[[s]]
      if (is.null(fo)) {
        stop("No Fragment object mapped to sample '", s, "'.")
      }
      sampleCells <- cellIds[meta[[sampleColumn]] == s]
      out[[s]] <- .readSampleFragmentsFromFragment(
        fragObj = fo,
        sampleCells = sampleCells,
        cellCol = cellCol,
        verbose = verbose
      )
    }
  }

  emptySamples <- samples[lengths(out) == 0L]
  if (length(emptySamples) > 0 && dropEmptySamples) {
    warning(
      "Dropping sample(s) with no fragments: ",
      paste(emptySamples, collapse = ", ")
    )
    out <- out[lengths(out) > 0L]
    samples <- names(out)
  }

  if (length(out) == 0) {
    stop("No fragments remained after reading Seurat/Signac input.")
  }

  missingMcol <- vapply(out, function(gr) {
    !(cellCol %in% colnames(GenomicRanges::mcols(gr)))
  }, logical(1))
  if (any(missingMcol)) {
    stop(
      "Fragments missing column '", cellCol, "' for sample(s): ",
      paste(names(out)[missingMcol], collapse = ", ")
    )
  }

  list(
    ATACFragments = GenomicRanges::GRangesList(out),
    cellColData = meta
  )
}

#' @rdname callOpenTiles-methods
#' @aliases callOpenTiles,Seurat-method
#' @noRd
.callOpenTiles_Seurat <- function(
    ATACFragments,
    cellColData = NULL,
    blackList,
    genome,
    cellPopLabel,
    cellPopulations = "ALL",
    sampleColumn = "Sample",
    studySignal = NULL,
    generalizeStudySignal = FALSE,
    cellCol = "RG",
    TxDb,
    OrgDb,
    outDir,
    numCores = 30,
    verbose = FALSE,
    force = FALSE,
    peakModel = NULL,
    returnClass = c("legacy", "mocha")) {
  returnClass <- match.arg(returnClass)
  if (!methods::is(blackList, "GRanges")) {
    stop("Invalid blackList. blackList must be a GRanges.")
  }
  if (missing(genome) || is.null(genome)) {
    stop(
      "genome must be provided for Seurat input (e.g. 'hg38' or a BSgenome package name)."
    )
  }
  if (missing(outDir) || is.null(outDir) || !nzchar(outDir)) {
    stop("outDir must be a non-empty directory path for Seurat input.")
  }

  parsed <- .seuratToSampleLevelFragments(
    seuratObj = ATACFragments,
    assay = NULL,
    sampleColumn = sampleColumn,
    cellPopLabel = cellPopLabel,
    cellCol = cellCol,
    fragmentPathColumn = NULL,
    dropEmptySamples = TRUE,
    verbose = verbose
  )

  if (!is.null(cellColData)) {
    parsed$cellColData <- cellColData
  }

  if (is.null(studySignal) && !isTRUE(generalizeStudySignal)) {
    if (!("nFrags" %in% colnames(parsed$cellColData))) {
      if ("nCount_ATAC" %in% colnames(parsed$cellColData)) {
        parsed$cellColData$nFrags <- parsed$cellColData$nCount_ATAC
      }
    }
  }

  .callOpenTiles_default(
    ATACFragments = parsed$ATACFragments,
    cellColData = parsed$cellColData,
    blackList = blackList,
    genome = genome,
    cellPopLabel = cellPopLabel,
    cellPopulations = cellPopulations,
    sampleColumn = sampleColumn,
    studySignal = studySignal,
    generalizeStudySignal = generalizeStudySignal,
    cellCol = cellCol,
    TxDb = TxDb,
    OrgDb = OrgDb,
    outDir = outDir,
    numCores = numCores,
    verbose = verbose,
    force = force,
    peakModel = peakModel,
    returnClass = returnClass
  )
}
