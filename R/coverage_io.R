#' @keywords internal
#' @noRd

.require_rtracklayer <- function() {
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop(
      "Package 'rtracklayer' is required to read structured coverage tracks. ",
      "Please install 'rtracklayer' to proceed."
    )
  }
}

.sanitize_track_filename <- function(sampleId) {
  x <- gsub("%", "%25", sampleId, fixed = TRUE)
  x <- gsub(" ", "%20", x, fixed = TRUE)
  x <- gsub(".", "%2E", x, fixed = TRUE)
  x
}

.unsanitize_track_filename <- function(safeName) {
  x <- gsub("%2E", ".", safeName, fixed = TRUE)
  x <- gsub("%20", " ", x, fixed = TRUE)
  x <- gsub("%25", "%", x, fixed = TRUE)
  x
}

.validate_track_path_label <- function(label, labelType = c("sample", "cellPop")) {
  labelType <- match.arg(labelType)
  if (length(label) != 1L || is.na(label) || !nzchar(label)) {
    stop(labelType, " label must be a non-empty string.")
  }
  if (label %in% c(".", "..")) {
    stop("Invalid ", labelType, " label '", label, "'. Please rename.")
  }
  if (grepl("[/\\\\:*?\"<>|\\x00]", label, perl = TRUE)) {
    stop(
      "Invalid ", labelType, " label '", label,
      "': contains characters not allowed in track paths ",
      "(/ \\ : * ? \" < > |). Please rename."
    )
  }
  invisible(label)
}

.validate_track_sample_ids <- function(sampleIds) {
  vapply(sampleIds, .validate_track_path_label, FUN.VALUE = character(1L), labelType = "sample")
  safeIds <- vapply(sampleIds, .sanitize_track_filename, character(1L))
  if (any(duplicated(safeIds))) {
    dupSafe <- unique(safeIds[duplicated(safeIds)])
    offenders <- unique(sampleIds[safeIds %in% dupSafe])
    stop(
      "Sample IDs collide after filename encoding: ",
      paste(offenders, collapse = ", "),
      ". Please rename samples so encoded track filenames are unique."
    )
  }
  invisible(safeIds)
}

.mocha_coverage_cache <- new.env(parent = emptyenv())

.clear_coverage_cache <- function() {
  rm(list = ls(envir = .mocha_coverage_cache), envir = .mocha_coverage_cache)
  invisible(NULL)
}

.coverage_cache_key <- function(outDir, cellPop, trackFolder, samples) {
  sampleKey <- if (is.null(samples)) {
    ""
  } else {
    paste(samples, collapse = "\001")
  }
  paste(
    normalizePath(outDir, winslash = "/"),
    cellPop,
    trackFolder,
    sampleKey,
    sep = "::"
  )
}

.has_structured_coverage <- function(outDir, cellPop) {
  structuredDir <- .coverage_cell_pop_dir(outDir, cellPop)
  dir.exists(structuredDir) && (
    dir.exists(.coverage_track_type_dir(outDir, cellPop, "Accessibility")) ||
      dir.exists(.coverage_track_type_dir(outDir, cellPop, "Insertions"))
  )
}

.has_legacy_coverage <- function(outDir, cellPop) {
  file.exists(.coverage_legacy_path(outDir, cellPop))
}

.mocha_warned_mixed_layout <- new.env(parent = emptyenv())

.warn_mixed_coverage_layout <- function(outDir, cellPop) {
  if (.has_structured_coverage(outDir, cellPop) && .has_legacy_coverage(outDir, cellPop)) {
    warnKey <- paste0(normalizePath(outDir, winslash = "/"), "::", cellPop)
    if (!exists(warnKey, envir = .mocha_warned_mixed_layout, inherits = FALSE)) {
      assign(warnKey, TRUE, envir = .mocha_warned_mixed_layout)
      warning(
        "Both structured tracks and legacy RDS coverage exist for '",
        cellPop,
        "'. Using structured tracks. Consider removing ",
        .coverage_legacy_path(outDir, cellPop),
        " or running migrateCoverageToTracks(removeLegacy = TRUE).",
        call. = FALSE
      )
    }
  }
}

.expected_coverage_sample_ids <- function(outDir, cellPop, sampleIds = NULL) {
  if (!is.null(sampleIds)) {
    return(sampleIds)
  }
  manifest <- .read_coverage_manifest(outDir)
  if (!is.null(manifest$cellPopulationSamples[[cellPop]])) {
    return(manifest$cellPopulationSamples[[cellPop]])
  }
  NULL
}

.structured_track_files_complete <- function(outDir, cellPop, sampleIds) {
  trackTypes <- c("Accessibility", "Insertions")
  for (trackType in trackTypes) {
    for (sampleId in sampleIds) {
      trackPath <- .coverage_track_file_path(outDir, cellPop, trackType, sampleId)
      if (!file.exists(trackPath)) {
        return(FALSE)
      }
    }
  }
  TRUE
}

.structured_tracks_complete <- function(outDir, cellPop, sampleIds = NULL) {
  if (is.null(sampleIds) || length(sampleIds) == 0) {
    sampleIds <- .expected_coverage_sample_ids(outDir, cellPop)
  }
  if (is.null(sampleIds) || length(sampleIds) == 0) {
    return(FALSE)
  }
  .structured_track_files_complete(outDir, cellPop, sampleIds)
}

.prune_orphaned_track_files <- function(outDir, cellPop, sampleIds, verbose = FALSE) {
  expectedSafe <- vapply(sampleIds, .sanitize_track_filename, character(1L))
  trackTypes <- c("Accessibility", "Insertions")
  removed <- character(0)
  for (trackType in trackTypes) {
    typeDir <- .coverage_track_type_dir(outDir, cellPop, trackType)
    if (!dir.exists(typeDir)) {
      next
    }
    bwFiles <- list.files(typeDir, pattern = "\\.bw$", full.names = TRUE)
    for (bwFile in bwFiles) {
      safeName <- sub("\\.bw$", "", basename(bwFile))
      if (!safeName %in% expectedSafe) {
        file.remove(bwFile)
        removed <- c(removed, bwFile)
      }
    }
  }
  if (verbose && length(removed) > 0) {
    message(
      "Removed ", length(removed), " orphaned track file(s) for cell population ",
      cellPop, "."
    )
  }
  invisible(removed)
}

.coverage_layout_dir <- function(outDir) {
  file.path(outDir, "tracks")
}

.coverage_manifest_path <- function(outDir) {
  file.path(.coverage_layout_dir(outDir), "manifest.rds")
}

.coverage_cell_pop_dir <- function(outDir, cellPop) {
  file.path(.coverage_layout_dir(outDir), cellPop)
}

.coverage_track_type_dir <- function(outDir, cellPop, trackType) {
  file.path(.coverage_cell_pop_dir(outDir, cellPop), trackType)
}

.coverage_track_file_path <- function(outDir, cellPop, trackType, sampleId) {
  safeId <- .sanitize_track_filename(sampleId)
  file.path(
    .coverage_track_type_dir(outDir, cellPop, trackType),
    paste0(safeId, ".bw")
  )
}

.coverage_legacy_path <- function(outDir, cellPop) {
  file.path(outDir, paste0(cellPop, "_CoverageFiles.RDS"))
}

.detect_coverage_layout <- function(outDir, cellPop) {
  if (.has_structured_coverage(outDir, cellPop)) {
    return("structured")
  }
  if (.has_legacy_coverage(outDir, cellPop)) {
    return("legacy")
  }
  "missing"
}

.coverage_tracks_exist <- function(outDir, cellPop, sampleIds = NULL) {
  layout <- .detect_coverage_layout(outDir, cellPop)
  if (layout == "legacy") {
    return(TRUE)
  }
  if (layout == "missing") {
    return(FALSE)
  }

  expectedSamples <- .expected_coverage_sample_ids(outDir, cellPop, sampleIds)
  if (is.null(expectedSamples) || length(expectedSamples) == 0) {
    accDir <- .coverage_track_type_dir(outDir, cellPop, "Accessibility")
    insDir <- .coverage_track_type_dir(outDir, cellPop, "Insertions")
    return(
      length(list.files(accDir, pattern = "\\.bw$")) > 0 &&
        length(list.files(insDir, pattern = "\\.bw$")) > 0
    )
  }

  .structured_track_files_complete(outDir, cellPop, expectedSamples)
}

.read_coverage_manifest <- function(outDir) {
  path <- .coverage_manifest_path(outDir)
  if (file.exists(path)) {
    return(readRDS(path))
  }
  list(
    version = 1L,
    format = "bigwig",
    sample_name_map = list(),
    cellPopulations = character(0)
  )
}

.write_coverage_manifest <- function(outDir, manifest) {
  tracksDir <- .coverage_layout_dir(outDir)
  if (!dir.exists(tracksDir)) {
    dir.create(tracksDir, recursive = TRUE, showWarnings = FALSE)
  }
  saveRDS(manifest, .coverage_manifest_path(outDir))
}

.update_coverage_manifest <- function(outDir, cellPop, sampleIds) {
  .validate_track_sample_ids(sampleIds)
  manifest <- .read_coverage_manifest(outDir)
  for (sid in sampleIds) {
    safeId <- .sanitize_track_filename(sid)
    manifest$sample_name_map[[safeId]] <- sid
  }
  if (is.null(manifest$cellPopulationSamples)) {
    manifest$cellPopulationSamples <- list()
  }
  manifest$cellPopulationSamples[[cellPop]] <- sampleIds
  if (!cellPop %in% manifest$cellPopulations) {
    manifest$cellPopulations <- c(manifest$cellPopulations, cellPop)
  }
  manifest$version <- 1L
  manifest$format <- "bigwig"
  .write_coverage_manifest(outDir, manifest)
  invisible(manifest)
}

.restore_sample_name <- function(safeName, manifest) {
  mapped <- manifest$sample_name_map[[safeName]]
  if (!is.null(mapped)) {
    return(mapped)
  }
  decoded <- .unsanitize_track_filename(safeName)
  if (!identical(decoded, safeName)) {
    return(decoded)
  }
  safeName
}

.normalize_track_granges_score <- function(gr) {
  if (is.null(gr) || length(gr) == 0) {
    return(gr)
  }
  mcols <- S4Vectors::mcols(gr)
  if ("score" %in% colnames(mcols)) {
    return(gr)
  }
  scoreCol <- if ("Score" %in% colnames(mcols)) {
    "Score"
  } else if (ncol(mcols) >= 1) {
    colnames(mcols)[1]
  } else {
    return(gr)
  }
  S4Vectors::mcols(gr)$score <- mcols[[scoreCol]]
  gr
}

.ensure_bigwig_seqinfo <- function(gr) {
  if (length(gr) == 0) {
    return(gr)
  }
  chr <- as.character(GenomicRanges::seqnames(gr))
  maxEnds <- stats::setNames(
    as.integer(vapply(split(GenomicRanges::end(gr), chr), max, numeric(1L))),
    unique(chr)
  )
  sl <- GenomeInfoDb::seqlengths(gr)
  for (nm in names(maxEnds)) {
    if (is.na(sl[nm]) || sl[nm] < maxEnds[[nm]]) {
      sl[nm] <- maxEnds[[nm]]
    }
  }
  GenomeInfoDb::seqlengths(gr) <- sl
  gr
}

.filter_coverage_granges_list <- function(grList, samples = NULL, region = NULL) {
  if (is.null(grList)) {
    return(list())
  }
  if (methods::is(grList, "GRanges")) {
    grList <- stats::setNames(list(grList), "sample1")
  }
  if (!is.null(samples)) {
    keep <- intersect(names(grList), samples)
    grList <- grList[keep]
  }
  if (!is.null(region) && length(grList) > 0) {
    grList <- lapply(grList, function(gr) {
      if (length(gr) == 0) {
        return(gr)
      }
      plyranges::join_overlap_intersect(gr, region)
    })
  }
  grList
}

.writeCoverageTracks <- function(covFiles,
                                 outDir,
                                 cellPop,
                                 force = FALSE,
                                 verbose = FALSE) {
  trackTypes <- c("Accessibility", "Insertions")
  missingTypes <- setdiff(trackTypes, names(covFiles))
  if (length(missingTypes) > 0) {
    stop(
      "Coverage bundle must contain ",
      paste(trackTypes, collapse = " and "),
      ". Missing: ",
      paste(missingTypes, collapse = ", ")
    )
  }

  sampleIds <- names(covFiles[["Accessibility"]])
  if (length(sampleIds) == 0) {
    stop("No samples found in coverage bundle for cell population ", cellPop, ".")
  }
  .validate_track_path_label(cellPop, labelType = "cellPop")
  .validate_track_sample_ids(sampleIds)
  .prune_orphaned_track_files(outDir, cellPop, sampleIds, verbose = verbose)

  for (trackType in trackTypes) {
    typeDir <- .coverage_track_type_dir(outDir, cellPop, trackType)
    if (!dir.exists(typeDir)) {
      dir.create(typeDir, recursive = TRUE, showWarnings = FALSE)
    }
    grList <- covFiles[[trackType]]
    for (sampleId in names(grList)) {
      outFile <- .coverage_track_file_path(outDir, cellPop, trackType, sampleId)
      if (!force && file.exists(outFile)) {
        next
      }
      if (verbose) {
        message("Writing ", outFile)
      }
      grOut <- .ensure_bigwig_seqinfo(grList[[sampleId]])
      plyranges::write_bigwig(grOut, outFile)
    }
  }

  .update_coverage_manifest(outDir, cellPop, sampleIds)
  .clear_coverage_cache()
  invisible(TRUE)
}

.loadCoverageTracks <- function(outDir,
                                cellPop,
                                trackType = c("accessibility", "insertions"),
                                samples = NULL,
                                region = NULL,
                                verbose = FALSE) {
  trackType <- match.arg(trackType)
  layout <- .detect_coverage_layout(outDir, cellPop)
  trackFolder <- if (trackType == "accessibility") {
    "Accessibility"
  } else {
    "Insertions"
  }

  if (layout == "missing") {
    stop(
      "Coverage tracks for cell population '", cellPop, "' could not be found. ",
      "Expected structured layout under ",
      .coverage_cell_pop_dir(outDir, cellPop),
      " or legacy file ",
      .coverage_legacy_path(outDir, cellPop),
      "."
    )
  }

  if (layout == "legacy") {
    bundle <- readRDS(.coverage_legacy_path(outDir, cellPop))
    grList <- .selectCoverageFromBundle(bundle, coverage = trackType == "accessibility")
    return(.filter_coverage_granges_list(grList, samples = samples, region = region))
  }

  .warn_mixed_coverage_layout(outDir, cellPop)
  .require_rtracklayer()

  if (is.null(region)) {
    cacheKey <- .coverage_cache_key(outDir, cellPop, trackFolder, samples)
    if (exists(cacheKey, envir = .mocha_coverage_cache, inherits = FALSE)) {
      return(get(cacheKey, envir = .mocha_coverage_cache))
    }
  }

  manifest <- .read_coverage_manifest(outDir)
  trackDir <- .coverage_track_type_dir(outDir, cellPop, trackFolder)
  if (!dir.exists(trackDir)) {
    stop(
      "Structured coverage directory not found: ",
      trackDir
    )
  }

  bwFiles <- list.files(trackDir, pattern = "\\.bw$", full.names = TRUE)
  grList <- list()
  for (bwFile in bwFiles) {
    safeName <- sub("\\.bw$", "", basename(bwFile))
    sampleName <- .restore_sample_name(safeName, manifest)
    if (!is.null(samples) && !sampleName %in% samples) {
      next
    }
    if (verbose) {
      message("Reading ", bwFile)
    }
    if (!is.null(region)) {
      gr <- rtracklayer::import.bw(bwFile, format = "BigWig", which = region)
    } else {
      gr <- rtracklayer::import.bw(bwFile, format = "BigWig")
    }
    grList[[sampleName]] <- .normalize_track_granges_score(gr)
  }

  if (is.null(region)) {
    cacheKey <- .coverage_cache_key(outDir, cellPop, trackFolder, samples)
    assign(cacheKey, grList, envir = .mocha_coverage_cache)
  }

  grList
}

.copy_structured_accessibility_tracks <- function(outDir,
                                                   cellPop,
                                                   sampleIds,
                                                   destDir,
                                                   verbose = FALSE) {
  if (!dir.exists(destDir)) {
    dir.create(destDir, recursive = TRUE, showWarnings = FALSE)
  }
  exported <- list()
  for (sample in sampleIds) {
    src <- .coverage_track_file_path(outDir, cellPop, "Accessibility", sample)
    fileName <- gsub(" ", "__", paste(cellPop, gsub("\\.", "__", sample), sep = "__"))
    dest <- file.path(destDir, paste0(fileName, "_Coverage.bw"))
    if (verbose) {
      message("Copying ", src, " to ", dest)
    }
    file.copy(src, dest, overwrite = TRUE)
    exported[[sample]] <- dest
  }
  names(exported) <- paste(cellPop, gsub("\\.", "__", names(exported)), sep = "__")
  exported
}

.readCoverageBundle <- function(outDir,
                                cellPop,
                                samples = NULL,
                                region = NULL,
                                verbose = FALSE) {
  layout <- .detect_coverage_layout(outDir, cellPop)
  if (layout == "missing") {
    stop(
      "Coverage tracks for cell population '", cellPop, "' could not be found. ",
      "Expected structured layout under ",
      .coverage_cell_pop_dir(outDir, cellPop),
      " or legacy file ",
      .coverage_legacy_path(outDir, cellPop),
      "."
    )
  }

  if (layout == "legacy") {
    bundle <- readRDS(.coverage_legacy_path(outDir, cellPop))
    if (is.null(samples) && is.null(region)) {
      return(bundle)
    }
    if ("Accessibility" %in% names(bundle) || "Insertions" %in% names(bundle)) {
      if ("Accessibility" %in% names(bundle)) {
        bundle$Accessibility <- .filter_coverage_granges_list(
          bundle$Accessibility,
          samples = samples,
          region = region
        )
      }
      if ("Insertions" %in% names(bundle)) {
        bundle$Insertions <- .filter_coverage_granges_list(
          bundle$Insertions,
          samples = samples,
          region = region
        )
      }
      return(bundle)
    }
    return(.filter_coverage_granges_list(bundle, samples = samples, region = region))
  }

  .warn_mixed_coverage_layout(outDir, cellPop)

  list(
    Accessibility = .loadCoverageTracks(
      outDir,
      cellPop,
      trackType = "accessibility",
      samples = samples,
      region = region,
      verbose = verbose
    ),
    Insertions = .loadCoverageTracks(
      outDir,
      cellPop,
      trackType = "insertions",
      samples = samples,
      region = region,
      verbose = verbose
    )
  )
}

#' @title Migrate legacy coverage RDS bundles to structured bigWig tracks
#'
#' @description Converts `{cellPop}_CoverageFiles.RDS` bundles written by older
#'   MOCHA versions into the structured `tracks/` layout. Legacy files are kept
#'   unless `removeLegacy = TRUE`.
#'
#' @param MOCHAObj A MOCHA tileResults or Sample-Tile Matrix object with
#'   `metadata$Directory` set.
#' @param cellPopulations Character vector of cell populations to migrate, or
#'   `"ALL"` for all assays.
#' @param removeLegacy If TRUE, delete legacy RDS bundles after successful
#'   migration. Default is FALSE.
#' @param force If TRUE, overwrite existing structured track files. Default is
#'   FALSE.
#' @param verbose Display additional messages. Default is FALSE.
#'
#' @return The input MOCHA object with `metadata$CoverageLayout` set to
#'   `"tracks"`.
#' @export
#' @keywords exporting
migrateCoverageToTracks <- function(MOCHAObj,
                                    cellPopulations = "ALL",
                                    removeLegacy = FALSE,
                                    force = FALSE,
                                    verbose = FALSE) {
  outDir <- .validateCoverageDirectory(MOCHAObj, objectName = "MOCHAObj")
  cellNames <- getCellTypes(MOCHAObj)

  if (all(toupper(cellPopulations) == "ALL")) {
    cellPopulations <- cellNames
  }
  if (!all(cellPopulations %in% cellNames)) {
    stop("Some or all cell populations provided are not found.")
  }

  migrated <- character(0)
  for (cellPop in cellPopulations) {
    legacyFile <- .coverage_legacy_path(outDir, cellPop)
    layout <- .detect_coverage_layout(outDir, cellPop)

    expectedSamples <- if (file.exists(legacyFile)) {
      legacyBundle <- readRDS(legacyFile)
      names(.selectAccessibilityCoverage(legacyBundle))
    } else {
      .expected_coverage_sample_ids(outDir, cellPop)
    }

    isComplete <- .structured_tracks_complete(outDir, cellPop, expectedSamples)

    if (layout == "structured" && !force && isComplete) {
      if (verbose) {
        message("Structured tracks already complete for ", cellPop, ", skipping.")
      }
      migrated <- c(migrated, cellPop)
      next
    }

    if (!file.exists(legacyFile)) {
      if (layout == "structured" && !isComplete) {
        warning(
          "Structured coverage tracks for '", cellPop,
          "' are incomplete and no legacy RDS bundle exists at ",
          legacyFile, ". Re-run callOpenTiles(force = TRUE) to regenerate.",
          call. = FALSE
        )
      } else if (layout != "structured") {
        warning(
          "No legacy coverage bundle found for cell population '",
          cellPop,
          "' at ",
          legacyFile,
          ".",
          call. = FALSE
        )
      }
      next
    }

    if (verbose) {
      message("Migrating coverage tracks for ", cellPop)
    }
    bundle <- readRDS(legacyFile)
    writeForce <- force || (layout == "structured" && !isComplete)
    .writeCoverageTracks(
      covFiles = bundle,
      outDir = outDir,
      cellPop = cellPop,
      force = writeForce,
      verbose = verbose
    )
    migrated <- c(migrated, cellPop)

    if (removeLegacy) {
      file.remove(legacyFile)
    }
  }

  if (length(migrated) > 0) {
    MOCHAObj@metadata$CoverageLayout <- "tracks"
  }

  MOCHAObj
}
