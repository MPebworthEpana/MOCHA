#' @import MultiAssayExperiment
#' @import SummarizedExperiment
#'
#' @name MochaTileResults-class
#' @aliases MochaTileResults MochaTileResults-class
#' @exportClass MochaTileResults
#'
#' @description
#' Subclass of \code{\link[MultiAssayExperiment]{MultiAssayExperiment}} returned
#' optionally from \code{\link{callOpenTiles}} when \code{returnClass = "mocha"}.
NULL

#' @name MochaSampleTileMatrix-class
#' @aliases MochaSampleTileMatrix MochaSampleTileMatrix-class
#' @exportClass MochaSampleTileMatrix
#'
#' @description
#' Subclass of \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}}
#' returned optionally from \code{\link{getSampleTileMatrix}} when
#' \code{returnClass = "mocha"}.
NULL

.mocha_required_metadata <- function() {
  c("summarizedData", "Genome", "TxDb", "OrgDb", "Directory", "History")
}

.mocha_metadata_validity <- function(object) {
  meta <- S4Vectors::metadata(object)
  missing <- setdiff(.mocha_required_metadata(), names(meta))
  if (length(missing) > 0) {
    return(paste("Missing required metadata:", paste(missing, collapse = ", ")))
  }
  TRUE
}

setClass(
  "MochaTileResults",
  contains = "MultiAssayExperiment",
  validity = function(object) {
    .mocha_metadata_validity(object)
  }
)

setClass(
  "MochaSampleTileMatrix",
  contains = "RangedSummarizedExperiment",
  validity = function(object) {
    .mocha_metadata_validity(object)
  }
)

setAs(
  "MultiAssayExperiment",
  "MochaTileResults",
  function(from) {
    if (methods::is(from, "MochaTileResults")) {
      return(from)
    }
    methods::new("MochaTileResults", from)
  }
)

setAs(
  "RangedSummarizedExperiment",
  "MochaSampleTileMatrix",
  function(from) {
    if (methods::is(from, "MochaSampleTileMatrix")) {
      return(from)
    }
    methods::new("MochaSampleTileMatrix", from)
  }
)

#' @title Coerce to MochaTileResults
#' @param x A MultiAssayExperiment from \code{callOpenTiles}.
#' @return A \code{MochaTileResults} object.
#' @export
#' @keywords utils
asMochaTileResults <- function(x) {
  if (!methods::is(x, "MultiAssayExperiment")) {
    stop("`x` must be a MultiAssayExperiment from callOpenTiles().")
  }
  methods::as(x, "MochaTileResults")
}

#' @title Coerce to MochaSampleTileMatrix
#' @param x A RangedSummarizedExperiment from \code{getSampleTileMatrix}.
#' @return A \code{MochaSampleTileMatrix} object.
#' @export
#' @keywords utils
asMochaSTM <- function(x) {
  if (!methods::is(x, "RangedSummarizedExperiment")) {
    stop("`x` must be a RangedSummarizedExperiment from getSampleTileMatrix().")
  }
  methods::as(x, "MochaSampleTileMatrix")
}

#' @keywords internal
#' @noRd
.mocha_promote_return <- function(object, returnClass = c("legacy", "mocha")) {
  returnClass <- match.arg(returnClass)
  if (returnClass == "legacy") {
    return(object)
  }
  if (methods::is(object, "MultiAssayExperiment") &&
      !methods::is(object, "MochaTileResults")) {
    return(methods::as(object, "MochaTileResults"))
  }
  if (methods::is(object, "RangedSummarizedExperiment") &&
      !methods::is(object, "MochaSampleTileMatrix")) {
    return(methods::as(object, "MochaSampleTileMatrix"))
  }
  object
}

#' @keywords internal
#' @noRd
.mocha_repromote <- function(object) {
  if (methods::is(object, "MultiAssayExperiment")) {
    return(methods::as(object, "MochaTileResults"))
  }
  if (methods::is(object, "RangedSummarizedExperiment")) {
    return(methods::as(object, "MochaSampleTileMatrix"))
  }
  object
}

setMethod(
  "show",
  "MochaTileResults",
  function(object) {
    n_assays <- length(names(object))
    hist_len <- length(if (is.null(object@metadata$History)) list() else object@metadata$History)
    cat("MochaTileResults (MultiAssayExperiment)\n")
    cat("  cell populations:", n_assays, "\n")
    cat("  samples:", nrow(MultiAssayExperiment::colData(object)), "\n")
    cat("  directory:", object@metadata$Directory, "\n")
    cat("  history entries:", hist_len, "\n")
    callNextMethod()
  }
)

setMethod(
  "show",
  "MochaSampleTileMatrix",
  function(object) {
    n_assays <- length(SummarizedExperiment::assayNames(object))
    hist_len <- length(if (is.null(object@metadata$History)) list() else object@metadata$History)
    cat("MochaSampleTileMatrix (RangedSummarizedExperiment)\n")
    cat("  cell population assays:", n_assays, "\n")
    cat("  tiles:", nrow(object), "\n")
    cat("  samples:", ncol(object), "\n")
    cat("  directory:", object@metadata$Directory, "\n")
    cat("  history entries:", hist_len, "\n")
    callNextMethod()
  }
)
