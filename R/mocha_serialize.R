#' @keywords internal
#' @noRd
.mocha_coerce_loaded_object <- function(object) {
  cls <- class(object)
  if (methods::is(object, "MochaTileResults")) {
    return(object)
  }
  if (methods::is(object, "MochaSampleTileMatrix")) {
    return(object)
  }
  if (any(cls == "MochaTileResults") && methods::is(object, "MultiAssayExperiment")) {
    return(methods::as(object, "MochaTileResults"))
  }
  if (any(cls == "MochaSampleTileMatrix") && methods::is(object, "RangedSummarizedExperiment")) {
    return(methods::as(object, "MochaSampleTileMatrix"))
  }
  if ("MochaTileResults" %in% cls && !methods::isClass("MochaTileResults")) {
    message(
      "Loaded object was saved as MochaTileResults but this MOCHA build ",
      "does not define that class; coercing to MultiAssayExperiment."
    )
    return(methods::as(object, "MultiAssayExperiment"))
  }
  if ("MochaSampleTileMatrix" %in% cls && !methods::isClass("MochaSampleTileMatrix")) {
    message(
      "Loaded object was saved as MochaSampleTileMatrix but this MOCHA build ",
      "does not define that class; coercing to RangedSummarizedExperiment."
    )
    return(methods::as(object, "RangedSummarizedExperiment"))
  }
  object
}
