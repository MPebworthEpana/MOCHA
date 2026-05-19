#' @title Get cell population names from a MOCHA object
#' @param x A \code{MochaTileResults} or \code{MochaSampleTileMatrix} object.
#' @return Character vector of cell population names.
#' @export
#' @keywords utils
setGeneric("cellTypes", function(x) standardGeneric("cellTypes"))

#' @title Get open tiles from a MochaTileResults object
#' @param x A \code{MochaTileResults} object.
#' @param cellPopulations Cell populations to return, or \code{"all"}.
#' @param returnType \code{"GRangesList"} or \code{"data.frame"}.
#' @return Open tiles per cell population.
#' @export
#' @keywords utils
setGeneric(
  "openTiles",
  function(x, cellPopulations = "all", returnType = c("GRangesList", "data.frame")) {
    standardGeneric("openTiles")
  }
)

setMethod(
  "cellTypes",
  signature(x = "MochaTileResults"),
  function(x) {
    MOCHA::getCellTypes(x)
  }
)

setMethod(
  "cellTypes",
  signature(x = "MochaSampleTileMatrix"),
  function(x) {
    MOCHA::getCellTypes(x)
  }
)

setMethod(
  "openTiles",
  signature(x = "MochaTileResults"),
  function(x, cellPopulations = "all", returnType = c("GRangesList", "data.frame")) {
    MOCHA::getOpenTiles(x, cellPopulations = cellPopulations, returnType = returnType)
  }
)
