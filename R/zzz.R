#' @keywords internal
#' @noRd
.onLoad <- function(libname, pkgname) {
  if (requireNamespace("Seurat", quietly = TRUE)) {
    methods::setOldClass("Seurat")
    methods::setMethod(
      "callOpenTiles",
      signature(ATACFragments = "Seurat"),
      .callOpenTiles_Seurat
    )
  }
  invisible()
}
