#' @keywords internal
#' @noRd
.onLoad <- function(libname, pkgname) {
  if (requireNamespace("Seurat", quietly = TRUE)) {
    meth <- tryCatch(
      methods::getMethod(
        "callOpenTiles",
        signature = c(ATACFragments = "Seurat")
      ),
      error = function(e) NULL
    )
    if (is.null(meth)) {
      methods::setMethod(
        "callOpenTiles",
        signature(ATACFragments = "Seurat"),
        .callOpenTiles_Seurat
      )
    }
  }
  invisible()
}
