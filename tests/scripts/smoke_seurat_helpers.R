#' Load blacklist for Seurat ingest smoke tests.
#'
#' Uses MOCHA exampleBlackList when the lazy-loaded package data is available
#' (after \code{load_all()} or \code{library(MOCHA)}). Otherwise falls back to a
#' small chr1 GRanges that does not overlap the mini fragment fixture
#' (peaks around 760101--760130).
#' @noRd
smoke_load_blacklist <- function() {
  e <- new.env(parent = emptyenv())
  loaded <- tryCatch(
    {
      utils::data("exampleBlackList", package = "MOCHA", envir = e)
      TRUE
    },
    error = function(err) FALSE
  )
  if (loaded && exists("exampleBlackList", envir = e, inherits = FALSE)) {
    bl <- e$exampleBlackList
    if (inherits(bl, "GRanges") && length(bl) > 0L) {
      message("Using MOCHA::exampleBlackList (", length(bl), " ranges)")
      return(bl)
    }
  }

  message("Using built-in mini blacklist for smoke test (chr1:1-750100)")
  GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1L, end = 750100L),
    strand = "*"
  )
}

#' Report missing optional dependencies for callOpenTiles smoke.
#' @return Character vector of missing package names (empty if all present).
#' @noRd
smoke_callOpenTiles_missing_pkgs <- function() {
  required <- c(
    "BSgenome.Hsapiens.UCSC.hg19",
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "org.Hs.eg.db"
  )
  required[!vapply(required, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
}
