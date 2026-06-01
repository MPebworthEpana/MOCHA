# Included from vignette setup chunks (not a knitr chunk file).
# Loads tutorial R code via read_chunk; skipped during knitr::purl / R CMD build tangle.

vignette_read_tutorial <- function(...) {
  tangling <- isTRUE(getOption("knitr.tangle", FALSE))
  if (requireNamespace("knitr", quietly = TRUE)) {
    tangling <- tangling || isTRUE(knitr::opts_knit$get("tangle"))
  }
  if (tangling) {
    return(invisible(NULL))
  }
  files <- unlist(list(...))
  for (chunk_file in files) {
    chunk_path <- system.file("tutorials", chunk_file, package = "MOCHA")
    if (!nzchar(chunk_path)) {
      chunk_path <- file.path("..", "inst", "tutorials", chunk_file)
    }
    if (!file.exists(chunk_path)) {
      stop("Tutorial chunk file not found: ", chunk_file)
    }
    knitr::read_chunk(chunk_path)
  }
  invisible(NULL)
}
