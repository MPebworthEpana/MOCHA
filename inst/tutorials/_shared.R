# Shared helpers for MOCHA tutorial scripts (inst/tutorials/).
# Sourced by run_tutorials.R and vignette setup; not a knitr chunk file.

tutorial_repo_root <- function() {
  normalizePath(
    file.path(dirname(getwd()), ".."),
    winslash = "/",
    mustWork = FALSE
  )
}

tutorial_dir <- function() {
  p <- system.file("tutorials", package = "MOCHA")
  if (!nzchar(p)) {
    candidates <- c(
      file.path(getwd(), "inst", "tutorials"),
      file.path(dirname(sys.frame(1)$ofile %||% ""), "."),
      normalizePath(
        file.path(
          Sys.getenv("MOCHA_REPO_ROOT", unset = normalizePath(getwd(), mustWork = FALSE)),
          "inst",
          "tutorials"
        ),
        mustWork = FALSE
      )
    )
    p <- candidates[dir.exists(candidates)][1]
  }
  if (is.na(p) || !nzchar(p) || !dir.exists(p)) {
    stop("Could not locate inst/tutorials/")
  }
  normalizePath(p, winslash = "/", mustWork = TRUE)
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x)) y else x

tutorial_read_chunks <- function(...) {
  files <- unlist(list(...))
  for (f in files) {
    p <- system.file("tutorials", f, package = "MOCHA")
    if (!nzchar(p)) {
      p <- file.path(tutorial_dir(), f)
    }
    if (!file.exists(p)) {
      stop("Tutorial chunk file not found: ", f)
    }
    knitr::read_chunk(p)
  }
}

tutorial_check_deps <- function() {
  has_deps <- all(
    requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE),
    requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE),
    requireNamespace("org.Hs.eg.db", quietly = TRUE)
  )
  if (!has_deps) {
    stop(
      "Install suggested packages: BSgenome.Hsapiens.UCSC.hg19, ",
      "TxDb.Hsapiens.UCSC.hg38.knownGene, org.Hs.eg.db"
    )
  }
  invisible(TRUE)
}

tutorial_synthetic_stm <- function(n_tiles = 40L, n_per_group = 4L) {
  n_samples <- n_per_group * 2L
  tile_names <- paste0(
    "chr1:",
    seq(10000, by = 500, length.out = n_tiles),
    "-",
    seq(10499, by = 500, length.out = n_tiles)
  )
  set.seed(42)
  counts <- matrix(
    stats::runif(n_tiles * n_samples, min = 8, max = 16),
    nrow = n_tiles,
    dimnames = list(tile_names, paste0("sample_", seq_len(n_samples)))
  )
  gr <- MOCHA::StringsToGRanges(tile_names)
  S4Vectors::mcols(gr)$C2 <- TRUE
  sample_names <- colnames(counts)
  cd <- data.frame(
    Sample = sample_names,
    PassQC = TRUE,
    GroupA = rep(c("A", "B"), each = n_per_group),
    stringsAsFactors = FALSE
  )
  frag_counts <- rep(50000, length(sample_names))
  cell_counts <- rep(500, length(sample_names))
  names(frag_counts) <- sample_names
  names(cell_counts) <- sample_names
  summarizedData <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      CellCounts = matrix(cell_counts, nrow = 1, dimnames = list("C2", sample_names)),
      FragmentCounts = matrix(frag_counts, nrow = 1, dimnames = list("C2", sample_names))
    )
  )
  SummarizedExperiment::SummarizedExperiment(
    assays = list(C2 = counts),
    rowRanges = gr,
    colData = cd,
    metadata = list(
      Genome = "BSgenome.Hsapiens.UCSC.hg19",
      summarizedData = summarizedData
    )
  )
}

tutorial_synthetic_dropout_stm <- function() {
  n_bio_closed <- 20L
  n_tech_dropout <- 20L
  n_per_group <- 5L
  n_tiles <- n_bio_closed + n_tech_dropout
  n_samples <- n_per_group * 2L
  tile_names <- paste0(
    "chr1:",
    seq(10000, by = 500, length.out = n_tiles),
    "-",
    seq(10499, by = 500, length.out = n_tiles)
  )
  sample_names <- paste0("sample_", seq_len(n_samples))
  group_a <- sample_names[seq_len(n_per_group)]
  group_b <- sample_names[(n_per_group + 1L):n_samples]
  counts <- matrix(
    0,
    nrow = n_tiles,
    ncol = n_samples,
    dimnames = list(tile_names, sample_names)
  )
  for (i in seq_len(n_bio_closed)) {
    counts[i, group_b] <- stats::runif(length(group_b), min = 12, max = 20)
  }
  set.seed(31)
  for (i in seq(n_bio_closed + 1L, n_tiles)) {
    counts[i, group_b] <- stats::runif(length(group_b), min = 14, max = 22)
    counts[i, group_a] <- 0
  }
  for (i in seq(n_bio_closed + 1L, n_tiles)) {
    counts[i, group_a[1:2]] <- stats::runif(2, min = 1, max = 3)
  }
  gr <- MOCHA::StringsToGRanges(tile_names)
  S4Vectors::mcols(gr)$C2 <- TRUE
  frag_counts <- rep(50000, n_samples)
  cell_counts <- rep(500, n_samples)
  names(frag_counts) <- sample_names
  names(cell_counts) <- sample_names
  frag_counts[group_a] <- 500
  cell_counts[group_a] <- 10
  cd <- data.frame(
    Sample = sample_names,
    PassQC = TRUE,
    GroupA = rep(c("A", "B"), each = n_per_group),
    FragNumber = frag_counts,
    CellCounts = cell_counts,
    row.names = sample_names,
    stringsAsFactors = FALSE
  )
  summarizedData <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      CellCounts = matrix(cell_counts, nrow = 1, dimnames = list("C2", sample_names)),
      FragmentCounts = matrix(frag_counts, nrow = 1, dimnames = list("C2", sample_names))
    )
  )
  SummarizedExperiment::SummarizedExperiment(
    assays = list(C2 = counts),
    rowRanges = gr,
    colData = cd,
    metadata = list(
      Genome = "BSgenome.Hsapiens.UCSC.hg19",
      summarizedData = summarizedData
    )
  )
}

tutorial_try_run <- function(expr, label = "step") {
  tryCatch(
    force(expr),
    error = function(e) {
      message("Tutorial ", label, " skipped: ", conditionMessage(e))
      invisible(NULL)
    }
  )
}
