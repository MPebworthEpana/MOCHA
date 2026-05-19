make_synthetic_dropout_tsam <- function(
    n_bio_closed = 20L,
    n_tech_dropout = 20L,
    n_per_group = 5L) {
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

  counts <- matrix(0, nrow = n_tiles, ncol = n_samples, dimnames = list(tile_names, sample_names))

  # Biologically closed in group A, open in group B
  for (i in seq_len(n_bio_closed)) {
    counts[i, group_b] <- stats::runif(length(group_b), min = 12, max = 20)
  }

  # Open in B; zeros in A are technical (low depth), not biological
  set.seed(31)
  for (i in seq(n_bio_closed + 1L, n_tiles)) {
    counts[i, group_b] <- stats::runif(length(group_b), min = 14, max = 22)
    counts[i, group_a] <- 0
  }
  # Sparse non-zeros in group A for tech tiles (still mostly technical zeros)
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
    metadata = list(summarizedData = summarizedData)
  )
}
