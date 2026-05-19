make_threshold_tile_results <- function(n_samples = 5L, n_tiles = 30L) {
  build_sample_gr <- function(peaks) {
    n <- length(peaks)
    gr <- GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges = IRanges::IRanges(
        start = seq(1, by = 100, length.out = n),
        end = seq(50, by = 100, length.out = n)
      )
    )
    S4Vectors::mcols(gr)$peak <- peaks
    gr
  }
  set.seed(101)
  # Make the first n_tiles/3 tiles "core" (peak in most samples), the rest sparse.
  core <- floor(n_tiles / 3)
  per_sample_peaks <- lapply(seq_len(n_samples), function(s) {
    p <- logical(n_tiles)
    p[seq_len(core)] <- stats::runif(core) > 0.1
    p[(core + 1L):n_tiles] <- stats::runif(n_tiles - core) > 0.7
    p
  })
  grl <- lapply(per_sample_peaks, build_sample_gr)
  names(grl) <- paste0("sample", seq_len(n_samples))
  re <- RaggedExperiment::RaggedExperiment(grl)
  colData_df <- S4Vectors::DataFrame(
    Sample = names(grl),
    row.names = names(grl)
  )
  MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(C2 = re),
    colData = colData_df
  )
}

test_that("suggestConsensusThreshold returns a data.frame with one row per cell population", {
  mae <- make_threshold_tile_results()
  rec <- MOCHA::suggestConsensusThreshold(mae, method = "kneedle")
  expect_s3_class(rec, "data.frame")
  expect_true(all(c("CellPopulation", "Reproducibility", "PeakNumber", "Method") %in% colnames(rec)))
  expect_equal(nrow(rec), 1L)
  expect_true(rec$Reproducibility >= 0 && rec$Reproducibility <= 1)

  rec2 <- MOCHA::suggestConsensusThreshold(mae, method = "second_derivative")
  expect_equal(nrow(rec2), 1L)
  expect_equal(rec2$Method, "second_derivative")
})

test_that(".kneedle_index picks the elbow of a synthetic L-curve", {
  # L-shape: stays high then drops sharply at x = 0.6
  x <- seq(0, 1, by = 0.05)
  y <- ifelse(x < 0.6, log10(1000), log10(50))
  idx <- MOCHA:::.kneedle_index(x, y)
  expect_equal(x[idx], 0.6, tolerance = 0.1)
})

test_that(".second_derivative_index finds the steep transition", {
  x <- seq(0, 1, length.out = 21L)
  y <- ifelse(x < 0.5, 3, 1)
  idx <- MOCHA:::.second_derivative_index(y)
  expect_true(abs(x[idx] - 0.5) <= 0.1)
})

test_that(".suggest_threshold_from_curve returns a 1-row data.frame for ungrouped curves", {
  df <- data.frame(
    Reproducibility = seq(0, 1, by = 0.1),
    PeakNumber = c(1000, 990, 980, 950, 900, 200, 150, 120, 100, 80, 50)
  )
  rec <- MOCHA:::.suggest_threshold_from_curve(df, method = "kneedle")
  expect_s3_class(rec, "data.frame")
  expect_equal(nrow(rec), 1L)
  expect_true(rec$Reproducibility >= 0 && rec$Reproducibility <= 1)
})

test_that(".suggest_threshold_from_curve handles grouped curves", {
  df <- rbind(
    data.frame(Reproducibility = seq(0, 1, by = 0.1),
               PeakNumber = c(1000, 990, 980, 950, 900, 200, 150, 120, 100, 80, 50),
               GroupName = "G1"),
    data.frame(Reproducibility = seq(0, 1, by = 0.1),
               PeakNumber = c(500, 490, 480, 470, 460, 450, 100, 80, 60, 40, 20),
               GroupName = "G2")
  )
  rec <- MOCHA:::.suggest_threshold_from_curve(df, method = "kneedle", groupColumn = "GroupName")
  expect_equal(nrow(rec), 2L)
  expect_setequal(rec$GroupName, c("G1", "G2"))
})

test_that("TwoPartPaired returns a non-significant p-value under the null", {
  set.seed(7)
  n <- 30
  v0 <- pmax(rnorm(n, 2, 0.5), 0)
  v1 <- pmax(rnorm(n, 2, 0.5), 0)
  res <- MOCHA:::TwoPartPaired(
    data = c(v1, v0),
    group = c(rep(1, n), rep(0, n)),
    pair_id = c(seq_len(n), seq_len(n))
  )
  expect_gt(res$pvalue, 0.05)
})

test_that("TwoPartPaired detects a paired shift in the non-zero component", {
  set.seed(11)
  n <- 25
  v0 <- pmax(rnorm(n, 2, 0.3), 0.5)
  v1 <- v0 + 1.0
  res <- MOCHA:::TwoPartPaired(
    data = c(v1, v0),
    group = c(rep(1, n), rep(0, n)),
    pair_id = c(seq_len(n), seq_len(n))
  )
  expect_lt(res$pvalue, 0.01)
})

test_that("TwoPartPaired detects a paired shift in dropout (McNemar)", {
  set.seed(13)
  n <- 20
  # Group 1 mostly has values; group 0 has many zeros at the matched indices
  v1 <- rep(2, n)
  v0 <- c(rep(0, 15), rep(2, 5))
  res <- MOCHA:::TwoPartPaired(
    data = c(v1, v0),
    group = c(rep(1, n), rep(0, n)),
    pair_id = c(seq_len(n), seq_len(n))
  )
  expect_lt(res$pvalue, 0.01)
})

test_that(".twoPart_polr returns p ≤ 1 and detects a large shift", {
  skip_if_not_installed("MASS")
  set.seed(17)
  v0 <- c(rep(0, 15), rnorm(10, 2, 0.2))
  v1 <- c(rep(0, 3), rnorm(22, 5, 0.5))
  res <- MOCHA:::.twoPart_polr(c(v1, v0), c(rep(1, 25), rep(0, 25)))
  expect_true(res$pvalue >= 0 && res$pvalue <= 1)
  expect_lt(res$pvalue, 0.05)
})

test_that("estimate_differential_accessibility supports method = 'paired_wilcoxon'", {
  set.seed(19)
  n <- 20
  v1 <- pmax(rnorm(n, 3, 0.4), 0)
  v0 <- pmax(rnorm(n, 3, 0.4), 0)
  res <- MOCHA:::estimate_differential_accessibility(
    tile_values = c(v1, v0),
    group = c(rep(1, n), rep(0, n)),
    method = "paired_wilcoxon",
    pair_id = c(seq_len(n), seq_len(n))
  )
  expect_s3_class(res, "data.frame")
  expect_true(all(c("P_value", "TestStatistic", "Log2FC_C", "MeanDiff") %in% colnames(res)))
})

test_that("getDifferentialAccessibleTiles rejects method = 'paired_wilcoxon' without pairColumn", {
  expect_error(
    MOCHA::getDifferentialAccessibleTiles(
      SampleTileObj = make_synthetic_sample_tile_matrix(),
      cellPopulation = "C2",
      groupColumn = "GroupA",
      foreground = "A",
      background = "B",
      method = "paired_wilcoxon"
    ),
    "pairColumn"
  )
})
