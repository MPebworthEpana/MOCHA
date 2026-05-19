#' @title Run differential accessibility on a set of tiles
#'
#' @description \code{estimate_differential_accessibility} internal function for computing differential accessibility on a set of tiles
#'
#'
#' @param tile_values vector of intensity values extracted from the sampleTileMatrix
#' @param group: vector of group indices (0,1)
#'
#' @return a table indicating the results for differential accessibility
#'         which includes the following pieces of information
#'         - Tile = Tile ID
#'         - P_value = P-value from two-part wilcoxon test
#'         - TestStatistic = statistic from two-part wilcoxon test
#'         - HL_EffectSize = hodges_lehmann effect size,
#'         - Case_mu = Median of nonzero values in case condition,
#'         - Case_rho = proportion of zeros in case condition ,
#'         - Control_mu = Median of nonzero values in
#'          control condition
#'         - Control_rho = proportion of zeros in control condition ,
#'
#'
#' @details The technical details of the algorithm are found in XX.
#'
#' @references XX
#'
#' @noRd
estimate_differential_accessibility <- function(tile_values,
                                                group,
                                                method = c("wilcoxon", "paired_wilcoxon", "polr"),
                                                pair_id = NULL) {
  method <- match.arg(method)
  data_vec <- as.numeric(tile_values)

  ## conduct the chosen test
  two_part_results <- switch(
    method,
    "wilcoxon" = TwoPart(data_vec, group = group, test = "wilcoxon", point.mass = 0),
    "paired_wilcoxon" = {
      if (is.null(pair_id)) {
        stop("`pair_id` is required when method = 'paired_wilcoxon'.")
      }
      TwoPartPaired(data_vec, group = group, pair_id = pair_id, test = "wilcoxon", point.mass = 0)
    },
    "polr" = .twoPart_polr(data_vec, group = group)
  )

  ## filter non-zero values for group1
  nonzero_dx <- data_vec[group == 1]
  nonzero_dx <- nonzero_dx[nonzero_dx != 0]

  ## filter non-zero values for group0
  nonzero_control <- data_vec[group == 0]
  nonzero_control <- nonzero_control[nonzero_control != 0]

  # df <- data.frame(
  #   Vals = data_vec,
  #   Group = group,
  #   Group_label = ifelse(group == 1, "case", "Control")
  # )

  ## create all pairwise combinations
  pairwise_matrix <- data.table::as.data.table(expand.grid(nonzero_dx, nonzero_control))
  pairwise_matrix$diff <- pairwise_matrix[, 1] - pairwise_matrix[, 2]

  hodges_lehmann <- stats::median(pairwise_matrix$diff)
  meanDiff <- calculateMeanDiff(data_vec, group)
  ## Create Final Results Matrix
  res <- data.frame(
    P_value = two_part_results$pvalue,
    TestStatistic = two_part_results$statistic,
    Log2FC_C = hodges_lehmann,
    MeanDiff = meanDiff,
    Case_mu = stats::median(nonzero_dx),
    Case_rho = mean(data_vec[group == 1] == 0),
    Control_mu = stats::median(nonzero_control),
    Control_rho = mean(data_vec[group == 0] == 0)
  )

  res
}

calculateMeanDiff <- function(tile_values, group) {
  a <- log2(tile_values[which(group == 1)] + 1)
  b <- log2(tile_values[which(group == 0)] + 1)

  mean_diff <- mean(a) - mean(b)
  mean_diff
}

# Proportional-odds (cumulative-logit) test for the two-group comparison.
# Bins the response into ordered tertiles (zero, low-nonzero, high-nonzero)
# and fits MASS::polr; returns the Wald chi-square on the group coefficient.
.twoPart_polr <- function(data_vec, group) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for method = 'polr'. Install MASS or pick a different method.")
  }
  if (length(unique(group)) < 2L) {
    return(list(statistic = 0, pvalue = 1))
  }
  zeros <- data_vec == 0
  non_zero_vals <- data_vec[!zeros]
  if (length(non_zero_vals) < 2L) {
    return(list(statistic = 0, pvalue = 1))
  }
  split_at <- stats::median(non_zero_vals)
  bins <- factor(
    ifelse(zeros, "zero",
           ifelse(data_vec <= split_at, "low", "high")),
    levels = c("zero", "low", "high"),
    ordered = TRUE
  )
  if (length(unique(bins)) < 2L) {
    return(list(statistic = 0, pvalue = 1))
  }
  df <- data.frame(Y = bins, G = factor(group))
  fit <- tryCatch(
    suppressWarnings(MASS::polr(Y ~ G, data = df, Hess = TRUE)),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(statistic = 0, pvalue = 1))
  }
  co <- summary(fit)$coefficients
  # First row is the group coefficient
  est <- co[1L, "Value"]
  se <- co[1L, "Std. Error"]
  if (!is.finite(se) || se <= 0) {
    return(list(statistic = 0, pvalue = 1))
  }
  z2 <- (est / se)^2
  list(statistic = z2, pvalue = 1 - stats::pchisq(z2, df = 1))
}
