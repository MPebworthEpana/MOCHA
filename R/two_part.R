#' @title \code{TwoPart}
#'
#' @description \code{TwoPart} conducts a 2-part Wilcoxon test.
#'
#'
#' @param data is a vector containing the values for the 2-part test
#' @param group is a vector containing the two groups for comparison (in binary, 0,1)
#' @param test which type of test to use, "t.test" or "wilcoxon"
#' @param point.mass point mass
#'
#' @return ans is a list containing the test statistic & p-value corresponding to the
#'         two-part Wilcoxon test.
#'
#' @details The code from the two part test is taken from Taylor & pollard (2009) (reference below) which illustrates how to create a composite test from zero-inflated data.
#'
#' @references Taylor, Sandra, and Katherine Pollard. "Hypothesis tests for point-mass mixture data #' with application toomics data with many zero values." Statistical applications in genetics and #' molecular biology 8.1 (2009).
#'
#' @noRd

# Function for calculating two-part statistics
TwoPart <- function(data, group, test = "wilcoxon", point.mass = 0) {
  Index1 <- c(group == 1)
  Group1 <- data[Index1]
  Group0 <- data[!Index1]
  n1 <- length(Group1)
  n2 <- length(Group0)
  obs <- c(n1, n2)
  success <- c(sum(Group1 != point.mass), sum(Group0 != point.mass))
  pointmass <- obs - success

  if (sum(success) == 0) {
    T2 <- 0
    B2 <- 0
  } else if ((success[1] == 0) | (success[2] == 0)) {
    T2 <- 0
    B2 <- stats::prop.test(pointmass, obs)$statistic
  } else if ((success[1] == 1) | (success[2] == 1)) {
    T2 <- 0
    B2 <- stats::prop.test(pointmass, obs)$statistic
  } else {
    uniq1 <- length(unique(Group1[Group1 != point.mass]))
    uniq2 <- length(unique(Group0[Group0 != point.mass]))
    if ((uniq1 < 2) & (uniq2 < 2)) {
      T2 <- 0
      if (sum(pointmass) == 0) {
        B2 <- 0
      } else {
        B2 <- stats::prop.test(pointmass, obs)$statistic
      }
    } else if (sum(pointmass) == 0) {
      B2 <- 0
      if (test == "t.test") {
        T2 <- stats::t.test(data ~ group)$statistic^2
      }
      if (test == "wilcoxon") {
        W <- stats::wilcox.test(data ~ group, exact = FALSE)$statistic
        mu <- (n1 * n2) / 2
        sigma <- sqrt((n1 * n2 * (n1 + n2 + 1)) / 12)
        T2 <- ((abs(W - mu) - 0.5) / sigma)^2
      }
    } else {
      B2 <- stats::prop.test(pointmass, obs)$statistic
      contIndex <- data != point.mass
      cont <- data[contIndex]
      cGroup <- group[contIndex]
      n1c <- sum(cGroup == 1)
      n2c <- sum(cGroup == 0)

      if (test == "t.test") {
        T2 <- stats::t.test(cont ~ cGroup)$statistic^2
      }
      if (test == "wilcoxon") {
        W <- stats::wilcox.test(cont ~ cGroup, exact = FALSE)$statistic
        mu <- (n1c * n2c) / 2
        sigma <- sqrt((n1c * n2c * (n1c + n2c + 1)) / 12)
        T2 <- ((abs(W - mu) - 0.5) / sigma)^2
      }
    }
  }
  X2 <- B2 + T2
  if ((T2 == 0) | (B2 == 0)) {
    X2pv <- 1 - stats::pchisq(X2, 1)
  } else {
    X2pv <- 1 - stats::pchisq(X2, 2)
  }
  ans <- list(statistic = X2, pvalue = X2pv)
  return(ans)
}


#' @title \code{TwoPartPaired}
#'
#' @description \code{TwoPartPaired} is the paired analogue of \code{TwoPart}:
#'   it combines a McNemar statistic on the binary (zero / non-zero) component
#'   with a paired Wilcoxon signed-rank (or paired t-test) statistic on the
#'   continuous component, using only pairs where both samples are non-zero.
#'   The two squared standardized statistics are summed into a chi-square
#'   with 1 or 2 degrees of freedom depending on which components are
#'   non-degenerate, matching the construction used by \code{TwoPart}.
#'
#' @param data Numeric vector of values.
#' @param group Vector indicating the two groups for comparison (binary 0/1).
#' @param pair_id Vector of pair identifiers. Observations sharing a value
#'   are treated as paired across the two groups.
#' @param test Continuous-component test, \code{"wilcoxon"} (default) or
#'   \code{"t.test"}.
#' @param point.mass Numeric value treated as the point mass. Default 0.
#'
#' @return A list with \code{statistic} (the combined chi-square) and
#'   \code{pvalue}.
#'
#' @noRd
TwoPartPaired <- function(data, group, pair_id, test = "wilcoxon", point.mass = 0) {
  if (length(data) != length(group) || length(data) != length(pair_id)) {
    stop("`data`, `group`, and `pair_id` must be the same length.")
  }
  pair_id <- as.character(pair_id)
  g <- group == 1
  # Build paired vectors: for each pair_id with one obs per group, align.
  pairs <- intersect(
    unique(pair_id[g]),
    unique(pair_id[!g])
  )
  if (length(pairs) < 2L) {
    return(list(statistic = 0, pvalue = 1))
  }
  v1 <- data[g][match(pairs, pair_id[g])]
  v0 <- data[!g][match(pairs, pair_id[!g])]

  # Binary component: McNemar on discordant pairs (zero in one group, non-zero
  # in the other). The two off-diagonal counts give the McNemar statistic.
  z1 <- v1 == point.mass
  z0 <- v0 == point.mass
  b01 <- sum(!z1 & z0) # non-zero in group 1, zero in group 0
  b10 <- sum(z1 & !z0)
  if ((b01 + b10) > 0L) {
    B2 <- ((abs(b01 - b10) - 1)^2) / (b01 + b10) # continuity-corrected McNemar
  } else {
    B2 <- 0
  }

  # Continuous component: paired test on pairs where both are non-zero.
  both <- !z1 & !z0
  n_both <- sum(both)
  if (n_both >= 2L) {
    if (test == "t.test") {
      T2 <- tryCatch(
        stats::t.test(v1[both], v0[both], paired = TRUE)$statistic^2,
        error = function(e) 0
      )
    } else {
      w <- tryCatch(
        stats::wilcox.test(v1[both], v0[both], paired = TRUE, exact = FALSE),
        error = function(e) NULL
      )
      if (is.null(w)) {
        T2 <- 0
      } else {
        # Standardize the signed-rank statistic to a z² using its null moments.
        n <- n_both
        mu <- n * (n + 1) / 4
        sigma <- sqrt(n * (n + 1) * (2 * n + 1) / 24)
        T2 <- if (sigma > 0) ((abs(unname(w$statistic) - mu) - 0.5) / sigma)^2 else 0
      }
    }
  } else {
    T2 <- 0
  }

  X2 <- B2 + T2
  if (B2 == 0 || T2 == 0) {
    pv <- 1 - stats::pchisq(X2, df = 1)
  } else {
    pv <- 1 - stats::pchisq(X2, df = 2)
  }
  list(statistic = X2, pvalue = pv)
}
