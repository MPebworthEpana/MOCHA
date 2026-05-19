#' Default genomic tile width (bp) for peak calling and training.
#' @keywords internal
#' @noRd
.MOCHA_DEFAULT_TILE_SIZE <- 500L

#' MOCHA peak-calling model object
#'
#' Container for logistic-regression coefficients and Youden-threshold models
#' produced by \code{\link{trainPeakModel}} or the bundled defaults.
#'
#' @param finalModelObject List with \code{Loess} and \code{Linear} sublists,
#'   each containing \code{Intercept}, \code{Total}, and \code{Max} fitted
#'   models over cell count.
#' @param youden_threshold A \code{loess} object predicting optimal probability
#'   cutoffs from cell count (\code{OptimalCutpoint ~ Ncells}).
#' @param tileSize Width of genomic tiles in base pairs used during training.
#' @param trainingMedian Median fragments per cell in the training data; used
#'   as the numerator when calibrating study signal (\code{trainingMedian / studySignal}).
#' @param method List describing threshold and smoothing methods.
#' @param training Optional list of training metadata (cell counts, etc.).
#' @keywords internal
#' @noRd
.newPeakModel <- function(finalModelObject,
                          youden_threshold,
                          tileSize,
                          trainingMedian,
                          method = list(),
                          training = list()) {
  structure(
    list(
      finalModelObject = finalModelObject,
      youden_threshold = youden_threshold,
      tileSize = as.integer(tileSize),
      trainingMedian = as.numeric(trainingMedian),
      method = method,
      training = training
    ),
    class = "MOCHAPeakModel"
  )
}

#' @keywords internal
#' @noRd
.validatePeakModel <- function(x) {
  if (!inherits(x, "MOCHAPeakModel")) {
    stop("peakModel must be a MOCHAPeakModel object from trainPeakModel().")
  }
  req <- c("finalModelObject", "youden_threshold", "tileSize", "trainingMedian")
  missing <- setdiff(req, names(x))
  if (length(missing) > 0) {
    stop("Invalid MOCHAPeakModel: missing ", paste(missing, collapse = ", "))
  }
  if (!is.numeric(x$tileSize) || length(x$tileSize) != 1L || x$tileSize < 1L) {
    stop("MOCHAPeakModel$tileSize must be a positive integer.")
  }
  if (!is.numeric(x$trainingMedian) || length(x$trainingMedian) != 1L ||
      !is.finite(x$trainingMedian) || x$trainingMedian <= 0) {
    stop("MOCHAPeakModel$trainingMedian must be a positive finite number.")
  }
  fmo <- x$finalModelObject
  for (fam in c("Loess", "Linear")) {
    if (!fam %in% names(fmo)) {
      stop("finalModelObject must contain a '", fam, "' component.")
    }
    for (coefName in c("Intercept", "Total", "Max")) {
      if (!coefName %in% names(fmo[[fam]])) {
        stop("finalModelObject$", fam, " must contain '", coefName, "'.")
      }
      fit <- fmo[[fam]][[coefName]]
      if (!inherits(fit, c("loess", "lm"))) {
        stop(
          "finalModelObject$", fam, "$", coefName,
          " must be a loess or lm object."
        )
      }
    }
  }
  if (!inherits(x$youden_threshold, c("loess", "lm"))) {
    stop(
      "youden_threshold must be a loess or lm object (OptimalCutpoint ~ Ncells)."
    )
  }
  invisible(x)
}

#' Default bundled peak-calling model (500 bp tiles, COVID19 NK training).
#' @keywords internal
#' @noRd
.defaultPeakModel <- function() {
  .newPeakModel(
    finalModelObject = MOCHA::finalModelObject,
    youden_threshold = MOCHA::youden_threshold,
    tileSize = .MOCHA_DEFAULT_TILE_SIZE,
    trainingMedian = 3668,
    method = list(thresh = "youden", smooth = "auto"),
    training = list(source = "bundled")
  )
}

#' @export
#' @method print MOCHAPeakModel
print.MOCHAPeakModel <- function(x, ...) {
  x <- .validatePeakModel(x)
  cat("MOCHA peak-calling model\n")
  cat("  tileSize:         ", x$tileSize, " bp\n", sep = "")
  cat("  trainingMedian:   ", round(x$trainingMedian, 1), " fragments/cell\n", sep = "")
  if (length(x$method) > 0) {
    threshLab <- if (is.null(x$method$thresh)) "unknown" else x$method$thresh
    smoothLab <- if (is.null(x$method$smooth)) "unknown" else x$method$smooth
    cat("  threshold method: ", threshLab, "\n", sep = "")
    cat("  smooth method:    ", smoothLab, "\n", sep = "")
  }
  if (!is.null(x$training$cellSubsetSizes)) {
    cat(
      "  trained cell counts: ",
      paste(x$training$cellSubsetSizes, collapse = ", "),
      "\n",
      sep = ""
    )
  }
  if (!is.null(x$training$source)) {
    cat("  source:           ", x$training$source, "\n", sep = "")
  }
  invisible(x)
}
