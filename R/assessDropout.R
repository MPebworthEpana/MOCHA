#' @title Estimate a logistic dropout model for a cell population
#'
#' @description \code{estimateDropoutModel} fits a per-cell-population logistic
#'   model relating the probability of an observed zero to tile signal and sample
#'   depth, to distinguish technical dropout from biological closure.
#'
#' @param TSAM_Object A SummarizedExperiment from \code{\link{getSampleTileMatrix}}.
#' @param cellPopulation Character vector of assay (cell population) names, or
#'   \code{"all"} for every intensity assay (excluding existing
#'   \code{DropoutProb_*} assays).
#' @param additionalCovariates Optional character vector of \code{colData} column
#'   names to include in the model.
#' @param depthCols Character vector of depth covariate names. Values are read
#'   from \code{colData} when present; otherwise from
#'   \code{metadata(TSAM)$summarizedData} via fragment and cell count tables.
#'   \code{FragmentCounts} is aliased to \code{FragNumber}.
#' @param minNonzero Minimum number of non-zero observations required to include
#'   a tile in model fitting. Default is 5.
#' @param verbose Logical. Default is \code{FALSE}.
#'
#' @return An S3 object of class \code{DropoutModel} containing fitted models,
#'   formulas, depth column names, and per-cell-population diagnostics (including
#'   leave-one-sample-out AUROC with train-only tile means per fold).
#'
#' @examples
#' \dontrun{
#' model <- MOCHA::estimateDropoutModel(SampleTileMatrices, cellPopulation = "CD16 Mono")
#' }
#'
#' @export
#' @keywords downstream
estimateDropoutModel <- function(TSAM_Object,
                                 cellPopulation = "all",
                                 additionalCovariates = NULL,
                                 depthCols = c("FragNumber", "CellCounts"),
                                 minNonzero = 5,
                                 verbose = FALSE) {
  cellPops <- .resolve_dropout_cell_populations(TSAM_Object, cellPopulation)
  if (length(cellPops) == 0) {
    stop("No cell populations found for dropout modeling.")
  }

  models <- list()
  diagnostics <- list()
  formulaStr <- .dropout_formula_string(depthCols, additionalCovariates)

  for (cp in cellPops) {
    if (verbose) {
      message("Fitting dropout model for ", cp)
    }
    mat <- getCellPopMatrix(TSAM_Object, cp, NAtoZero = FALSE)
    depthDF <- .get_sample_depth_meta(TSAM_Object, cp, depthCols)
    longDF <- .build_dropout_long_data(mat, depthDF, additionalCovariates, TSAM_Object, minNonzero)
    explicitRequest <- !(length(cellPopulation) == 1L && tolower(cellPopulation) == "all")
    if (nrow(longDF) < 20) {
      if (explicitRequest) {
        stop(
          "Insufficient data to fit dropout model for ", cp,
          " (need at least 20 tile-sample rows after filtering).",
          call. = FALSE
        )
      }
      warning("Insufficient data to fit dropout model for ", cp, call. = FALSE)
      next
    }
    fit <- .fit_single_dropout_model(
      longDF, mat, depthDF, depthCols, additionalCovariates, TSAM_Object, minNonzero
    )
    models[[cp]] <- fit$model
    diagnostics[[cp]] <- fit$diagnostics
  }

  if (length(models) == 0) {
    stop("No dropout models could be fitted for the requested cell populations.")
  }

  structure(
    list(
      models = models,
      formula = formulaStr,
      depthCols = depthCols,
      additionalCovariates = additionalCovariates,
      diagnostics = diagnostics
    ),
    class = "DropoutModel"
  )
}

#' @title Assess technical vs biological zeros on a sample-tile matrix
#'
#' @description \code{assessDropout} fits a dropout model and stores predicted
#'   \code{P(zero)} per tile and sample in \code{DropoutProb_<CellPop>} assays.
#'   **Lower** values indicate a zero is more likely technical (unexpected,
#'   given tile signal and depth). Model objects are saved in
#'   \code{metadata(TSAM)$dropoutModels}.
#'
#' @param TSAM_Object A SummarizedExperiment from \code{\link{getSampleTileMatrix}}.
#' @param model Optional \code{DropoutModel} from \code{estimateDropoutModel}. If
#'   \code{NULL}, a model is estimated with default settings.
#' @param cellPopulation Cell population name(s) to assess, or \code{"all"} for
#'   every population in \code{model}. When estimating a new model, only these
#'   populations are fitted.
#' @param techThreshold Maximum \code{P(zero)} for classifying an observed zero as
#'   technical (unexpected). Default is 0.2.
#' @param bioThreshold Minimum \code{P(zero)} for classifying an observed zero as
#'   biological (expected closure). Default is 0.8.
#' @param verbose Logical. Default is \code{FALSE}.
#'
#' @return The input \code{TSAM_Object} with \code{DropoutProb_*} assays added.
#'
#' @examples
#' \dontrun{
#' TSAM_with_dropout <- MOCHA::assessDropout(SampleTileMatrices)
#' probs <- MOCHA::getDropoutProb(TSAM_with_dropout, "CD16 Mono")
#' }
#'
#' @export
#' @keywords downstream
assessDropout <- function(TSAM_Object,
                          model = NULL,
                          cellPopulation = "all",
                          techThreshold = 0.2,
                          bioThreshold = 0.8,
                          verbose = FALSE) {
  if (is.null(model)) {
    if (verbose) {
      message("No dropout model supplied; estimating with defaults.")
    }
    model <- estimateDropoutModel(
      TSAM_Object, cellPopulation = cellPopulation, verbose = verbose
    )
  }
  if (!inherits(model, "DropoutModel")) {
    stop("model must be a DropoutModel object from estimateDropoutModel().")
  }

  targetPops <- .resolve_dropout_cell_populations(TSAM_Object, cellPopulation)
  missingModel <- setdiff(targetPops, names(model$models))
  if (length(missingModel) > 0) {
    stop(
      "No dropout model available for cell population(s): ",
      paste(missingModel, collapse = ", "),
      ". Run estimateDropoutModel() for these populations or refit with ",
      "assessDropout(cellPopulation = ...).",
      call. = FALSE
    )
  }

  for (cp in targetPops) {
    refMat <- SummarizedExperiment::assay(TSAM_Object, cp)
    mat <- getCellPopMatrix(TSAM_Object, cp, dropSamples = FALSE, NAtoZero = FALSE)
    depthDF <- .get_sample_depth_meta(TSAM_Object, cp, model$depthCols)
    probMat <- .predict_dropout_matrix(mat, depthDF, model$models[[cp]], model$depthCols,
      model$additionalCovariates, TSAM_Object)
    probFull <- matrix(NA_real_, nrow = nrow(refMat), ncol = ncol(refMat), dimnames = dimnames(refMat))
    probFull[rownames(probMat), colnames(probMat)] <- probMat
    assayName <- .dropout_prob_assay_name(cp)
    SummarizedExperiment::assay(TSAM_Object, assayName) <- probFull
  }

  md <- S4Vectors::metadata(TSAM_Object)
  md$dropoutModels <- model
  md$dropoutThresholds <- list(techThreshold = techThreshold, bioThreshold = bioThreshold)
  S4Vectors::metadata(TSAM_Object) <- md

  TSAM_Object
}

#' @title Get predicted P(zero) dropout scores for a cell population
#'
#' @param TSAM_Object A SummarizedExperiment with \code{DropoutProb_*} assays from
#'   \code{\link{assessDropout}}.
#' @param cellPopulation Cell population (intensity assay) name.
#'
#' @return A matrix of predicted \code{P(zero)} for each tile and sample (lower
#'   values on observed zeros indicate more likely technical dropout).
#'
#' @export
#' @keywords utils
getDropoutProb <- function(TSAM_Object, cellPopulation) {
  assayName <- .dropout_prob_assay_name(cellPopulation)
  if (!assayName %in% SummarizedExperiment::assayNames(TSAM_Object)) {
    stop(
      "Dropout probability assay '", assayName,
      "' not found. Run assessDropout() first."
    )
  }
  SummarizedExperiment::assay(TSAM_Object, assayName)
}

#' @title Classify observed zeros as technical, biological, or ambiguous
#'
#' @param TSAM_Object A SummarizedExperiment with dropout probabilities from
#'   \code{\link{assessDropout}}.
#' @param cellPopulation Cell population name.
#' @param techThreshold See \code{\link{assessDropout}}.
#' @param bioThreshold See \code{\link{assessDropout}}.
#'
#' @return A matrix of the same dimensions as the intensity matrix with values
#'   \code{"technical"}, \code{"biological"}, \code{"ambiguous"}, or \code{NA}
#'   for non-zero observations.
#'
#' @export
#' @keywords utils
classifyZeros <- function(TSAM_Object,
                          cellPopulation,
                          techThreshold = 0.2,
                          bioThreshold = 0.8) {
  mat <- getCellPopMatrix(TSAM_Object, cellPopulation, dropSamples = FALSE, NAtoZero = FALSE)
  probs <- getDropoutProb(TSAM_Object, cellPopulation)
  probs <- probs[rownames(mat), colnames(mat), drop = FALSE]
  if (!identical(dim(mat), dim(probs))) {
    stop("Intensity and dropout probability matrices have mismatched dimensions.")
  }
  out <- matrix(NA_character_, nrow = nrow(mat), ncol = ncol(mat),
    dimnames = dimnames(mat))
  isZero <- !is.na(mat) & mat == 0
  out[isZero & probs <= techThreshold] <- "technical"
  out[isZero & probs >= bioThreshold] <- "biological"
  amb <- isZero & probs > techThreshold & probs < bioThreshold
  out[amb] <- "ambiguous"
  out
}

#' @title Plot dropout model diagnostics
#'
#' @param TSAM_Object A SummarizedExperiment with dropout results from
#'   \code{\link{assessDropout}}.
#' @param cellPopulation Cell population to plot.
#' @param type One of \code{"calibration"}, \code{"probHist"}, or \code{"perSample"}.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' \dontrun{
#' MOCHA::plotDropoutDiagnostics(TSAM_with_dropout, "CD16 Mono", type = "calibration")
#' }
#'
#' @export
#' @keywords plotting
plotDropoutDiagnostics <- function(TSAM_Object,
                                   cellPopulation,
                                   type = c("calibration", "probHist", "perSample")) {
  type <- match.arg(type)
  mat <- getCellPopMatrix(TSAM_Object, cellPopulation, dropSamples = FALSE, NAtoZero = FALSE)
  isZero <- !is.na(mat) & mat == 0

  if (type == "calibration") {
    md <- S4Vectors::metadata(TSAM_Object)
    if (is.null(md$dropoutModels) || !cellPopulation %in% names(md$dropoutModels$models)) {
      stop(
        "Calibration requires metadata(TSAM)$dropoutModels from assessDropout(). ",
        "Re-run assessDropout() before plotDropoutDiagnostics(type = 'calibration').",
        call. = FALSE
      )
    }
    depthDF <- .get_sample_depth_meta(TSAM_Object, cellPopulation, md$dropoutModels$depthCols)
    pZeroMat <- .predict_pzero_matrix(
      mat, depthDF, md$dropoutModels$models[[cellPopulation]],
      md$dropoutModels$depthCols, md$dropoutModels$additionalCovariates, TSAM_Object
    )
    depthDF_plot <- .get_sample_depth_meta(TSAM_Object, cellPopulation, c("FragNumber", "CellCounts"))
    depthCol <- if ("FragNumber" %in% colnames(depthDF_plot)) "FragNumber" else colnames(depthDF_plot)[1]
    sampleDepth <- depthDF_plot[[depthCol]]
    names(sampleDepth) <- rownames(depthDF_plot)
    depthBins <- cut(sampleDepth[colnames(mat)], breaks = 5)
    df <- data.frame(
      obs_rate = as.vector(isZero),
      pred_rate = as.vector(pZeroMat),
      depthBin = rep(as.character(depthBins), each = nrow(mat)),
      stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$pred_rate), , drop = FALSE]
    cal <- stats::aggregate(
      list(obs_rate = df$obs_rate, pred_rate = df$pred_rate),
      by = list(depthBin = df$depthBin),
      FUN = mean,
      na.rm = TRUE
    )
    return(
      ggplot2::ggplot(cal, ggplot2::aes(x = pred_rate, y = obs_rate)) +
        ggplot2::geom_point(size = 3) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        ggplot2::labs(
          title = paste("Dropout calibration by depth bin:", cellPopulation),
          x = "Mean predicted P(zero)", y = "Observed zero rate"
        ) +
        ggplot2::theme_bw()
    )
  }

  probs <- getDropoutProb(TSAM_Object, cellPopulation)
  probs <- probs[rownames(mat), colnames(mat), drop = FALSE]

  if (type == "probHist") {
    df <- data.frame(pZero = probs[isZero])
    return(
      ggplot2::ggplot(df, ggplot2::aes(x = pZero)) +
        ggplot2::geom_histogram(bins = 30, fill = "steelblue", color = "white") +
        ggplot2::labs(
          title = paste("P(zero) for observed zeros (lower = more technical):", cellPopulation),
          x = "P(zero)", y = "Count"
        ) +
        ggplot2::theme_bw()
    )
  }

  depthDF <- .get_sample_depth_meta(TSAM_Object, cellPopulation, c("FragNumber", "CellCounts"))
  depthCol <- if ("FragNumber" %in% colnames(depthDF)) "FragNumber" else colnames(depthDF)[1]
  sampleDepth <- depthDF[[depthCol]]
  names(sampleDepth) <- rownames(depthDF)

  if (type == "perSample") {
    cls <- classifyZeros(TSAM_Object, cellPopulation)
    df <- data.frame(
      Sample = colnames(mat),
      fracTechnical = vapply(seq_len(ncol(mat)), function(j) {
        z <- isZero[, j]
        if (!any(z)) {
          return(NA_real_)
        }
        mean(cls[z, j] == "technical", na.rm = TRUE)
      }, numeric(1)),
      depth = sampleDepth[colnames(mat)]
    )
    return(
      ggplot2::ggplot(df, ggplot2::aes(x = depth, y = fracTechnical, label = Sample)) +
        ggplot2::geom_point() +
        ggplot2::labs(
          title = paste("Fraction of zeros classified technical:", cellPopulation),
          x = depthCol, y = "Fraction technical"
        ) +
        ggplot2::theme_bw()
    )
  }

}

#' @title Compute group zero rates with optional dropout adjustment
#' @noRd
.compute_adjusted_zero_rates <- function(intensityMat,
                                         dropoutProb,
                                         group,
                                         adjustment = c("none", "biological_only", "weighted"),
                                         bioThreshold = 0.8) {
  adjustment <- match.arg(adjustment)
  isZero <- intensityMat == 0
  isZero[is.na(intensityMat)] <- FALSE

  group1 <- which(group == 1)
  group0 <- which(group == 0)

  if (adjustment == "none" || is.null(dropoutProb)) {
    matFilled <- intensityMat
    matFilled[is.na(matFilled)] <- 0
    isZeroLegacy <- matFilled == 0
    zero_A <- rowMeans(isZeroLegacy[, group1, drop = FALSE])
    zero_B <- rowMeans(isZeroLegacy[, group0, drop = FALSE])
    return(list(zero_A = zero_A, zero_B = zero_B))
  }

  if (!identical(dim(isZero), dim(dropoutProb))) {
    stop("dropoutProb dimensions do not match the intensity matrix.")
  }

  if (adjustment == "biological_only") {
    bioZero <- isZero & (is.na(dropoutProb) | dropoutProb >= bioThreshold)
    zero_A <- rowMeans(bioZero[, group1, drop = FALSE])
    zero_B <- rowMeans(bioZero[, group0, drop = FALSE])
  } else if (adjustment == "weighted") {
    w <- ifelse(isZero, dropoutProb, 0)
    w[is.na(w)] <- 0
    zero_A <- rowSums(w[, group1, drop = FALSE]) / length(group1)
    zero_B <- rowSums(w[, group0, drop = FALSE]) / length(group0)
  }

  list(zero_A = zero_A, zero_B = zero_B)
}

#' @noRd
.has_valid_dropout_assay <- function(TSAM_Object, cellPop) {
  assayName <- .dropout_prob_assay_name(cellPop)
  if (!assayName %in% SummarizedExperiment::assayNames(TSAM_Object)) {
    return(FALSE)
  }
  intMat <- getCellPopMatrix(TSAM_Object, cellPop, dropSamples = FALSE, NAtoZero = FALSE)
  probMat <- SummarizedExperiment::assay(TSAM_Object, assayName)
  probMat <- probMat[rownames(intMat), colnames(intMat), drop = FALSE]
  if (!identical(dim(probMat), dim(intMat))) {
    return(FALSE)
  }
  md <- S4Vectors::metadata(TSAM_Object)
  if (is.null(md$dropoutModels) || !cellPop %in% names(md$dropoutModels$models)) {
    return(FALSE)
  }
  TRUE
}

#' @title Ensure dropout probability assay exists on TSAM
#' @noRd
.ensure_dropout_assay <- function(TSAM_Object, cellPop, verbose = FALSE) {
  if (.has_valid_dropout_assay(TSAM_Object, cellPop)) {
    return(TSAM_Object)
  }
  assayName <- .dropout_prob_assay_name(cellPop)
  if (verbose) {
    message(
      "No valid ", assayName, " assay found; running assessDropout() for ", cellPop, "."
    )
  } else {
    message(
      "No dropout probability assay found for ", cellPop,
      "; running assessDropout() for ", cellPop, "."
    )
  }
  model <- estimateDropoutModel(
    TSAM_Object, cellPopulation = cellPop, verbose = verbose
  )
  TSAM_Object <- assessDropout(
    TSAM_Object, model = model, cellPopulation = cellPop, verbose = verbose
  )
  if (!.has_valid_dropout_assay(TSAM_Object, cellPop)) {
    stop(
      "Dropout-adjusted filtering requires a fitted DropoutProb assay for '",
      cellPop, "' (", assayName, "). Model fitting failed or was skipped. ",
      "Run assessDropout(SampleTileObj, cellPopulation = '", cellPop,
      "') successfully, or set dropoutAdjustment = 'none'.",
      call. = FALSE
    )
  }
  TSAM_Object
}

#' @noRd
.dropout_prob_assay_name <- function(cellPopulation) {
  paste0("DropoutProb_", gsub(" ", "_", cellPopulation))
}

#' @noRd
.resolve_dropout_cell_populations <- function(TSAM_Object, cellPopulation) {
  allAssays <- SummarizedExperiment::assayNames(TSAM_Object)
  intensityAssays <- allAssays[!grepl("^DropoutProb_", allAssays)]
  if (length(cellPopulation) == 1 && tolower(cellPopulation) == "all") {
    return(intensityAssays)
  }
  missing <- setdiff(cellPopulation, intensityAssays)
  if (length(missing) > 0) {
    stop(
      "cellPopulation(s) not found in TSAM assays: ",
      paste(missing, collapse = ", ")
    )
  }
  cellPopulation
}

#' @noRd
.get_sample_depth_meta <- function(TSAM_Object, cellPopulation, depthCols) {
  cd <- as.data.frame(SummarizedExperiment::colData(TSAM_Object))
  mat <- getCellPopMatrix(TSAM_Object, cellPopulation, dropSamples = FALSE, NAtoZero = FALSE)
  sampleIds <- colnames(mat)
  out <- data.frame(row.names = sampleIds, stringsAsFactors = FALSE)

  alias <- list(FragNumber = "FragmentCounts", CellCounts = "CellCounts")
  cdKey <- if ("Sample" %in% colnames(cd)) cd$Sample else rownames(cd)

  for (col in depthCols) {
    if (col %in% colnames(cd)) {
      vals <- cd[[col]]
      names(vals) <- cdKey
      out[[col]] <- vals[sampleIds]
    }
  }

  needFromMeta <- setdiff(depthCols, colnames(out))
  if (length(needFromMeta) > 0) {
    countTables <- tryCatch(
      .get_sample_celltype_count_tables(TSAM_Object),
      error = function(e) NULL
    )
    if (!is.null(countTables)) {
      cellRow <- cellPopulation
      if (!cellRow %in% rownames(countTables$CellCounts)) {
        cellRow <- grep(cellPopulation, rownames(countTables$CellCounts), value = TRUE)[1]
      }
      for (col in needFromMeta) {
        wideName <- if (col %in% names(alias)) alias[[col]] else col
        if (wideName %in% names(countTables)) {
          wide <- countTables[[wideName]]
          if (!is.null(cellRow) && !is.na(cellRow) && cellRow %in% rownames(wide)) {
            bioSamples <- colnames(wide)
            vals <- as.numeric(wide[cellRow, , drop = TRUE])
            names(vals) <- bioSamples
            mapIdx <- match(sampleIds, bioSamples)
            if (any(is.na(mapIdx)) && "Sample" %in% colnames(cd)) {
              mapIdx <- match(cd$Sample[match(sampleIds, cdKey)], bioSamples)
            }
            out[[col]] <- vals[mapIdx]
          }
        }
      }
    }
  }

  for (col in depthCols) {
    if (!col %in% colnames(out) || all(is.na(out[[col]]))) {
      out[[col]] <- 1
      warning(
        "Depth covariate '", col, "' not found; using constant 1 for dropout model.",
        call. = FALSE
      )
    }
    out[[col]] <- pmax(as.numeric(out[[col]]), 1e-6)
  }

  out
}

#' @noRd
.build_dropout_long_data <- function(mat,
                                     depthDF,
                                     additionalCovariates,
                                     TSAM_Object,
                                     minNonzero,
                                     muByTile = NULL) {
  if (is.null(muByTile)) {
    muByTile <- .tile_mu_from_matrix(mat)
  }

  cd <- as.data.frame(SummarizedExperiment::colData(TSAM_Object))
  rows <- list()
  for (j in seq_len(ncol(mat))) {
    sampleId <- colnames(mat)[j]
    depthRow <- if (sampleId %in% rownames(depthDF)) {
      depthDF[sampleId, , drop = FALSE]
    } else {
      depthDF[1, , drop = FALSE]
    }
    for (i in seq_len(nrow(mat))) {
      val <- mat[i, j]
      if (is.na(val)) {
        next
      }
      mu_t <- muByTile[i]
      if (is.na(mu_t)) {
        next
      }
      row <- data.frame(
        tile = rownames(mat)[i],
        sample = sampleId,
        isZero = as.integer(val == 0),
        log_mu = log(mu_t + 1),
        stringsAsFactors = FALSE
      )
      for (dc in colnames(depthRow)) {
        row[[dc]] <- as.numeric(depthRow[[dc]])
      }
      if (!is.null(additionalCovariates)) {
        sampIdx <- match(sampleId, if ("Sample" %in% colnames(cd)) cd$Sample else rownames(cd))
        for (ac in additionalCovariates) {
          row[[ac]] <- cd[[ac]][sampIdx]
        }
      }
      rows[[length(rows) + 1]] <- row
    }
  }

  if (length(rows) == 0) {
    return(data.frame())
  }
  longDF <- do.call(rbind, rows)

  keepTiles <- names(which(table(longDF$tile[longDF$isZero == 0]) >= minNonzero))
  longDF <- longDF[longDF$tile %in% keepTiles, , drop = FALSE]
  longDF
}

#' @noRd
.dropout_formula_string <- function(depthCols, additionalCovariates) {
  depthTerms <- paste0("log_", depthCols, "_sample", collapse = " + ")
  rhs <- paste(c("log_mu", depthTerms, additionalCovariates), collapse = " + ")
  paste("isZero ~", rhs)
}

#' Per-tile mean intensity among non-zero, non-NA observations.
#' @noRd
.tile_mu_from_matrix <- function(mat) {
  apply(mat, 1, function(x) {
    nz <- x[!is.na(x) & x > 0]
    if (length(nz) == 0) {
      return(NA_real_)
    }
    mean(nz)
  })
}

#' @noRd
.add_log_depth_cols <- function(df, depthCols) {
  for (dc in depthCols) {
    df[[paste0("log_", dc, "_sample")]] <- log(df[[dc]] + 1e-6)
  }
  df
}

#' @noRd
.dropout_glm_formula <- function(depthCols, additionalCovariates) {
  rhs <- c("log_mu", paste0("log_", depthCols, "_sample"), additionalCovariates)
  stats::as.formula(paste("isZero ~", paste(rhs, collapse = " + ")))
}

#' @noRd
.fit_dropout_glm <- function(longDF, depthCols, additionalCovariates) {
  longDF <- .add_log_depth_cols(longDF, depthCols)
  fml <- .dropout_glm_formula(depthCols, additionalCovariates)
  model <- tryCatch(
    stats::glm(fml, data = longDF, family = stats::binomial()),
    error = function(e) NULL
  )
  if (is.null(model)) {
    stop("glm dropout model fitting failed.")
  }
  model
}

#' @noRd
.predict_pzero <- function(model, newdata, depthCols) {
  nd <- .add_log_depth_cols(newdata, depthCols)
  stats::predict(model, newdata = nd, type = "response")
}

#' Leave-one-sample-out AUROC with train-only tile means per fold.
#' @noRd
.compute_holdout_auroc <- function(mat,
                                   depthDF,
                                   depthCols,
                                   additionalCovariates,
                                   TSAM_Object,
                                   minNonzero) {
  samples <- colnames(mat)
  if (length(samples) < 2L) {
    return(NA_real_)
  }
  cd <- as.data.frame(SummarizedExperiment::colData(TSAM_Object))
  oof_y <- integer(0)
  oof_pred <- numeric(0)

  for (s in samples) {
    train_cols <- setdiff(samples, s)
    train_mat <- mat[, train_cols, drop = FALSE]
    mu_train <- .tile_mu_from_matrix(train_mat)
    train_long <- .build_dropout_long_data(
      train_mat, depthDF, additionalCovariates, TSAM_Object, minNonzero,
      muByTile = mu_train
    )
    if (nrow(train_long) < 20L) {
      next
    }
    model <- tryCatch(
      .fit_dropout_glm(train_long, depthCols, additionalCovariates),
      error = function(e) NULL
    )
    if (is.null(model)) {
      next
    }

    j <- match(s, colnames(mat))
    sampleId <- colnames(mat)[j]
    depthRow <- if (sampleId %in% rownames(depthDF)) {
      depthDF[sampleId, , drop = FALSE]
    } else {
      depthDF[1, , drop = FALSE]
    }
    for (i in seq_len(nrow(mat))) {
      val <- mat[i, j]
      if (is.na(val)) {
        next
      }
      mu_t <- mu_train[i]
      if (is.na(mu_t)) {
        next
      }
      row <- data.frame(
        log_mu = log(mu_t + 1),
        stringsAsFactors = FALSE
      )
      for (dc in depthCols) {
        row[[dc]] <- as.numeric(depthRow[[dc]])
      }
      if (!is.null(additionalCovariates)) {
        sampIdx <- match(sampleId, if ("Sample" %in% colnames(cd)) cd$Sample else rownames(cd))
        for (ac in additionalCovariates) {
          row[[ac]] <- cd[[ac]][sampIdx]
        }
      }
      oof_pred <- c(oof_pred, .predict_pzero(model, row, depthCols))
      oof_y <- c(oof_y, as.integer(val == 0))
    }
  }

  if (length(oof_pred) == 0L || length(unique(oof_y)) < 2L) {
    return(NA_real_)
  }
  .auroc_01(oof_y, oof_pred)
}

#' @noRd
.fit_single_dropout_model <- function(longDF,
                                      mat,
                                      depthDF,
                                      depthCols,
                                      additionalCovariates,
                                      TSAM_Object,
                                      minNonzero) {
  model <- .fit_dropout_glm(longDF, depthCols, additionalCovariates)
  holdout_auc <- .compute_holdout_auroc(
    mat, depthDF, depthCols, additionalCovariates, TSAM_Object, minNonzero
  )

  list(
    model = model,
    diagnostics = list(
      aic = stats::AIC(model),
      deviance = stats::deviance(model),
      auroc = holdout_auc,
      auroc_method = "leave_one_sample_out_train_only_features",
      n = nrow(longDF)
    )
  )
}

#' @noRd
.predict_pzero_matrix <- function(mat, depthDF, model, depthCols, additionalCovariates, TSAM_Object) {
  muByTile <- .tile_mu_from_matrix(mat)
  muByTile[is.na(muByTile)] <- 1

  cd <- as.data.frame(SummarizedExperiment::colData(TSAM_Object))
  pzeroMat <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))

  for (j in seq_len(ncol(mat))) {
    sampleId <- colnames(mat)[j]
    depthRow <- if (sampleId %in% rownames(depthDF)) depthDF[sampleId, , drop = FALSE] else depthDF[1, , drop = FALSE]
    newdat <- data.frame(log_mu = log(muByTile + 1), stringsAsFactors = FALSE)
    for (dc in depthCols) {
      newdat[[dc]] <- as.numeric(depthRow[[dc]])
    }
    if (!is.null(additionalCovariates)) {
      sampIdx <- match(sampleId, if ("Sample" %in% colnames(cd)) cd$Sample else rownames(cd))
      for (ac in additionalCovariates) {
        newdat[[ac]] <- cd[[ac]][sampIdx]
      }
    }
    pzeroMat[, j] <- .predict_pzero(model, newdat, depthCols)
    pzeroMat[is.na(mat[, j]), j] <- NA_real_
  }

  pzeroMat
}

#' @noRd
.predict_dropout_matrix <- function(mat, depthDF, model, depthCols, additionalCovariates, TSAM_Object) {
  pzeroMat <- .predict_pzero_matrix(mat, depthDF, model, depthCols, additionalCovariates, TSAM_Object)
  probMat <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  isObsZero <- !is.na(mat) & mat == 0
  probMat[isObsZero] <- pzeroMat[isObsZero]
  probMat[is.na(mat)] <- NA_real_
  probMat
}

#' @noRd
.auroc_01 <- function(y, scores) {
  y <- as.integer(y)
  if (length(unique(y)) < 2) {
    return(NA_real_)
  }
  pos <- scores[y == 1]
  neg <- scores[y == 0]
  npos <- length(pos)
  nneg <- length(neg)
  if (npos == 0 || nneg == 0) {
    return(NA_real_)
  }
  ranks <- rank(scores)
  sumRanksPos <- sum(ranks[y == 1])
  (sumRanksPos - npos * (npos + 1) / 2) / (npos * nneg)
}
