#' @title Plot to determine the reproducibility threshold
#'
#' @description \code{plotConsensus} Extracts the peak reproducibility and generates a
#' 					heuristic plots that can be used to determine the reproducibility threshold
#' 					used within getSampleTileMatrix.
#' @param tileObject A MultiAssayExperiment object from callOpenTiles,
#' @param cellPopulations the cell populations you want to visualize.
#' @param groupColumn Optional parameter, same as in getSampleTileMatrix, which defines whether you
#' 				  want to plot reproducibility within each
#' @param returnPlotList Instead of one plot with all celltypes/conditions, it returns a list of plots for each cell types
#' @param returnDFs Instead of a plot, returns a data.frame of the reproducibility across samples.
#' 						If set to false, then it plots the data.frame instead of returning it.
#' @param showSuggested Logical. If \code{TRUE}, overlay a dashed vertical
#'   line at the threshold recommended by
#'   \code{suggestConsensusThreshold()} for each cell population. Default
#'   \code{FALSE}.
#' @param suggestMethod Method passed to
#'   \code{suggestConsensusThreshold()} when \code{showSuggested = TRUE}.
#'   One of \code{"kneedle"} (default) or \code{"second_derivative"}.
#' @param numCores Number of cores to multithread over.
#'
#' @return A data.frame of reproducibility, or plots
#'
#'
#' @export
#' @keywords plotting

plotConsensus <- function(tileObject,
                          cellPopulations = "All",
                          groupColumn = NULL,
                          returnPlotList = FALSE,
                          returnDFs = FALSE,
                          showSuggested = FALSE,
                          suggestMethod = c("kneedle", "second_derivative"),
                          numCores = 1) {
  Reproducibility <- PeakNumber <- groups <- GroupName <- NULL
  suggestMethod <- match.arg(suggestMethod)

  if (all(tolower(cellPopulations) == "all")) {
    subTileResults <- tileObject
    cellPopulations <- names(tileObject)
  } else {
    if (all(cellPopulations %in% names(tileObject))) {
      subTileResults <- tileObject[names(tileObject) %in% cellPopulations]
    } else {
      stop(paste(
        "All of `cellPopulations` must present in tileResults.",
        "Check `names(tileResults)` for possible cell populations."
      ))
    }
  }

  sampleData <- SummarizedExperiment::colData(tileObject)


  iterList <- lapply(names(subTileResults), function(x) {
    list(subTileResults[[x]], sampleData, groupColumn, returnPlotList)
  })

  # return(iterList)
  # cl <- parallel::makeCluster(numCores)

  alldf <- pbapply::pblapply(cl = numCores, X = iterList, cellTypeDF)

  names(alldf) <- names(subTileResults)

  if (returnDFs) {
    return(alldf)
  }

  suggestionsByPop <- if (showSuggested) {
    lapply(alldf, function(df) .suggest_threshold_from_curve(df, method = suggestMethod, groupColumn = groupColumn))
  } else {
    NULL
  }

  if (returnPlotList) {
    if (!is.null(groupColumn)) {
      allPlots2 <- lapply(seq_along(alldf), function(x) {
        p <- ggplot2::ggplot(alldf[[x]], ggplot2::aes(x = Reproducibility, y = PeakNumber, group = GroupName, color = GroupName)) +
          ggplot2::geom_point() +
          ggplot2::ggtitle(names(alldf)[x]) +
          ggplot2::scale_y_continuous(trans = "log2") +
          ggplot2::ylab("Peak Number") +
          ggplot2::theme_bw()
        if (showSuggested) {
          p <- p + ggplot2::geom_vline(
            data = suggestionsByPop[[x]],
            mapping = ggplot2::aes(xintercept = Reproducibility, color = GroupName),
            linetype = "dashed"
          )
        }
        p
      })

      return(allPlots2)
    } else {
      allPlots <- lapply(seq_along(alldf), function(x) {
        p <- ggplot2::ggplot(alldf[[x]], ggplot2::aes(x = Reproducibility, y = PeakNumber)) +
          ggplot2::geom_point() +
          ggplot2::ggtitle(names(alldf)[x]) +
          ggplot2::scale_y_continuous(trans = "log2") +
          ggplot2::ylab("Peak Number") +
          ggplot2::theme_bw()
        if (showSuggested && !is.null(suggestionsByPop[[x]])) {
          p <- p + ggplot2::geom_vline(
            xintercept = suggestionsByPop[[x]]$Reproducibility,
            linetype = "dashed"
          )
        }
        p
      })
      return(allPlots)
    }
  } else {
    combinedDF <- do.call("rbind", alldf)
    combinedDF$CellPop <- gsub("\\.", "", gsub("[0-9]{1,3}", "", rownames(combinedDF)))

    if (!is.null(groupColumn)) {
      combinedDF$groups <- paste(combinedDF$CellPop, combinedDF$GroupName, sep = "_")
    } else {
      combinedDF$groups <- combinedDF$CellPop
    }

    p1 <- ggplot2::ggplot(combinedDF, ggplot2::aes(x = Reproducibility, y = PeakNumber, group = groups, color = groups)) +
      ggplot2::geom_line() +
      ggplot2::scale_y_continuous(trans = "log10") +
      ggplot2::ylab("Peak Number") +
      ggplot2::theme_bw()

    if (showSuggested) {
      sugDF <- do.call(rbind, lapply(seq_along(suggestionsByPop), function(i) {
        s <- suggestionsByPop[[i]]
        if (is.null(s)) return(NULL)
        s$CellPop <- names(alldf)[i]
        s$groups <- if (!is.null(groupColumn)) paste(s$CellPop, s$GroupName, sep = "_") else s$CellPop
        s
      }))
      if (!is.null(sugDF) && nrow(sugDF) > 0) {
        p1 <- p1 + ggplot2::geom_vline(
          data = sugDF,
          mapping = ggplot2::aes(xintercept = Reproducibility, color = groups),
          linetype = "dashed"
        )
      }
    }

    return(p1)
  }
}


#' @title Suggest a reproducibility threshold for consensus tile selection
#'
#' @description \code{suggestConsensusThreshold} examines the
#'   reproducibility-vs-peak-number curve produced internally by
#'   \code{plotConsensus()} and returns a recommended threshold per cell
#'   population (and optionally per group within population). Two automated
#'   methods are available: \code{"kneedle"} finds the point farthest from the
#'   secant joining the endpoints (the classic knee detection algorithm);
#'   \code{"second_derivative"} returns the reproducibility value at which the
#'   numerical second difference of \code{log10(PeakNumber)} is maximised in
#'   magnitude (the inflection in steepness).
#'
#' @param tileObject A MultiAssayExperiment object from \code{callOpenTiles()}.
#' @param cellPopulations Cell populations to evaluate (default \code{"all"}).
#' @param groupColumn Optional grouping column from \code{colData(tileObject)}.
#'   If supplied, one threshold is returned per (population, group).
#' @param method One of \code{"kneedle"} (default) or \code{"second_derivative"}.
#' @param numCores Cores passed through to \code{plotConsensus()}.
#'
#' @return A \code{data.frame} with columns \code{CellPopulation},
#'   \code{Reproducibility} (suggested threshold), \code{PeakNumber} (peaks
#'   retained at that threshold), \code{Method}, and \code{GroupName} when a
#'   grouping is requested.
#'
#' @export
#' @keywords plotting
suggestConsensusThreshold <- function(tileObject,
                                      cellPopulations = "all",
                                      groupColumn = NULL,
                                      method = c("kneedle", "second_derivative"),
                                      numCores = 1) {
  method <- match.arg(method)
  dfs <- plotConsensus(
    tileObject,
    cellPopulations = cellPopulations,
    groupColumn = groupColumn,
    returnDFs = TRUE,
    numCores = numCores
  )
  out <- do.call(rbind, lapply(names(dfs), function(pop) {
    rec <- .suggest_threshold_from_curve(dfs[[pop]], method = method, groupColumn = groupColumn)
    if (is.null(rec) || nrow(rec) == 0) {
      return(NULL)
    }
    rec$CellPopulation <- pop
    rec$Method <- method
    rec
  }))
  if (is.null(out)) {
    return(data.frame(
      CellPopulation = character(),
      Reproducibility = numeric(),
      PeakNumber = numeric(),
      Method = character(),
      stringsAsFactors = FALSE
    ))
  }
  rownames(out) <- NULL
  cols <- c("CellPopulation", "Reproducibility", "PeakNumber", "Method")
  if ("GroupName" %in% colnames(out)) {
    cols <- c("CellPopulation", "GroupName", "Reproducibility", "PeakNumber", "Method")
  }
  out[, cols, drop = FALSE]
}

# Suggest a single (or per-group) reproducibility threshold from one cell
# population's curve. Returns a data.frame with Reproducibility and PeakNumber
# (plus GroupName when grouped). Returns NULL when the curve is too short.
.suggest_threshold_from_curve <- function(df, method = "kneedle", groupColumn = NULL) {
  if (is.null(df) || !nrow(df)) {
    return(NULL)
  }
  pick <- function(sub) {
    if (nrow(sub) < 3L) {
      return(NULL)
    }
    x <- sub$Reproducibility
    y <- log10(pmax(sub$PeakNumber, 1))
    idx <- switch(
      method,
      "kneedle" = .kneedle_index(x, y),
      "second_derivative" = .second_derivative_index(y)
    )
    if (is.null(idx) || is.na(idx)) {
      return(NULL)
    }
    data.frame(
      Reproducibility = sub$Reproducibility[idx],
      PeakNumber = sub$PeakNumber[idx],
      stringsAsFactors = FALSE
    )
  }
  if (!is.null(groupColumn) && "GroupName" %in% colnames(df)) {
    do.call(rbind, lapply(split(df, df$GroupName), function(sub) {
      rec <- pick(sub)
      if (is.null(rec)) return(NULL)
      rec$GroupName <- unique(sub$GroupName)
      rec
    }))
  } else {
    pick(df)
  }
}

# Kneedle: normalise x and y to [0,1], find the point of max perpendicular
# distance from the line joining the first and last points. Works on
# monotone-decreasing curves (more reproducibility → fewer peaks).
.kneedle_index <- function(x, y) {
  n <- length(x)
  if (n < 3L) return(NA_integer_)
  rng_x <- range(x); rng_y <- range(y)
  if (diff(rng_x) == 0 || diff(rng_y) == 0) return(NA_integer_)
  xn <- (x - rng_x[1]) / diff(rng_x)
  yn <- (y - rng_y[1]) / diff(rng_y)
  # Secant between endpoints: y = m*x + b
  m <- (yn[n] - yn[1]) / (xn[n] - xn[1])
  b <- yn[1] - m * xn[1]
  # Perpendicular distance from point to line ax + by + c = 0 → m*x - y + b = 0
  d <- abs(m * xn - yn + b) / sqrt(m^2 + 1)
  which.max(d)
}

# Numerical second-difference index. Returns the position of the maximum |d²y|.
.second_derivative_index <- function(y) {
  n <- length(y)
  if (n < 3L) return(NA_integer_)
  d2 <- diff(y, differences = 2L)
  # d2 is aligned with y[2 : (n-1)]; offset by 1 to map back to original index.
  which.max(abs(d2)) + 1L
}


cellTypeDF <- function(list1 = NULL, peaksExperiment, sampleData, groupColumn, returnPlotList = FALSE) {
  if (!is.null(list1)) {
    peaksExperiment <- list1[[1]]
    sampleData <- list1[[2]]
    groupColumn <- list1[[3]]
    returnPlotList <- list1[[4]]
  }
  samplePeakMat <- RaggedExperiment::compactAssay(
    peaksExperiment,
    i = "peak"
  )

  # Identify any samples that had less than 5 cells, and therefore no peak calls
  emptySamples <- apply(samplePeakMat, 2, function(x) all(is.na(x) | !x))

  # Now we'll filter out those samples, and replace all NAs with zeros.
  samplePeakMat <- samplePeakMat[, !emptySamples, drop = FALSE]
  samplePeakMat[is.na(samplePeakMat)] <- FALSE
  sampleData <- sampleData[rownames(sampleData) %in% names(which(!emptySamples)), ]

  if (ncol(samplePeakMat) == 0L) {
    return(data.frame(Reproducibility = numeric(), PeakNumber = numeric()))
  }

  if (is.null(groupColumn)) {

    # The number of samples for reproducibility will be all samples with peak calls
    nSamples <- sum(!emptySamples)
    reproducibility_perc <- seq(0, 1, by = 1 / nSamples)

    # Count the number of TRUE peaks in each row and divide by nSamples
    peakReps <- rowSums(samplePeakMat) / nSamples
    RepPeaks <- sapply(reproducibility_perc, function(x) sum(peakReps >= x))
    TruePeaks <- data.frame("Reproducibility" = reproducibility_perc, "PeakNumber" = RepPeaks)
  } else {

    # Get consensus peaks for groupings of samples
    groups <- unique(sampleData[[groupColumn]])
    TruePeaks <- data.frame()
    for (group in groups) {
      # Filter sample-peak matrix to samples in this group
      samplesInGroupDF <- sampleData[sampleData[[groupColumn]] == group, ]
      samplesInGroup <- rownames(samplesInGroupDF)
      groupSamplePeakMat <- samplePeakMat[, samplesInGroup, drop = FALSE]

      reproducibility_perc <- seq(0, 1, by = 1 / length(samplesInGroup))

      peakReps_tmp <- rowSums(groupSamplePeakMat) / length(samplesInGroup)
      RepPeaks_tmp <- sapply(reproducibility_perc, function(x) sum(peakReps_tmp >= x))
      TruePeaks_tmp <- data.frame("Reproducibility" = reproducibility_perc, "PeakNumber" = RepPeaks_tmp)

      TruePeaks_tmp$GroupName <- rep(group, length(RepPeaks_tmp))

      TruePeaks <- rbind(TruePeaks, TruePeaks_tmp)
    }
  }

  return(TruePeaks)
}
