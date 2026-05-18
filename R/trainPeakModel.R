#' Retrain MOCHA peak-calling models at a chosen tile size
#'
#' Fits cell-count-specific logistic regression models and Youden-optimal
#' probability thresholds following the MOCHA publication Methods (NK-cell
#' training on MACS2 pseudo-bulk labels). The returned \code{MOCHAPeakModel}
#' can be passed to \code{\link{callOpenTiles}} via \code{peakModel}.
#'
#' @param ATACFragments A \code{GRangesList} of single-cell ATAC fragments for
#'   one training cell population (typically pooled across samples).
#' @param cellColData Cell metadata; used to compute the training median
#'   fragments per cell when column \code{nFrags} is present.
#' @param blackList \code{GRanges} of regions to exclude when tiling.
#' @param groundTruthPeaks \code{GRanges} of accessible regions for training
#'   labels. Required unless \code{callPeaks = TRUE}.
#' @param tileSize Width of genomic tiles in base pairs. Default \code{500L}.
#' @param cellCol Column in fragment \code{mcols} with cell barcodes. Default
#'   \code{"RG"}.
#' @param sampleColumn Unused for training labels; reserved for compatibility.
#' @param cellSubsetSizes Integer vector of cell counts at which to fit models.
#' @param replicatesFn Function of cell count \code{n} returning the number of
#'   training replicates. Default: 10 below 50k cells, 5 otherwise.
#' @param threshMethod Method for optimal probability cutoffs per cell count:
#'   \code{"youden"} (requires \pkg{cutpointr}), \code{"f1"}, or \code{"topleft"}.
#' @param smoothMethod Coefficient smoothing: \code{"auto"} (LOESS + linear, as
#'   in \code{make_prediction}), \code{"loess"}, or \code{"linear"}.
#' @param callPeaks If \code{TRUE}, run MACS2 on pooled pseudobulk fragments to
#'   derive \code{groundTruthPeaks}. Requires \code{macs2} on \code{PATH} and
#'   \pkg{rtracklayer}.
#' @param macs2Args Named list of MACS2 arguments passed to
#'   \code{macs2 callpeak} (defaults mirror the MOCHA paper / ArchR settings).
#' @param trimPeakBp Base pairs trimmed from each end of broad peaks before
#'   overlay on tiles. Default \code{75L}.
#' @param numCores Cores for parallel replicate fitting. Default \code{1L}.
#' @param seed Random seed for cell subsampling.
#' @param verbose Print progress messages.
#'
#' @return A \code{MOCHAPeakModel} object for use with \code{peakModel} in
#'   \code{\link{callOpenTiles}}.
#'
#' @details Optional dependencies:
#' \itemize{
#'   \item \pkg{cutpointr} when \code{threshMethod = "youden"}
#'   \item \pkg{rtracklayer} and a \code{macs2} binary when \code{callPeaks = TRUE}
#' }
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("cutpointr", quietly = TRUE)) {
#'   model <- MOCHA::trainPeakModel(
#'     ATACFragments = MOCHA::exampleFragments[[1]],
#'     cellColData = MOCHA::exampleCellColData,
#'     blackList = MOCHA::exampleBlackList,
#'     groundTruthPeaks = myPeaks,
#'     tileSize = 250L,
#'     cellSubsetSizes = c(50, 100, 500),
#'     replicatesFn = function(n) 2L,
#'     numCores = 2
#'   )
#'   tiles <- MOCHA::callOpenTiles(
#'     ATACFragments = MOCHA::exampleFragments,
#'     cellColData = MOCHA::exampleCellColData,
#'     blackList = MOCHA::exampleBlackList,
#'     genome = "hg19",
#'     TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene",
#'     OrgDb = "org.Hs.eg.db",
#'     outDir = tempdir(),
#'     cellPopLabel = "Clusters",
#'     cellPopulations = "C2",
#'     peakModel = model,
#'     numCores = 1
#'   )
#' }
#' }
#'
#' @export
trainPeakModel <- function(ATACFragments,
                           cellColData,
                           blackList,
                           groundTruthPeaks = NULL,
                           tileSize = .MOCHA_DEFAULT_TILE_SIZE,
                           cellCol = "RG",
                           sampleColumn = "Sample",
                           cellSubsetSizes = c(
                             5, 10, 25, 50, 100, 250, 500, 1000,
                             2500, 5000, 10000, 25000, 50000,
                             100000, 170000
                           ),
                           replicatesFn = function(n) {
                             if (n < 50000) 10L else 5L
                           },
                           threshMethod = c("youden", "f1", "topleft"),
                           smoothMethod = c("auto", "loess", "linear"),
                           callPeaks = FALSE,
                           macs2Args = list(
                             g = "hs",
                             f = "BED",
                             nolambda = TRUE,
                             shift = -75,
                             extsize = 150,
                             broad = TRUE
                           ),
                           trimPeakBp = 75L,
                           numCores = 1L,
                           seed = 1L,
                           verbose = FALSE) {
  threshMethod <- match.arg(threshMethod)
  smoothMethod <- match.arg(smoothMethod)
  sampleColumn # reserved

  if (threshMethod == "youden" &&
      !requireNamespace("cutpointr", quietly = TRUE)) {
    stop(
      "Package 'cutpointr' is required when threshMethod = 'youden'. ",
      "Install cutpointr or use threshMethod = 'f1' or 'topleft'."
    )
  }
  if (callPeaks) {
    if (Sys.which("macs2") == "") {
      stop(
        "MACS2 binary not found on PATH. Install MACS2, pass groundTruthPeaks, ",
        "or set callPeaks = FALSE."
      )
    }
    if (!requireNamespace("rtracklayer", quietly = TRUE)) {
      stop(
        "Package 'rtracklayer' is required when callPeaks = TRUE ",
        "(BED export for MACS2)."
      )
    }
  }
  if (is.null(groundTruthPeaks) && !callPeaks) {
    stop(
      "Either supply groundTruthPeaks or set callPeaks = TRUE with macs2 on PATH."
    )
  }

  if (methods::is(ATACFragments, "GRanges")) {
    ATACFragments <- GenomicRanges::GRangesList(ATACFragments)
  } else if (!inherits(ATACFragments, "GRangesList") && !is.list(ATACFragments)) {
    stop("ATACFragments must be a GRanges, GRangesList, or list of GRanges.")
  } else {
    ATACFragments <- GenomicRanges::GRangesList(ATACFragments)
  }
  if (!methods::is(blackList, "GRanges")) {
    stop("blackList must be a GRanges object.")
  }

  pooledFrags <- if (length(ATACFragments) == 1L) {
    ATACFragments[[1]]
  } else {
    plyranges::bind_ranges(ATACFragments)
  }
  if (length(pooledFrags) == 0) {
    stop("ATACFragments contain no fragments.")
  }
  if (!(cellCol %in% colnames(GenomicRanges::mcols(pooledFrags)))) {
    stop("Fragments must contain column '", cellCol, "' with cell barcodes.")
  }

  trainingMedian <- .computeTrainingMedian(ATACFragments, cellColData, cellCol)

  if (callPeaks) {
    if (verbose) {
      message("Running MACS2 on pooled pseudobulk fragments...")
    }
    groundTruthPeaks <- .runMACS2Wrapper(
      ATACFragments = ATACFragments,
      macs2Args = macs2Args,
      verbose = verbose
    )
  }
  if (!methods::is(groundTruthPeaks, "GRanges")) {
    stop("groundTruthPeaks must be a GRanges object.")
  }

  if (verbose) {
    message("Building ", tileSize, " bp candidate tiles...")
  }
  candidateBins <- determine_dynamic_range(
    AllFragmentsList = pooledFrags,
    blackList = blackList,
    binSize = tileSize,
    doBin = FALSE
  )
  if (length(candidateBins) == 0) {
    stop("No candidate tiles after blacklist filtering.")
  }

  if (verbose) {
    message("Labeling tiles from ground-truth peaks...")
  }
  labelled <- .labelTilesFromPeaks(
    bins = candidateBins,
    peaks = groundTruthPeaks,
    fragMat = pooledFrags,
    trimPeakBp = trimPeakBp
  )

  nCellsAvail <- length(unique(as.character(GenomicRanges::mcols(pooledFrags)[[cellCol]])))
  cellSubsetSizes <- as.integer(cellSubsetSizes)
  cellSubsetSizes <- cellSubsetSizes[cellSubsetSizes >= 5L]
  cellSubsetSizes <- unique(cellSubsetSizes[cellSubsetSizes <= nCellsAvail])
  if (length(cellSubsetSizes) == 0) {
    stop(
      "No valid cellSubsetSizes (need at least 5 cells and <= ",
      nCellsAvail, " available)."
    )
  }

  if (verbose) {
    message(
      "Fitting LRMs at cell counts: ",
      paste(cellSubsetSizes, collapse = ", ")
    )
  }

  fitResults <- .fitLRMsPerCellCount(
    pooledFrags = pooledFrags,
    labelled = labelled,
    cellSubsetSizes = cellSubsetSizes,
    replicatesFn = replicatesFn,
    cellCol = cellCol,
    seed = seed,
    numCores = numCores,
    verbose = verbose
  )

  if (verbose) {
    message("Estimating optimal thresholds (", threshMethod, ")...")
  }
  threshDF <- .estimateThresholds(fitResults, method = threshMethod)

  if (verbose) {
    message("Smoothing coefficients and thresholds...")
  }
  finalModelObject <- .smoothCoefficients(fitResults$coefDF, method = smoothMethod)
  youden_threshold <- .smoothThresholds(threshDF, method = smoothMethod)

  .newPeakModel(
    finalModelObject = finalModelObject,
    youden_threshold = youden_threshold,
    tileSize = tileSize,
    trainingMedian = trainingMedian,
    method = list(thresh = threshMethod, smooth = smoothMethod),
    training = list(
      cellSubsetSizes = cellSubsetSizes,
      nCellsAvail = nCellsAvail,
      nTiles = nrow(labelled),
      nAccessible = sum(labelled$accessible),
      source = "trainPeakModel"
    )
  )
}

#' @keywords internal
#' @noRd
.computeTrainingMedian <- function(ATACFragments, cellColData, cellCol) {
  if (!is.null(cellColData) && "nFrags" %in% colnames(cellColData)) {
    return(stats::median(cellColData$nFrags, na.rm = TRUE))
  }
  pooled <- if (length(ATACFragments) == 1L) {
    ATACFragments[[1]]
  } else {
    plyranges::bind_ranges(ATACFragments)
  }
  cells <- unique(as.character(GenomicRanges::mcols(pooled)[[cellCol]]))
  tab <- table(as.character(GenomicRanges::mcols(pooled)[[cellCol]]))
  stats::median(as.integer(tab))
}

#' @keywords internal
#' @noRd
.runMACS2Wrapper <- function(ATACFragments, macs2Args, tmpDir = NULL, verbose = FALSE) {
  if (is.null(tmpDir)) {
    tmpDir <- tempfile("mocha_macs2_")
  }
  dir.create(tmpDir, recursive = TRUE, showWarnings = FALSE)

  pooled <- if (length(ATACFragments) == 1L) {
    ATACFragments[[1]]
  } else {
    plyranges::bind_ranges(ATACFragments)
  }
  bedDf <- as.data.frame(pooled)[, c("seqnames", "start", "end")]
  bedDf$start <- bedDf$start - 1L
  bedFile <- file.path(tmpDir, "pseudobulk.bed")
  utils::write.table(
    bedDf,
    bedFile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  outPrefix <- file.path(tmpDir, "mocha_macs2")
  args <- c(
    "callpeak",
    "-t", bedFile,
    "-f", macs2Args$f %||% "BED",
    "-g", macs2Args$g %||% "hs",
    "-n", outPrefix,
    "--outdir", tmpDir
  )
  if (isTRUE(macs2Args$nolambda)) {
    args <- c(args, "--nolambda")
  }
  if (!is.null(macs2Args$shift)) {
    args <- c(args, "--shift", as.character(macs2Args$shift))
  }
  if (!is.null(macs2Args$extsize)) {
    args <- c(args, "--extsize", as.character(macs2Args$extsize))
  }
  if (isTRUE(macs2Args$broad)) {
    args <- c(args, "--broad")
  }

  if (verbose) {
    message("MACS2: ", paste(args, collapse = " "))
  }
  macs2Out <- system2("macs2", args, stdout = TRUE, stderr = TRUE)
  exitCode <- attr(macs2Out, "status")
  if (is.null(exitCode)) {
    exitCode <- 0L
  }
  if (exitCode != 0L) {
    stop("MACS2 callpeak failed:\n", paste(macs2Out, collapse = "\n"))
  }

  peakFile <- paste0(outPrefix, "_peaks.broadPeak")
  if (!file.exists(peakFile)) {
    peakFile <- paste0(outPrefix, "_peaks.narrowPeak")
  }
  if (!file.exists(peakFile)) {
    stop(
      "MACS2 did not produce expected peak files under ", tmpDir, ". ",
      "Check MACS2 output."
    )
  }

  peaks <- rtracklayer::import(peakFile, format = "BED")
  if (length(peaks) == 0) {
    stop("MACS2 returned no peaks.")
  }
  peaks
}

#' @keywords internal
#' @noRd
.labelTilesFromPeaks <- function(bins, peaks, fragMat, trimPeakBp) {
  trimPeakBp <- as.integer(trimPeakBp)
  w <- IRanges::width(peaks)
  if (any(w > 2L * trimPeakBp, na.rm = TRUE)) {
    peaks <- GenomicRanges::trim(
      IRanges::narrow(peaks, start = trimPeakBp + 1L, end = -trimPeakBp)
    )
    peaks <- peaks[IRanges::width(peaks) > 0L]
  }

  tileID <- paste0(
    as.character(GenomicRanges::seqnames(bins)), ":",
    GenomicRanges::start(bins), "-", GenomicRanges::end(bins)
  )

  peakOl <- GenomicRanges::findOverlaps(bins, peaks, minoverlap = 1L)
  accessibleIDs <- unique(tileID[S4Vectors::queryHits(peakOl)])

  fragOl <- GenomicRanges::findOverlaps(bins, fragMat, minoverlap = 1L)
  fragIDs <- unique(tileID[S4Vectors::queryHits(fragOl)])

  if (length(fragIDs) == 0) {
    stop("No tiles overlap training fragments.")
  }

  keep <- tileID %in% fragIDs
  data.frame(
    seqnames = as.character(GenomicRanges::seqnames(bins))[keep],
    start = GenomicRanges::start(bins)[keep],
    end = GenomicRanges::end(bins)[keep],
    strand = as.character(GenomicRanges::strand(bins))[keep],
    accessible = tileID[keep] %in% accessibleIDs,
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
#' @noRd
.subsampleAndFit <- function(pooledFrags,
                             labelled,
                             nCells,
                             repIdx,
                             cellCol,
                             seed) {
  set.seed(seed + repIdx * 1000L + nCells)
  allCells <- unique(as.character(GenomicRanges::mcols(pooledFrags)[[cellCol]]))
  nSample <- min(as.integer(nCells), length(allCells))
  sampled <- sample(allCells, nSample, replace = FALSE)

  subFrags <- pooledFrags[
    as.character(GenomicRanges::mcols(pooledFrags)[[cellCol]]) %in% sampled
  ]
  if (length(subFrags) == 0) {
    return(NULL)
  }

  bins <- GenomicRanges::GRanges(
    seqnames = labelled$seqnames,
    ranges = IRanges::IRanges(start = labelled$start, end = labelled$end),
    strand = labelled$strand
  )

  totalFrags <- as.integer(length(subFrags))
  counts <- calculate_intensities(
    fragMat = subFrags,
    candidatePeaks = bins,
    totalFrags = totalFrags,
    cellCol = cellCol,
    verbose = FALSE
  )
  labelKey <- paste0(
    labelled$seqnames, ":", labelled$start, "-", labelled$end
  )
  counts$accessible <- labelled$accessible[match(counts$tileID, labelKey)]

  trainDF <- counts[counts$TotalIntensity > 0, ]
  if (nrow(trainDF) < 20L || length(unique(trainDF$accessible)) < 2L) {
    return(NULL)
  }

  fit <- stats::glm(
    accessible ~ TotalIntensity + maxIntensity,
    data = trainDF,
    family = stats::binomial("logit")
  )
  coefs <- stats::coef(fit)
  pred <- stats::predict(fit, type = "response")

  list(
    NumCells = nSample,
    Intercept = unname(coefs["(Intercept)"]),
    Total = unname(coefs["TotalIntensity"]),
    Max = unname(coefs["maxIntensity"]),
    pred = pred,
    obs = trainDF$accessible,
    rep = repIdx
  )
}

#' @keywords internal
#' @noRd
.fitLRMsPerCellCount <- function(pooledFrags,
                                 labelled,
                                 cellSubsetSizes,
                                 replicatesFn,
                                 cellCol,
                                 seed,
                                 numCores,
                                 verbose) {
  jobs <- lapply(cellSubsetSizes, function(n) {
    nrep <- as.integer(replicatesFn(n))
    seq_len(max(1L, nrep))
  })
  names(jobs) <- as.character(cellSubsetSizes)

  allRuns <- list()
  for (nChr in names(jobs)) {
    n <- as.integer(nChr)
    for (repIdx in jobs[[nChr]]) {
      allRuns <- append(allRuns, list(list(n = n, rep = repIdx)))
    }
  }

  runOne <- function(job) {
    .subsampleAndFit(
      pooledFrags = pooledFrags,
      labelled = labelled,
      nCells = job$n,
      repIdx = job$rep,
      cellCol = cellCol,
      seed = seed
    )
  }

  if (numCores > 1L && length(allRuns) > 1L) {
    cl <- parallel::makeCluster(numCores)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    parallel::clusterEvalQ(cl, library(MOCHA))
    ns <- asNamespace("MOCHA")
    parallel::clusterExport(
      cl,
      varlist = c(".subsampleAndFit", "calculate_intensities"),
      envir = ns
    )
    parallel::clusterExport(
      cl,
      varlist = c("pooledFrags", "labelled", "cellCol", "seed", "allRuns"),
      envir = environment()
    )
    raw <- parallel::parLapply(cl, allRuns, function(job) {
      .subsampleAndFit(
        pooledFrags = pooledFrags,
        labelled = labelled,
        nCells = job$n,
        repIdx = job$rep,
        cellCol = cellCol,
        seed = seed
      )
    })
  } else {
    raw <- lapply(allRuns, runOne)
  }

  raw <- raw[!vapply(raw, is.null, logical(1))]
  if (length(raw) == 0) {
    stop("All training replicates failed. Check fragment counts and labels.")
  }

  coefList <- lapply(raw, function(r) {
    data.frame(
      NumCells = r$NumCells,
      Intercept = r$Intercept,
      Total = r$Total,
      Max = r$Max,
      stringsAsFactors = FALSE
    )
  })
  coefDF <- do.call(rbind, coefList)
  coefDF <- stats::aggregate(
    cbind(Intercept, Total, Max) ~ NumCells,
    data = coefDF,
    FUN = stats::median
  )

  list(coefDF = coefDF, fitRuns = raw)
}

#' @keywords internal
#' @noRd
.estimateThresholds <- function(fitResults, method) {
  runs <- fitResults$fitRuns
  threshRows <- lapply(runs, function(r) {
    th <- .optimalCutpoint(r$pred, r$obs, method = method)
    data.frame(
      Ncells = r$NumCells,
      OptimalCutpoint = th,
      rep = r$rep,
      stringsAsFactors = FALSE
    )
  })
  threshAll <- do.call(rbind, threshRows)
  stats::aggregate(
    OptimalCutpoint ~ Ncells,
    data = threshAll,
    FUN = stats::median
  )
}

#' @keywords internal
#' @noRd
.optimalCutpoint <- function(pred, truth, method) {
  if (length(unique(truth)) < 2L) {
    return(0.5)
  }
  if (method == "youden" && requireNamespace("cutpointr", quietly = TRUE)) {
    df <- data.frame(pred = pred, class = truth)
    posClass <- if (is.logical(truth)) {
      TRUE
    } else {
      levels(truth)[2]
    }
    cp <- tryCatch(
      cutpointr::cutpointr(
        data = df,
        x = "pred",
        class = "class",
        pos_class = posClass,
        direction = ">=",
        method = cutpointr::maximize_metric,
        metric = cutpointr::youden
      ),
      error = function(e) NULL
    )
    if (!is.null(cp) && "optimal_cutpoint" %in% names(cp)) {
      oc <- as.numeric(cp$optimal_cutpoint[1])
      if (is.finite(oc)) {
        return(oc)
      }
    }
  }

  grid <- sort(unique(stats::na.omit(pred)))
  if (length(grid) < 2L) {
    return(0.5)
  }
  best <- 0.5
  bestScore <- -Inf
  pos <- truth
  nPos <- sum(pos)
  nNeg <- sum(!pos)
  for (th in grid) {
    tp <- sum(pred >= th & pos)
    fp <- sum(pred >= th & !pos)
    fn <- sum(pred < th & pos)
    tn <- sum(pred < th & !pos)
    if (method == "youden") {
      sens <- if (nPos > 0) tp / nPos else 0
      spec <- if (nNeg > 0) tn / nNeg else 0
      score <- sens + spec - 1
    } else if (method == "f1") {
      prec <- if (tp + fp > 0) tp / (tp + fp) else 0
      rec <- if (nPos > 0) tp / nPos else 0
      score <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
    } else {
      fpr <- if (nNeg > 0) fp / nNeg else 0
      tpr <- if (nPos > 0) tp / nPos else 0
      score <- -(fpr^2 + (1 - tpr)^2)
    }
    if (score > bestScore) {
      bestScore <- score
      best <- th
    }
  }
  best
}

#' @keywords internal
#' @noRd
.smoothCoefficients <- function(coefDF, method = "auto") {
  smoothVar <- function(varName) {
    df <- data.frame(
      NumCells = coefDF$NumCells,
      value = coefDF[[varName]]
    )
    loessFit <- stats::loess(value ~ NumCells, data = df, span = 0.75)
    linearFit <- stats::lm(value ~ NumCells, data = df)
    if (method == "loess") {
      return(list(loess = loessFit, linear = loessFit))
    }
    if (method == "linear") {
      return(list(loess = linearFit, linear = linearFit))
    }
    list(loess = loessFit, linear = linearFit)
  }

  i <- smoothVar("Intercept")
  t <- smoothVar("Total")
  m <- smoothVar("Max")

  list(
    Loess = list(Intercept = i$loess, Total = t$loess, Max = m$loess),
    Linear = list(Intercept = i$linear, Total = t$linear, Max = m$linear)
  )
}

#' @keywords internal
#' @noRd
.smoothThresholds <- function(threshDF, method = "auto") {
  df <- data.frame(
    Ncells = threshDF$Ncells,
    OptimalCutpoint = threshDF$OptimalCutpoint
  )
  df <- df[is.finite(df$OptimalCutpoint) & is.finite(df$Ncells), , drop = FALSE]
  if (nrow(df) < 2L) {
    stop("Too few finite threshold estimates to smooth.")
  }
  if (nrow(df) < 4L || method == "linear") {
    return(stats::lm(OptimalCutpoint ~ Ncells, data = df))
  }
  if (method == "loess") {
    return(stats::loess(OptimalCutpoint ~ Ncells, data = df, span = 0.75))
  }
  stats::loess(OptimalCutpoint ~ Ncells, data = df, span = 0.75)
}

#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x
