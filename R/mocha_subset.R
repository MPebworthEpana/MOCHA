#' @keywords internal
#' @noRd
.subset_mocha_object <- function(Object,
                                 subsetBy,
                                 groupList,
                                 removeNA = TRUE,
                                 subsetPeaks = TRUE,
                                 verbose = FALSE) {
  was_mocha_mae <- methods::is(Object, "MochaTileResults")
  was_mocha_stm <- methods::is(Object, "MochaSampleTileMatrix")

  summarizedData <- S4Vectors::metadata(Object)$summarizedData
  sampleData <- SummarizedExperiment::colData(Object)

  if (!(subsetBy %in% colnames(sampleData)) & !grepl("celltype", tolower(subsetBy))) {
    stop(
      "Variable given in subsetBy is not in the colData of the input Object.",
      "subsetBy must either be 'celltype', or a column name within colData(Object)."
    )
  }

  if (grepl("celltype", tolower(subsetBy))) {
    if (methods::is(Object, "MultiAssayExperiment")) {
      if ((subsetBy %in% colnames(sampleData)) & grepl("celltype", tolower(subsetBy))) {
        if (verbose) {
          warning(
            "subsetBy is set to 'celltype', but that is also a column name ",
            "within the colData of the input Object. The object will be filtered",
            " by cell type annotation, not by colData of the input Object."
          )
        }
      }

      if (!all(groupList %in% names(Object))) {
        stop("groupList includes celltypes not found within Object.")
      }

      newObject <- MultiAssayExperiment::subsetByAssay(Object, groupList)
      newObject@metadata$summarizedData <- summarizedData[groupList, ]
      if (was_mocha_mae) {
        newObject <- methods::as(newObject, "MochaTileResults")
      }
      return(newObject)
    }

    if (methods::is(Object, "RangedSummarizedExperiment")) {
      if (!all(groupList %in% names(SummarizedExperiment::assays(Object)))) {
        stop("groupList includes celltypes not found within Object.")
      }

      dropoutAssays <- vapply(groupList, .dropout_prob_assay_name, character(1))
      keepAssays <- c(groupList, dropoutAssays)
      keepIdx <- which(names(SummarizedExperiment::assays(Object)) %in% keepAssays)
      SummarizedExperiment::assays(Object) <- SummarizedExperiment::assays(Object)[keepIdx]
      Object@metadata$summarizedData <- summarizedData[groupList, ]

      if (subsetPeaks) {
        rowMeta <- GenomicRanges::mcols(SummarizedExperiment::rowRanges(Object))[, groupList, drop = FALSE]
        if (!is.null(dim(rowMeta))) {
          rowMeta <- rowSums(as.data.frame(rowMeta)) > 0
        }
        Object <- Object[rowMeta, ]
      }

      if (was_mocha_stm) {
        Object <- methods::as(Object, "MochaSampleTileMatrix")
      }
      return(Object)
    }
  }

  if (subsetBy %in% colnames(sampleData)) {
    if (!all(groupList %in% unique(sampleData[[subsetBy]]))) {
      stop(
        stringr::str_interp(
          "groupList includes names not found within the column '${subsetBy}'"
        ),
        " in the sample metadata. (see `colData(Object)`). "
      )
    }

    if (removeNA) {
      keepSamples <- rownames(sampleData)[which(sampleData[[subsetBy]] %in% groupList)]
    } else {
      keepSamples <- rownames(sampleData)[which(
        sampleData[[subsetBy]] %in% groupList | is.na(sampleData[[subsetBy]])
      )]
    }

    Object <- Object[, keepSamples]
    Object@metadata$summarizedData <- summarizedData[, keepSamples]

    if (was_mocha_mae) {
      Object <- methods::as(Object, "MochaTileResults")
    } else if (was_mocha_stm) {
      Object <- methods::as(Object, "MochaSampleTileMatrix")
    }
    return(Object)
  }

  stop("subsetBy not recognized.")
}

#' @keywords internal
#' @noRd
.subset_mocha_by_sample_ids <- function(Object, samples) {
  was_mocha_mae <- methods::is(Object, "MochaTileResults")
  was_mocha_stm <- methods::is(Object, "MochaSampleTileMatrix")
  summarizedData <- S4Vectors::metadata(Object)$summarizedData
  sampleData <- SummarizedExperiment::colData(Object)
  missing <- setdiff(samples, rownames(sampleData))
  if (length(missing) > 0) {
    stop(
      "These samples were not found in colData: ",
      paste(missing, collapse = ", ")
    )
  }
  Object <- Object[, samples]
  Object@metadata$summarizedData <- summarizedData[, samples, drop = FALSE]
  if (was_mocha_mae) {
    Object <- methods::as(Object, "MochaTileResults")
  } else if (was_mocha_stm) {
    Object <- methods::as(Object, "MochaSampleTileMatrix")
  }
  Object
}

#' @keywords internal
#' @noRd
.mocha_subset_dispatch <- function(x,
                                   cells = NULL,
                                   samples = NULL,
                                   subsetBy = NULL,
                                   groupList = NULL,
                                   subsetPeaks = TRUE,
                                   removeNA = TRUE,
                                   verbose = FALSE) {
  if (!is.null(cells)) {
    return(.subset_mocha_object(
      x,
      subsetBy = "celltype",
      groupList = cells,
      subsetPeaks = subsetPeaks,
      removeNA = removeNA,
      verbose = verbose
    ))
  }
  if (!is.null(samples)) {
    return(.subset_mocha_by_sample_ids(x, samples))
  }
  if (!is.null(subsetBy) && !is.null(groupList)) {
    return(.subset_mocha_object(
      x,
      subsetBy = subsetBy,
      groupList = groupList,
      subsetPeaks = subsetPeaks,
      removeNA = removeNA,
      verbose = verbose
    ))
  }
  stop("Provide `cells`, `samples`, or both `subsetBy` and `groupList`.")
}

setMethod(
  "subset",
  signature(x = "MochaTileResults"),
  function(x,
           subset = NULL,
           select,
           cells = NULL,
           samples = NULL,
           subsetBy = NULL,
           groupList = NULL,
           removeNA = TRUE,
           verbose = FALSE,
           ...) {
    .mocha_subset_dispatch(
      x,
      cells = cells,
      samples = samples,
      subsetBy = subsetBy,
      groupList = groupList,
      subsetPeaks = FALSE,
      removeNA = removeNA,
      verbose = verbose
    )
  }
)

setMethod(
  "subset",
  signature(x = "MochaSampleTileMatrix"),
  function(x,
           subset = NULL,
           select,
           cells = NULL,
           samples = NULL,
           subsetBy = NULL,
           groupList = NULL,
           subsetPeaks = TRUE,
           removeNA = TRUE,
           verbose = FALSE,
           ...) {
    .mocha_subset_dispatch(
      x,
      cells = cells,
      samples = samples,
      subsetBy = subsetBy,
      groupList = groupList,
      subsetPeaks = subsetPeaks,
      removeNA = removeNA,
      verbose = verbose
    )
  }
)

setMethod(
  "[",
  signature(x = "MochaTileResults", i = "ANY", j = "ANY"),
  function(x, i, j, ..., drop = TRUE) {
    dots <- list(...)
    k <- if ("k" %in% names(dots)) dots[["k"]] else NULL
    out <- if (!is.null(k)) {
      callNextMethod(x, i, j, k = k, drop = drop)
    } else {
      callNextMethod(x, i, j, drop = drop)
    }
    if (!is.null(S4Vectors::metadata(x)$summarizedData)) {
      sd <- S4Vectors::metadata(x)$summarizedData
      if (!missing(j) && !is.null(j)) {
        sd <- sd[, colnames(out), drop = FALSE]
      }
      if (!is.null(k)) {
        sd <- sd[k, , drop = FALSE]
      }
      out@metadata$summarizedData <- sd
    }
    methods::as(out, "MochaTileResults")
  }
)

setMethod(
  "[",
  signature(x = "MochaSampleTileMatrix", i = "ANY", j = "ANY"),
  function(x, i, j, ..., drop = TRUE) {
    out <- callNextMethod()
    if (!missing(j) && !is.null(j) && !is.null(S4Vectors::metadata(x)$summarizedData)) {
      out@metadata$summarizedData <- S4Vectors::metadata(x)$summarizedData[
        ,
        colnames(out),
        drop = FALSE
      ]
    }
    methods::as(out, "MochaSampleTileMatrix")
  }
)
