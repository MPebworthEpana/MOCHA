#' Import SnapATAC2 AnnData / AnnDataSet for MOCHA peak calling
#'
#' Converts a SnapATAC2 `AnnData` or `AnnDataSet` object (or path to on-disk
#' `.h5ad`) into MOCHA-native inputs for [callOpenTiles()]. Requires optional
#' Python packages via **reticulate** (`anndata` preferred; `snapatac2` for
#' `AnnDataSet` objects).
#'
#' @param x A reticulate Python object (`AnnData` / `AnnDataSet`), a character
#'   path to an `.h5ad` file, or a named list of paths (names used as sample IDs
#'   when `sampleColumn` is absent from `obs`).
#' @param cellPopLabel Column in `obs` containing cell population / cluster labels.
#' @param sampleColumn Column in `obs` with biological sample IDs. Default `NULL`:
#'   uses `"sample"` for `AnnDataSet`, otherwise must be supplied for `AnnData`.
#' @param cellCol Name of the metadata column stamped onto each fragment GRanges
#'   (default `"RG"`, matching MOCHA / ArchR conventions).
#' @param fragmentsLayer Which `.obsm` layer to use: `"auto"`, `"fragment_paired"`,
#'   or `"fragment_single"`.
#' @param chromPrefix Chromosome name harmonization: `"auto"` (no change),
#'   `"add"` (prepend `chr` if missing), `"strip"`, or `"none"`.
#' @param pythonModule Python module for I/O: `"auto"`, `"anndata"`, or
#'   `"snapatac2"`.
#' @param verbose Logical; print progress messages.
#'
#' @return A named list with elements suitable for splatting into
#'   [callOpenTiles()]:
#'   \describe{
#'     \item{ATACFragments}{Sample-level \code{GRangesList}.}
#'     \item{cellColData}{\code{data.frame} with rownames = cell barcodes.}
#'     \item{cellCol}{Character; column name on fragment metadata.}
#'     \item{seqinfo}{\code{Seqinfo} from SnapATAC2 \code{reference_sequences}.}
#'   }
#'   You must still provide \code{genome} and \code{blackList} to
#'   [callOpenTiles()].
#'
#' @details
#' SnapATAC2 stores fragments in `.obsm` as a CSR-like sparse matrix (rows =
#' cells, columns = global genomic start positions, values = fragment length).
#' Decode offsets use `uns['reference_sequences']`. See SnapATAC2 documentation:
#' \url{https://scverse.org/SnapATAC2/}.
#'
#' @keywords core
#' @export
import_snap_atac <- function(
    x,
    cellPopLabel,
    sampleColumn = NULL,
    cellCol = "RG",
    fragmentsLayer = c("auto", "fragment_paired", "fragment_single"),
    chromPrefix = c("auto", "add", "strip", "none"),
    pythonModule = c("auto", "anndata", "snapatac2"),
    verbose = FALSE) {
  fragmentsLayer <- match.arg(fragmentsLayer)
  chromPrefix <- match.arg(chromPrefix)
  pythonModule <- match.arg(pythonModule)

  if (missing(cellPopLabel) || !nzchar(cellPopLabel)) {
    stop("'cellPopLabel' must be a non-empty string.", call. = FALSE)
  }

  intermediate <- if (is.character(x) || is.list(x)) {
    .extract_snap_from_path(
      x,
      pythonModule = pythonModule,
      fragmentsLayer = fragmentsLayer,
      verbose = verbose
    )
  } else if (.is_reticulate_object(x)) {
    .extract_snap_from_object(
      x,
      pythonModule = pythonModule,
      fragmentsLayer = fragmentsLayer,
      verbose = verbose
    )
  } else {
    stop(
      "'x' must be a SnapATAC2 AnnData/AnnDataSet object (reticulate proxy), ",
      "a file path, or a named list of paths.",
      call. = FALSE
    )
  }

  objectType <- intermediate$object_type
  if (is.null(sampleColumn)) {
    if (identical(objectType, "AnnDataSet")) {
      sampleColumn <- "sample"
    } else {
      stop(
        "'sampleColumn' is required for AnnData input. ",
        "Provide the column in obs with biological sample IDs.",
        call. = FALSE
      )
    }
  }

  bundle <- .assembleMOCHAFromSnapATAC2(
    intermediate = intermediate,
    cellPopLabel = cellPopLabel,
    sampleColumn = sampleColumn,
    cellCol = cellCol,
    fragmentsLayer = fragmentsLayer,
    chromPrefix = chromPrefix,
    verbose = verbose
  )

  .validate_snap_import_bundle(
    bundle,
    cellPopLabel = cellPopLabel,
    sampleColumn = sampleColumn,
    cellCol = cellCol
  )

  bundle
}

# --- reticulate / Python guards ---------------------------------------------

#' @noRd
.require_reticulate <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop(
      "Package 'reticulate' is required for import_snap_atac(). ",
      "Install with install.packages('reticulate').",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' @noRd
.is_reticulate_object <- function(x) {
  inherits(x, "python.builtin.object") ||
    inherits(x, "environment") && isTRUE(tryCatch(reticulate::py_has_attr(x, "obs"), error = function(e) FALSE))
}

#' @noRd
.resolve_python_module <- function(pythonModule = c("auto", "anndata", "snapatac2")) {
  pythonModule <- match.arg(pythonModule)
  .require_reticulate()
  if (pythonModule == "anndata") {
    if (!reticulate::py_module_available("anndata")) {
      stop("Python package 'anndata' is not available.", call. = FALSE)
    }
    return("anndata")
  }
  if (pythonModule == "snapatac2") {
    if (!reticulate::py_module_available("snapatac2")) {
      stop("Python package 'snapatac2' is not available.", call. = FALSE)
    }
    return("snapatac2")
  }
  if (reticulate::py_module_available("anndata")) {
    return("anndata")
  }
  if (reticulate::py_module_available("snapatac2")) {
    return("snapatac2")
  }
  stop(
    "No supported Python module found. Install 'anndata' (preferred) or ",
    "'snapatac2' in the active Python environment used by reticulate.",
    call. = FALSE
  )
}

#' @noRd
.py_class_name <- function(x) {
  reticulate::py_to_r(x$`__class__`$`__name__`)
}

#' @noRd
.detect_snap_object_type <- function(x) {
  cls <- .py_class_name(x)
  if (grepl("AnnDataSet", cls, fixed = TRUE)) {
    return("AnnDataSet")
  }
  if (grepl("AnnData", cls, fixed = TRUE)) {
    return("AnnData")
  }
  stop(
    "Expected a SnapATAC2 AnnData or AnnDataSet object; got Python class '",
    cls, "'.",
    call. = FALSE
  )
}

# --- extraction from Python objects -----------------------------------------

#' @noRd
.extract_snap_from_object <- function(
    x,
    pythonModule = "auto",
    fragmentsLayer = "auto",
    verbose = FALSE) {
  .require_reticulate()
  objectType <- .detect_snap_object_type(x)
  if (verbose) {
    message("Detected SnapATAC2 object type: ", objectType)
  }
  .extract_snap_intermediate(
    x,
    object_type = objectType,
    fragmentsLayer = fragmentsLayer,
    verbose = verbose
  )
}

#' @noRd
.extract_snap_from_path <- function(
    x,
    pythonModule = "auto",
    fragmentsLayer = "auto",
    verbose = FALSE) {
  .require_reticulate()
  mod <- .resolve_python_module(pythonModule)

  if (is.list(x)) {
    if (is.null(names(x)) || any(!nzchar(names(x)))) {
      stop("Named list of paths must have non-empty names (used as sample IDs).", call. = FALSE)
    }
    adatas <- lapply(x, function(p) .load_snap_path(p, mod = mod))
    return(.merge_snap_adatas(
      adatas,
      sample_ids = names(x),
      fragmentsLayer = fragmentsLayer,
      verbose = verbose
    ))
  }

  if (!is.character(x) || length(x) != 1L || !nzchar(x)) {
    stop("'x' must be a single file path or a named list of paths.", call. = FALSE)
  }
  if (!file.exists(x)) {
    stop("File not found: ", x, call. = FALSE)
  }

  obj <- .load_snap_path(x, mod = mod)
  objectType <- .detect_snap_object_type(obj)
  if (identical(objectType, "AnnDataSet") && identical(mod, "anndata")) {
    stop(
      "AnnDataSet paths require the snapatac2 Python package. ",
      "Install snapatac2 or pass the loaded object instead of a path.",
      call. = FALSE
    )
  }
  .extract_snap_intermediate(
    obj,
    object_type = objectType,
    fragmentsLayer = fragmentsLayer,
    verbose = verbose
  )
}

#' @noRd
.load_snap_path <- function(path, mod = "anndata") {
  if (identical(mod, "anndata")) {
    ad <- reticulate::import("anndata", delay_load = TRUE)
    return(ad$read_h5ad(path))
  }
  snap <- reticulate::import("snapatac2", delay_load = TRUE)
  if (grepl("\\.h5ad$", path, ignore.case = TRUE)) {
    return(snap$read(path))
  }
  snap$read_dataset(path)
}

#' @noRd
.merge_snap_adatas <- function(
    adatas,
    sample_ids,
    fragmentsLayer = "auto",
    verbose = FALSE) {
  if (length(adatas) < 1L) {
    stop("No AnnData objects to merge.", call. = FALSE)
  }
  ref <- .py_reference_sequences(adatas[[1L]])
  layer_name <- .py_resolve_fragments_layer(adatas[[1L]], fragmentsLayer = fragmentsLayer)

  obs_list <- vector("list", length(adatas))
  csr_parts <- vector("list", length(adatas))
  row_offset <- 0L
  obs_names <- character()

  for (i in seq_along(adatas)) {
    ad <- adatas[[i]]
    df <- .py_obs_to_df(ad)
    nm <- .py_obs_names(ad)
    if (length(nm) != nrow(df)) {
      nm <- rownames(df)
    }
    if (is.null(nm) || !length(nm)) {
      nm <- paste0("cell_", seq_len(nrow(df)), "_", sample_ids[[i]])
    }
    nm <- make.unique(as.character(nm), sep = "_")
    rownames(df) <- nm
    df$`.mocha_sample` <- sample_ids[[i]]
    obs_list[[i]] <- df
    obs_names <- c(obs_names, nm)

    csr <- .py_extract_csr(ad, layer = layer_name, fragmentsLayer = "auto", verbose = FALSE)
    csr$rows <- csr$rows + row_offset
    csr_parts[[i]] <- csr
    row_offset <- row_offset + length(nm)
  }

  obs_df <- do.call(rbind, obs_list)
  rownames(obs_df) <- obs_names

  csr_rows <- unlist(lapply(csr_parts, `[[`, "rows"), use.names = FALSE)
  csr_cols <- unlist(lapply(csr_parts, `[[`, "cols"), use.names = FALSE)
  csr_values <- unlist(lapply(csr_parts, `[[`, "values"), use.names = FALSE)

  if (verbose) {
    message("Merged ", length(adatas), " AnnData files into ", nrow(obs_df), " cells.")
  }

  list(
    object_type = "AnnData",
    obs_df = obs_df,
    obs_names = obs_names,
    reference_sequences = ref,
    csr_rows = csr_rows,
    csr_cols = csr_cols,
    csr_values = csr_values,
    layer_name = layer_name,
    sample_column_override = ".mocha_sample"
  )
}

#' @noRd
.extract_snap_intermediate <- function(
    x,
    object_type,
    fragmentsLayer = "auto",
    verbose = FALSE) {
  obs_df <- .py_obs_to_df(x)
  obs_names <- .py_obs_names(x)
  if (length(obs_names) != nrow(obs_df)) {
    obs_names <- rownames(obs_df)
  }
  rownames(obs_df) <- obs_names

  ref <- .py_reference_sequences(x)
  csr <- .py_extract_csr(x, layer = NULL, fragmentsLayer = fragmentsLayer, verbose = verbose)

  out <- list(
    object_type = object_type,
    obs_df = obs_df,
    obs_names = obs_names,
    reference_sequences = ref,
    csr_rows = csr$rows,
    csr_cols = csr$cols,
    csr_values = csr$values,
    layer_name = csr$layer_name,
    sample_column_override = NULL
  )
  if (verbose) {
    message(
      "Extracted ", length(out$csr_values), " non-zero fragment entries from layer '",
      out$layer_name, "'."
    )
  }
  out
}

#' @noRd
.py_obs_to_df <- function(x) {
  obs <- x$obs
  df <- tryCatch(
    reticulate::py_to_r(obs),
    error = function(e) {
      pd <- reticulate::import("pandas", delay_load = TRUE)
      reticulate::py_to_r(pd$DataFrame(obs))
    }
  )
  if (!is.data.frame(df)) {
    df <- as.data.frame(df, stringsAsFactors = FALSE)
  }
  df
}

#' @noRd
.py_obs_names <- function(x) {
  nm <- tryCatch(
    reticulate::py_to_r(x$obs_names),
    error = function(e) character()
  )
  as.character(nm)
}

#' @noRd
.py_reference_sequences <- function(x) {
  uns <- x$uns
  if (is.null(uns)) {
    stop(
      "SnapATAC2 object is missing uns slot.",
      call. = FALSE
    )
  }
  ref <- tryCatch(
    uns[["reference_sequences"]],
    error = function(e) {
      tryCatch(
        uns$reference_sequences,
        error = function(e2) NULL
      )
    }
  )
  if (is.null(ref)) {
    stop(
      "Could not read uns['reference_sequences']. ",
      "See https://scverse.org/SnapATAC2/",
      call. = FALSE
    )
  }
  ref_r <- reticulate::py_to_r(ref)
  if (is.data.frame(ref_r)) {
    if (all(c("seqname", "length") %in% names(ref_r))) {
      return(data.frame(
        seqname = as.character(ref_r$seqname),
        length = as.integer(ref_r$length),
        stringsAsFactors = FALSE
      ))
    }
    if (ncol(ref_r) >= 2L) {
      return(data.frame(
        seqname = as.character(ref_r[[1L]]),
        length = as.integer(ref_r[[2L]]),
        stringsAsFactors = FALSE
      ))
    }
  }
  if (is.list(ref_r) && !is.null(names(ref_r))) {
    return(data.frame(
      seqname = names(ref_r),
      length = as.integer(unlist(ref_r, use.names = FALSE)),
      stringsAsFactors = FALSE
    ))
  }
  stop("Unrecognized format for reference_sequences.", call. = FALSE)
}

#' @noRd
.py_obsm_has_layer <- function(obsm, layer) {
  mat <- tryCatch(
    obsm[[layer]],
    error = function(e) NULL
  )
  !is.null(mat)
}

#' @noRd
.py_list_obsm_layers <- function(obsm) {
  keys <- tryCatch(
    reticulate::py_to_r(obsm$keys()),
    error = function(e) character()
  )
  keys <- as.character(keys)
  if (length(keys) > 0L && any(nzchar(keys))) {
    return(keys)
  }
  # AnnData AxisArrays KeysView may not coerce cleanly; probe known layers.
  candidates <- c("fragment_paired", "fragment_single")
  keys <- candidates[vapply(candidates, function(k) .py_obsm_has_layer(obsm, k), logical(1))]
  keys
}

#' @noRd
.py_resolve_fragments_layer <- function(x, fragmentsLayer = "auto") {
  obsm <- x$obsm
  if (is.null(obsm)) {
    stop("SnapATAC2 object has no obsm slot.", call. = FALSE)
  }
  keys <- .py_list_obsm_layers(obsm)
  if (fragmentsLayer != "auto") {
    if (!.py_obsm_has_layer(obsm, fragmentsLayer)) {
      stop(
        "Requested fragmentsLayer '", fragmentsLayer, "' not found in obsm. ",
        "Available: ", paste(keys, collapse = ", "),
        call. = FALSE
      )
    }
    return(fragmentsLayer)
  }
  if (.py_obsm_has_layer(obsm, "fragment_paired")) {
    return("fragment_paired")
  }
  if (.py_obsm_has_layer(obsm, "fragment_single")) {
    return("fragment_single")
  }
  stop(
    "No fragment layer found in obsm (expected 'fragment_paired' or ",
    "'fragment_single'). Import fragments with snapatac2.pp.import_fragments(). ",
    "Available obsm keys: ", paste(keys, collapse = ", "),
    call. = FALSE
  )
}

#' @noRd
.extract_csr_from_r_matrix <- function(mat) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    return(NULL)
  }
  if (inherits(mat, "dgRMatrix")) {
    p <- mat@p
    j <- mat@j
    x <- mat@x
    n_rows <- nrow(mat)
    rows <- integer()
    cols <- integer()
    values <- numeric()
    for (r in seq_len(n_rows)) {
      start <- p[r] + 1L
      end <- p[r + 1L]
      if (end >= start) {
        idx <- seq.int(start, end)
        rows <- c(rows, rep.int(r, length(idx)))
        cols <- c(cols, as.integer(j[idx]))
        values <- c(values, as.numeric(x[idx]))
      }
    }
    return(list(rows = rows, cols = cols, values = values))
  }
  if (inherits(mat, "dgCMatrix")) {
    return(.extract_csr_from_r_matrix(Matrix::t(mat)))
  }
  NULL
}

#' @noRd
.py_extract_csr <- function(x, layer = NULL, fragmentsLayer = "auto", verbose = FALSE) {
  layer_name <- if (is.null(layer)) {
    .py_resolve_fragments_layer(x, fragmentsLayer = fragmentsLayer)
  } else {
    layer
  }
  mat <- x$obsm[[layer_name]]
  if (is.null(mat)) {
    stop("obsm layer '", layer_name, "' is NULL.", call. = FALSE)
  }

  r_csr <- .extract_csr_from_r_matrix(mat)
  if (!is.null(r_csr)) {
    r_csr$layer_name <- layer_name
    return(r_csr)
  }

  sp <- reticulate::import("scipy.sparse", delay_load = TRUE)
  coo <- tryCatch(
    {
      m <- reticulate::py_get_attr(sp, "coo_matrix")(mat)
      list(
        rows = as.integer(reticulate::py_to_r(reticulate::py_get_attr(m, "row"))) + 1L,
        cols = as.integer(reticulate::py_to_r(reticulate::py_get_attr(m, "col"))),
        values = as.numeric(reticulate::py_to_r(reticulate::py_get_attr(m, "data")))
      )
    },
    error = function(e) {
      .py_extract_csr_manual(mat)
    }
  )
  coo$layer_name <- layer_name
  coo
}

#' @noRd
.py_extract_csr_manual <- function(mat) {
  indptr <- as.integer(reticulate::py_to_r(reticulate::py_get_attr(mat, "indptr")))
  indices <- as.integer(reticulate::py_to_r(reticulate::py_get_attr(mat, "indices")))
  data <- as.numeric(reticulate::py_to_r(reticulate::py_get_attr(mat, "data")))
  n_rows <- length(indptr) - 1L
  rows <- integer()
  cols <- integer()
  values <- numeric()
  for (r in seq_len(n_rows)) {
    start <- indptr[r] + 1L
    end <- indptr[r + 1L]
    if (end >= start) {
      idx <- seq.int(start, end)
      rows <- c(rows, rep.int(r, length(idx)))
      cols <- c(cols, indices[idx])
      values <- c(values, data[idx])
    }
  }
  list(rows = rows, cols = cols, values = values)
}

# --- pure-R assembly --------------------------------------------------------

#' @noRd
.build_reference_offsets <- function(reference_sequences) {
  ref <- reference_sequences
  if (nrow(ref) == 0L) {
    stop("reference_sequences is empty.", call. = FALSE)
  }
  offsets <- c(0L, cumsum(as.integer(ref$length)))
  list(
    seqnames = as.character(ref$seqname),
    lengths = as.integer(ref$length),
    offsets = offsets
  )
}

#' @noRd
.global_col_to_locus <- function(col_global, ref_info) {
  col_global <- as.integer(col_global)
  n_chr <- length(ref_info$seqnames)
  chrom_idx <- rep.int(NA_integer_, length(col_global))
  for (i in seq_len(n_chr)) {
    lo <- ref_info$offsets[i]
    hi <- ref_info$offsets[i + 1L] - 1L
    in_chr <- col_global >= lo & col_global <= hi
    chrom_idx[in_chr] <- i
  }
  if (any(is.na(chrom_idx))) {
    stop(
      "Some fragment column indices fall outside reference_sequences. ",
      "Check genome build consistency.",
      call. = FALSE
    )
  }
  start <- col_global - ref_info$offsets[chrom_idx] + 1L
  list(
    seqnames = ref_info$seqnames[chrom_idx],
    start = start
  )
}

#' @noRd
.normalize_chrom_names <- function(seqnames, chromPrefix = "auto") {
  chromPrefix <- match.arg(chromPrefix, c("auto", "add", "strip", "none"))
  if (chromPrefix %in% c("auto", "none")) {
    return(seqnames)
  }
  if (chromPrefix == "add") {
    needs <- !grepl("^chr", seqnames, ignore.case = TRUE)
    seqnames[needs] <- paste0("chr", seqnames[needs])
    return(seqnames)
  }
  sub("^chr", "", seqnames, ignore.case = TRUE)
}

#' @noRd
.decode_csr_to_granges <- function(
    csr_rows,
    csr_cols,
    csr_values,
    obs_names,
    layer_name,
    chromPrefix = "auto",
    ref_info) {
  if (length(csr_rows) == 0L) {
    return(GenomicRanges::GRanges())
  }
  locus <- .global_col_to_locus(csr_cols, ref_info)
  seqnames <- .normalize_chrom_names(locus$seqnames, chromPrefix = chromPrefix)
  is_single <- identical(layer_name, "fragment_single")

  if (is_single) {
    strand <- ifelse(csr_values > 0, "+", "-")
    width <- abs(csr_values)
  } else {
    strand <- rep.int("*", length(csr_values))
    width <- csr_values
  }
  width[width < 1] <- 1L
  end <- locus$start + as.integer(width) - 1L

  GenomicRanges::GRanges(
    seqnames = seqnames,
    ranges = IRanges::IRanges(start = locus$start, end = end),
    strand = strand,
    cell_id = obs_names[csr_rows]
  )
}

#' @noRd
.assembleMOCHAFromSnapATAC2 <- function(
    intermediate,
    cellPopLabel,
    sampleColumn,
    cellCol = "RG",
    fragmentsLayer = "auto",
    chromPrefix = "auto",
    verbose = FALSE) {
  obs_df <- intermediate$obs_df
  obs_names <- intermediate$obs_names

  if (!is.null(intermediate$sample_column_override)) {
    sampleColumn <- intermediate$sample_column_override
  }

  if (!(cellPopLabel %in% colnames(obs_df))) {
    stop(
      "cellPopLabel '", cellPopLabel, "' not found in obs. ",
      "Columns: ", paste(colnames(obs_df), collapse = ", "),
      call. = FALSE
    )
  }
  if (!(sampleColumn %in% colnames(obs_df))) {
    stop(
      "sampleColumn '", sampleColumn, "' not found in obs. ",
      "Columns: ", paste(colnames(obs_df), collapse = ", "),
      call. = FALSE
    )
  }

  layer_name <- intermediate$layer_name
  if (fragmentsLayer != "auto" && fragmentsLayer != layer_name) {
    stop(
      "Requested fragmentsLayer '", fragmentsLayer, "' but intermediate has '",
      layer_name, "'. Re-extract with matching layer.",
      call. = FALSE
    )
  }

  ref_info <- .build_reference_offsets(intermediate$reference_sequences)
  all_gr <- .decode_csr_to_granges(
    csr_rows = intermediate$csr_rows,
    csr_cols = intermediate$csr_cols,
    csr_values = intermediate$csr_values,
    obs_names = obs_names,
    layer_name = layer_name,
    chromPrefix = chromPrefix,
    ref_info = ref_info
  )

  if (length(all_gr) == 0L) {
    stop("No fragments decoded from SnapATAC2 object.", call. = FALSE)
  }

  sample_ids <- as.character(obs_df[[sampleColumn]])
  names(sample_ids) <- obs_names
  GenomicRanges::mcols(all_gr)[[cellCol]] <- GenomicRanges::mcols(all_gr)$cell_id
  GenomicRanges::mcols(all_gr)$cell_id <- NULL

  cellColData <- obs_df
  rownames(cellColData) <- obs_names

  unique_samples <- unique(sample_ids)
  frag_list <- stats::setNames(vector("list", length(unique_samples)), unique_samples)

  for (s in unique_samples) {
    cells_in_sample <- names(sample_ids)[sample_ids == s]
    frag_list[[s]] <- all_gr[GenomicRanges::mcols(all_gr)[[cellCol]] %in% cells_in_sample]
  }

  empty_samples <- names(frag_list)[vapply(frag_list, length, integer(1)) == 0L]
  if (length(empty_samples) > 0L) {
    stop(
      "No fragments for sample(s): ",
      paste(empty_samples, collapse = ", "),
      call. = FALSE
    )
  }

  ATACFragments <- GenomicRanges::GRangesList(frag_list)

  seqinfo <- GenomeInfoDb::Seqinfo(
    seqnames = ref_info$seqnames,
    seqlengths = ref_info$lengths
  )

  if (verbose) {
    message(
      "Assembled ", length(ATACFragments), " sample-level fragment sets (",
      sum(lengths(ATACFragments)), " total fragments)."
    )
  }

  list(
    ATACFragments = ATACFragments,
    cellColData = cellColData,
    cellCol = cellCol,
    seqinfo = seqinfo
  )
}

#' @noRd
.validate_snap_import_bundle <- function(
    bundle,
    cellPopLabel,
    sampleColumn,
    cellCol) {
  if (!is.list(bundle) || !all(c("ATACFragments", "cellColData", "cellCol") %in% names(bundle))) {
    stop("Invalid import bundle structure.", call. = FALSE)
  }
  if (!methods::is(bundle$ATACFragments, "GRangesList")) {
    stop("ATACFragments must be a GRangesList.", call. = FALSE)
  }
  if (is.null(names(bundle$ATACFragments)) || any(!nzchar(names(bundle$ATACFragments)))) {
    stop("ATACFragments must be a named GRangesList (sample IDs).", call. = FALSE)
  }
  cellColData <- bundle$cellColData
  if (is.null(rownames(cellColData)) || !any(nzchar(rownames(cellColData)))) {
    stop("cellColData must have rownames matching cell barcodes.", call. = FALSE)
  }
  if (!(cellPopLabel %in% colnames(cellColData))) {
    stop("cellColData must contain cellPopLabel column.", call. = FALSE)
  }
  if (!(sampleColumn %in% colnames(cellColData))) {
    stop("cellColData must contain sampleColumn.", call. = FALSE)
  }
  if (!(cellCol %in% colnames(GenomicRanges::mcols(bundle$ATACFragments[[1L]])))) {
    stop(
      "Fragment GRanges must contain metadata column '", cellCol, "'.",
      call. = FALSE
    )
  }

  meta_samples <- unique(as.character(cellColData[[sampleColumn]]))
  if (!all(names(bundle$ATACFragments) %in% meta_samples)) {
    stop(
      "ATACFragments sample names must match cellColData[[", sampleColumn, "]].",
      call. = FALSE
    )
  }

  frag_cells <- unique(as.character(GenomicRanges::mcols(bundle$ATACFragments[[1L]])[[cellCol]]))
  unknown <- setdiff(frag_cells, rownames(cellColData))
  if (length(unknown) > 0L) {
    stop(
      "Fragments contain cell IDs not in cellColData: ",
      paste(head(unknown, 5L), collapse = ", "),
      if (length(unknown) > 5L) paste0(" (and ", length(unknown) - 5L, " more)") else "",
      call. = FALSE
    )
  }

  invisible(bundle)
}
