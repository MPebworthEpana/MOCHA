#' Index a fragment TSV for Signac (bgzip + tabix).
#' @noRd
.index_fragment_tsv <- function(frag_tsv) {
  if (!nzchar(Sys.which("bgzip")) || !nzchar(Sys.which("tabix"))) {
    return(NULL)
  }
  frag_gz <- paste0(frag_tsv, ".gz")
  status_bgzip <- system2("bgzip", args = c("-f", shQuote(frag_tsv)), stdout = FALSE, stderr = FALSE)
  if (status_bgzip != 0 || !file.exists(frag_gz)) {
    return(NULL)
  }
  status_tabix <- system2(
    "tabix",
    args = c("-s", "1", "-b", "2", "-e", "3", shQuote(frag_gz)),
    stdout = FALSE,
    stderr = FALSE
  )
  if (status_tabix != 0 || !file.exists(paste0(frag_gz, ".tbi"))) {
    return(NULL)
  }
  frag_gz
}

#' Build cells_map (Seurat cell id -> file barcode) from meta and fragment lines.
#' @noRd
.build_cells_map <- function(meta, fragment_lines) {
  file_bcs <- vapply(
    strsplit(fragment_lines, "\t", fixed = TRUE),
    function(x) x[[4]],
    character(1)
  )
  cell_ids <- rownames(meta)
  n_cells <- length(cell_ids)

  # Signac Fragment@cells / CreateFragmentObject(cells=): names = Seurat cell IDs,
  # values = barcodes in the fragment file.
  if (length(file_bcs) == n_cells) {
    return(stats::setNames(file_bcs, cell_ids))
  }
  if (length(file_bcs) == 1L && n_cells > 1L) {
    return(stats::setNames(rep(file_bcs, n_cells), cell_ids))
  }
  if (n_cells <= length(file_bcs)) {
    return(stats::setNames(file_bcs[seq_len(n_cells)], cell_ids[seq_len(n_cells)]))
  }
  stats::setNames(
    rep(file_bcs, length.out = n_cells),
    cell_ids
  )
}

#' Path to bundled mini indexed fragments fixture (fallback when tabix missing).
#' @noRd
.mini_fragment_fixture_path <- function() {
  fixture_gz <- testthat::test_path("fixtures", "seurat", "mini_fragments.tsv.gz")
  if (file.exists(fixture_gz) && file.exists(paste0(fixture_gz, ".tbi"))) {
    return(fixture_gz)
  }
  NULL
}

#' Build a minimal Seurat + Signac object for MOCHA tests.
#' @noRd
make_test_seurat_chromatin <- function(
    fragment_lines,
    meta,
    assay_name = "ATAC",
    cells_map = NULL) {
  if (!requireNamespace("Seurat", quietly = TRUE) ||
      !requireNamespace("Signac", quietly = TRUE)) {
    return(NULL)
  }

  if (length(fragment_lines) == 0) {
    fragment_lines <- "chr1\t760101\t760110\tAAAC\t1"
  }

  frag_path <- .mini_fragment_fixture_path()
  temp_files <- character(0)

  if (!is.null(frag_path) && length(fragment_lines) == 0) {
    frag_tbl <- utils::read.table(
      gzfile(frag_path),
      sep = "\t",
      header = FALSE,
      stringsAsFactors = FALSE,
      comment.char = ""
    )
    fragment_lines <- apply(frag_tbl, 1, paste, collapse = "\t")
  }

  if (is.null(frag_path)) {
    frag_tsv <- tempfile(fileext = ".tsv")
    temp_files <- c(temp_files, frag_tsv, paste0(frag_tsv, ".gz"), paste0(frag_tsv, ".gz.tbi"))
    writeLines(fragment_lines, con = frag_tsv)
    frag_path <- .index_fragment_tsv(frag_tsv)
    if (is.null(frag_path)) {
      message(
        "make_test_seurat_chromatin: tabix/bgzip not available and no fixture found. ",
        "Install htslib in mocha-test or add tests/testthat/fixtures/seurat/mini_fragments.tsv.gz"
      )
      return(NULL)
    }
  }

  on.exit(unlink(temp_files[temp_files != frag_path], force = TRUE), add = TRUE)

  n_cells <- nrow(meta)
  # Signac CreateChromatinAssay requires >= 2 peaks (single-row matrices fail internally).
  peak_starts <- c(760101L, 760111L, 760121L)
  peak_ends <- c(760110L, 760120L, 760130L)
  n_peaks <- length(peak_starts)
  counts <- Matrix::sparseMatrix(
    i = rep(seq_len(n_peaks), each = n_cells),
    j = rep(seq_len(n_cells), times = n_peaks),
    x = 1,
    dims = c(n_peaks, n_cells)
  )
  rownames(counts) <- paste0(
    "chr1:", peak_starts, "-", peak_ends
  )
  colnames(counts) <- rownames(meta)

  tryCatch({
    if (is.null(cells_map)) {
      cells_map <- .build_cells_map(meta, fragment_lines)
    }

    frags <- Signac::CreateFragmentObject(
      path = frag_path,
      cells = cells_map
    )
    chrom_assay <- Signac::CreateChromatinAssay(
      counts = counts,
      sep = c(":", "-"),
      fragments = frags,
      min.cells = 0,
      min.features = 0
    )

    obj <- Seurat::CreateSeuratObject(
      counts = chrom_assay,
      assay = assay_name,
      min.cells = 0,
      min.features = 0
    )
    obj[[assay_name]] <- chrom_assay
    obj@meta.data <- cbind(obj@meta.data, meta)
    rownames(obj@meta.data) <- rownames(meta)
    obj
  }, error = function(e) {
    message("make_test_seurat_chromatin failed: ", conditionMessage(e))
    NULL
  })
}
