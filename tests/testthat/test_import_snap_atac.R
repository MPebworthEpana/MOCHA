# Pure-R tests for SnapATAC2 import (no Python required for assembler/decoder)

.make_snap_intermediate <- function(
    obs_df,
    csr_rows,
    csr_cols,
    csr_values,
    layer_name = "fragment_paired",
    reference_sequences = data.frame(
      seqname = c("chr1", "chr2"),
      length = c(1000L, 1000L),
      stringsAsFactors = FALSE
    )) {
  obs_names <- rownames(obs_df)
  list(
    object_type = "AnnData",
    obs_df = obs_df,
    obs_names = obs_names,
    reference_sequences = reference_sequences,
    csr_rows = as.integer(csr_rows),
    csr_cols = as.integer(csr_cols),
    csr_values = as.numeric(csr_values),
    layer_name = layer_name,
    sample_column_override = NULL
  )
}

test_that("decode paired fragments and group by sample", {
  obs_df <- data.frame(
    Sample = c("s1", "s1", "s2"),
    cluster = c("A", "A", "B"),
    row.names = c("c1", "c2", "c3"),
    stringsAsFactors = FALSE
  )
  # col 0 -> chr1 start 1; col 500 -> chr1 start 501; col 1000 -> chr2 start 1
  intermediate <- .make_snap_intermediate(
    obs_df = obs_df,
    csr_rows = c(1L, 1L, 2L, 3L),
    csr_cols = c(0L, 10L, 500L, 1000L),
    csr_values = c(50, 60, 100, 80),
    layer_name = "fragment_paired"
  )

  bundle <- MOCHA:::.assembleMOCHAFromSnapATAC2(
    intermediate = intermediate,
    cellPopLabel = "cluster",
    sampleColumn = "Sample",
    cellCol = "RG",
    fragmentsLayer = "auto",
    chromPrefix = "auto",
    verbose = FALSE
  )

  expect_s4_class(bundle$ATACFragments, "GRangesList")
  expect_equal(names(bundle$ATACFragments), c("s1", "s2"))
  expect_equal(length(bundle$ATACFragments[["s1"]]), 3L)
  expect_equal(length(bundle$ATACFragments[["s2"]]), 1L)
  expect_true("RG" %in% colnames(GenomicRanges::mcols(bundle$ATACFragments[[1L]])))
  expect_equal(as.character(GenomicRanges::mcols(bundle$ATACFragments[[1L]])$RG), c("c1", "c1", "c2"))
  expect_equal(as.character(GenomicRanges::seqnames(bundle$ATACFragments[[1L]])), c("chr1", "chr1", "chr1"))
})

test_that("decode single-end fragments with strand from sign", {
  obs_df <- data.frame(
    Sample = "s1",
    cluster = "A",
    row.names = "c1",
    stringsAsFactors = FALSE
  )
  intermediate <- .make_snap_intermediate(
    obs_df = obs_df,
    csr_rows = c(1L, 1L),
    csr_cols = c(0L, 5L),
    csr_values = c(50, -40),
    layer_name = "fragment_single"
  )

  bundle <- MOCHA:::.assembleMOCHAFromSnapATAC2(
    intermediate = intermediate,
    cellPopLabel = "cluster",
    sampleColumn = "Sample",
    cellCol = "RG",
    fragmentsLayer = "auto",
    chromPrefix = "auto",
    verbose = FALSE
  )

  gr <- bundle$ATACFragments[[1L]]
  expect_equal(as.character(GenomicRanges::strand(gr)), c("+", "-"))
  expect_equal(GenomicRanges::start(gr)[1L], 1L)
  expect_equal(GenomicRanges::end(gr)[1L], 50L)
})

test_that("chrom decoding at chromosome boundary", {
  obs_df <- data.frame(
    Sample = "s1",
    cluster = "A",
    row.names = "c1",
    stringsAsFactors = FALSE
  )
  # chr1 length 1000 -> global col 999 is last base chr1; col 1000 is first base chr2
  intermediate <- .make_snap_intermediate(
    obs_df = obs_df,
    csr_rows = c(1L, 1L),
    csr_cols = c(999L, 1000L),
    csr_values = c(10, 20),
    layer_name = "fragment_paired"
  )

  bundle <- MOCHA:::.assembleMOCHAFromSnapATAC2(
    intermediate = intermediate,
    cellPopLabel = "cluster",
    sampleColumn = "Sample",
    cellCol = "RG",
    fragmentsLayer = "auto",
    chromPrefix = "auto",
    verbose = FALSE
  )

  gr <- bundle$ATACFragments[[1L]]
  expect_equal(as.character(GenomicRanges::seqnames(gr)), c("chr1", "chr2"))
  expect_equal(GenomicRanges::start(gr), c(1000L, 1L))
})

test_that("chromPrefix add prepends chr", {
  obs_df <- data.frame(
    Sample = "s1",
    cluster = "A",
    row.names = "c1",
    stringsAsFactors = FALSE
  )
  intermediate <- .make_snap_intermediate(
    obs_df = obs_df,
    csr_rows = 1L,
    csr_cols = 0L,
    csr_values = 50,
    layer_name = "fragment_paired",
    reference_sequences = data.frame(
      seqname = "1",
      length = 2000L,
      stringsAsFactors = FALSE
    )
  )

  bundle <- MOCHA:::.assembleMOCHAFromSnapATAC2(
    intermediate = intermediate,
    cellPopLabel = "cluster",
    sampleColumn = "Sample",
    cellCol = "RG",
    fragmentsLayer = "auto",
    chromPrefix = "add",
    verbose = FALSE
  )

  expect_equal(as.character(GenomicRanges::seqnames(bundle$ATACFragments[[1L]])), "chr1")
})

test_that("errors when cellPopLabel missing from obs", {
  obs_df <- data.frame(
    Sample = "s1",
    row.names = "c1",
    stringsAsFactors = FALSE
  )
  intermediate <- .make_snap_intermediate(
    obs_df = obs_df,
    csr_rows = 1L,
    csr_cols = 0L,
    csr_values = 50
  )

  expect_error(
    MOCHA:::.assembleMOCHAFromSnapATAC2(
      intermediate = intermediate,
      cellPopLabel = "cluster",
      sampleColumn = "Sample",
      cellCol = "RG",
      fragmentsLayer = "auto",
      chromPrefix = "auto",
      verbose = FALSE
    ),
    regexp = "cellPopLabel"
  )
})

test_that("errors when no fragments decoded", {
  obs_df <- data.frame(
    Sample = "s1",
    cluster = "A",
    row.names = "c1",
    stringsAsFactors = FALSE
  )
  intermediate <- .make_snap_intermediate(
    obs_df = obs_df,
    csr_rows = integer(),
    csr_cols = integer(),
    csr_values = numeric()
  )

  expect_error(
    MOCHA:::.assembleMOCHAFromSnapATAC2(
      intermediate = intermediate,
      cellPopLabel = "cluster",
      sampleColumn = "Sample",
      cellCol = "RG",
      fragmentsLayer = "auto",
      chromPrefix = "auto",
      verbose = FALSE
    ),
    regexp = "No fragments decoded"
  )
})

test_that("import_snap_atac errors without reticulate for object input", {
  skip_if(requireNamespace("reticulate", quietly = TRUE), "reticulate installed")
  expect_error(
    MOCHA::import_snap_atac(
      x = structure(list(), class = "python.builtin.object"),
      cellPopLabel = "cluster"
    ),
    regexp = "reticulate"
  )
})

test_that("import_snap_atac requires sampleColumn for AnnData-like intermediate", {
  skip_if_not(requireNamespace("reticulate", quietly = TRUE))

  obs_df <- data.frame(
    cluster = "A",
    row.names = "c1",
    stringsAsFactors = FALSE
  )
  intermediate <- .make_snap_intermediate(
    obs_df = obs_df,
    csr_rows = 1L,
    csr_cols = 0L,
    csr_values = 50
  )

  # Bypass Python extraction by testing entry validation via mock would need mockery;
  # test assemble path validates sample column when override absent
  expect_error(
    MOCHA:::.assembleMOCHAFromSnapATAC2(
      intermediate = intermediate,
      cellPopLabel = "cluster",
      sampleColumn = "Sample",
      cellCol = "RG",
      fragmentsLayer = "auto",
      chromPrefix = "auto",
      verbose = FALSE
    ),
    regexp = "sampleColumn"
  )
})

test_that("import bundle validates fragment cell IDs against cellColData", {
  obs_df <- data.frame(
    Sample = "s1",
    cluster = "A",
    row.names = "c1",
    stringsAsFactors = FALSE
  )
  intermediate <- .make_snap_intermediate(
    obs_df = obs_df,
    csr_rows = 1L,
    csr_cols = 0L,
    csr_values = 50
  )
  bundle <- MOCHA:::.assembleMOCHAFromSnapATAC2(
    intermediate = intermediate,
    cellPopLabel = "cluster",
    sampleColumn = "Sample",
    cellCol = "RG",
    fragmentsLayer = "auto",
    chromPrefix = "auto",
    verbose = FALSE
  )
  GenomicRanges::mcols(bundle$ATACFragments[[1L]])$RG <- "unknown_cell"
  expect_error(
    MOCHA:::.validate_snap_import_bundle(
      bundle,
      cellPopLabel = "cluster",
      sampleColumn = "Sample",
      cellCol = "RG"
    ),
    regexp = "not in cellColData"
  )
})

.make_py_snap_adata <- function(
    rows,
    cols,
    data,
    n_obs,
    n_cols,
    sample,
    pop_label,
    pop_values,
    cell_names,
    ref_seq) {
  stopifnot(requireNamespace("reticulate", quietly = TRUE))
  ref_dict <- paste0("{", paste0("'", names(ref_seq), "': ", ref_seq, collapse = ", "), "}")
  reticulate::py_run_string(sprintf(
    "
import numpy as np
from scipy.sparse import coo_matrix
import pandas as pd
from anndata import AnnData

rows = np.array([%s], dtype=np.int32)
cols = np.array([%s], dtype=np.int32)
data = np.array([%s], dtype=np.float64)
fragment_paired = coo_matrix((data, (rows, cols)), shape=(%d, %d)).tocsr()

obs = pd.DataFrame({
    'Sample': [%s],
    '%s': [%s],
}, index=[%s])

uns = {'reference_sequences': %s}

py_test_adata = AnnData(
    obs=obs,
    obsm={'fragment_paired': fragment_paired},
    uns=uns,
)
",
    paste(rows, collapse = ", "),
    paste(cols, collapse = ", "),
    paste(data, collapse = ", "),
    as.integer(n_obs),
    as.integer(n_cols),
    paste0("'", sample, "'", collapse = ", "),
    pop_label,
    paste0("'", pop_values, "'", collapse = ", "),
    paste0("'", cell_names, "'", collapse = ", "),
    ref_dict
  ))
  reticulate::py$py_test_adata
}

test_that("integration: reticulate AnnData round-trip when Python available", {
  skip_if_not(requireNamespace("reticulate", quietly = TRUE))
  skip_if_not(reticulate::py_module_available("anndata"), "anndata not available")
  skip_if_not(reticulate::py_module_available("numpy"), "numpy not available")
  skip_if_not(reticulate::py_module_available("scipy"), "scipy not available")

  py_adata <- .make_py_snap_adata(
    rows = c(0L, 0L, 1L),
    cols = c(0L, 10L, 500L),
    data = c(50, 60, 100),
    n_obs = 2L,
    n_cols = 2000L,
    sample = c("s1", "s1"),
    pop_label = "cluster",
    pop_values = c("A", "A"),
    cell_names = c("c1", "c2"),
    ref_seq = c(chr1 = 1000L, chr2 = 1000L)
  )

  bundle <- MOCHA::import_snap_atac(
    py_adata,
    cellPopLabel = "cluster",
    sampleColumn = "Sample",
    verbose = FALSE
  )

  expect_s4_class(bundle$ATACFragments, "GRangesList")
  expect_equal(names(bundle$ATACFragments), "s1")
  expect_equal(nrow(bundle$cellColData), 2L)
})

test_that("integration: import_snap_atac output works with callOpenTiles sample-level path", {
  skip_if_not(requireNamespace("reticulate", quietly = TRUE))
  skip_if_not(reticulate::py_module_available("anndata"))
  skip_if_not(reticulate::py_module_available("numpy"))
  skip_if_not(reticulate::py_module_available("scipy"))
  skip_if_not(
    requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE) &&
      requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE) &&
      requireNamespace("org.Hs.eg.db", quietly = TRUE)
  )

  py_adata <- .make_py_snap_adata(
    rows = c(0L, 1L),
    cols = c(760100L, 760200L),
    data = c(50, 50),
    n_obs = 2L,
    n_cols = 249250621L,
    sample = c("sample1", "sample2"),
    pop_label = "cellPop",
    pop_values = c("t_cd8_temra", "t_cd8_temra"),
    cell_names = c("c1", "c2"),
    ref_seq = c(chr1 = 249250621L)
  )

  bundle <- MOCHA::import_snap_atac(
    py_adata,
    cellPopLabel = "cellPop",
    sampleColumn = "Sample",
    chromPrefix = "none",
    verbose = FALSE
  )

  expect_warning(
    tiles <- MOCHA::callOpenTiles(
      ATACFragments = bundle$ATACFragments,
      cellColData = bundle$cellColData,
      blackList = MOCHA::exampleBlackList,
      genome = "hg19",
      TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
      OrgDb = "org.Hs.eg.db",
      outDir = tempdir(),
      cellPopLabel = "cellPop",
      cellPopulations = "t_cd8_temra",
      cellCol = bundle$cellCol,
      studySignal = 10,
      numCores = 1
    ),
    regexp = NA
  )
  expect_s4_class(tiles, "MultiAssayExperiment")
})
