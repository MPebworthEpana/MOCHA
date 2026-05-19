skip_on_cran()

test_that(".readFragmentFile filters barcodes", {
  frag_path <- tempfile(fileext = ".tsv")
  writeLines(
    c(
      "chr1\t760101\t760110\tAAAC\t1",
      "chr1\t760111\t760120\tBBBB\t1"
    ),
    con = frag_path
  )
  on.exit(unlink(frag_path), add = TRUE)

  df <- MOCHA:::.readFragmentFile(frag_path, barcodes = "AAAC")
  expect_equal(nrow(df), 1L)
  expect_equal(df$barcode[[1]], "AAAC")
})

test_that("seuratToMOCHAInputs requires Seurat and Signac", {
  skip_if(requireNamespace("Seurat", quietly = TRUE) &&
    requireNamespace("Signac", quietly = TRUE),
  message = "Seurat/Signac installed; skip negative dependency test"
  )

  expect_error(
    MOCHA:::seuratToMOCHAInputs(
      seuratObj = structure(list(), class = "Seurat"),
      cellPopLabel = "cellPop"
    ),
    regexp = "required for Seurat/Signac"
  )
})

if (requireNamespace("Seurat", quietly = TRUE) &&
  requireNamespace("Signac", quietly = TRUE)) {
  test_that("seuratToMOCHAInputs converts fragments via Fragment@cells", {
    meta <- data.frame(
      Sample = c("sample1", "sample1", "sample2"),
      cellPop = c("TypeA", "TypeA", "TypeA"),
      row.names = c("c1", "c2", "c3"),
      stringsAsFactors = FALSE
    )

    lines <- c(
      "chr1\t760101\t760110\tAAAC\t1",
      "chr1\t760111\t760120\tAAAG\t1",
      "chr1\t760121\t760130\tTTTC\t1"
    )

    obj <- make_test_seurat_chromatin(
      fragment_lines = lines,
      meta = meta
    )
    skip_if(is.null(obj), "Could not build test Seurat object (indexed fragments required)")

    parsed <- MOCHA::seuratToMOCHAInputs(
      seuratObj = obj,
      cellPopLabel = "cellPop",
      verbose = FALSE
    )

    expect_type(parsed, "list")
    expect_s4_class(parsed$ATACFragments, "GRangesList")
    expect_equal(names(parsed$ATACFragments), c("sample1", "sample2"))
    expect_true("RG" %in% colnames(GenomicRanges::mcols(parsed$ATACFragments[[1]])))
    expect_true(all(c("c1", "c2") %in% unique(
      GenomicRanges::mcols(parsed$ATACFragments[["sample1"]])$RG
    )))
  })

  test_that("seuratToMOCHAInputs errors on missing metadata columns", {
    meta <- data.frame(
      Sample = rep("sample1", 2),
      row.names = c("c1", "c2"),
      stringsAsFactors = FALSE
    )
    obj <- make_test_seurat_chromatin(
      fragment_lines = c(
        "chr1\t760101\t760110\tAAAC\t1",
        "chr1\t760111\t760120\tAAAG\t1"
      ),
      meta = meta
    )
    skip_if(is.null(obj), "Could not build test Seurat object (indexed fragments required)")

    expect_error(
      MOCHA::seuratToMOCHAInputs(obj, cellPopLabel = "cellPop"),
      regexp = "cellPopLabel"
    )
  })

  test_that("seuratToMOCHAInputs errors on invalid assay", {
    meta <- data.frame(
      Sample = rep("sample1", 2),
      cellPop = rep("TypeA", 2),
      row.names = c("c1", "c2"),
      stringsAsFactors = FALSE
    )
    obj <- make_test_seurat_chromatin(
      fragment_lines = c(
        "chr1\t760101\t760110\tAAAC\t1",
        "chr1\t760111\t760120\tAAAG\t1"
      ),
      meta = meta
    )
    skip_if(is.null(obj), "Could not build test Seurat object (indexed fragments required)")

    rna_counts <- Matrix::Matrix(1, nrow = 1, ncol = 2, sparse = TRUE)
    rownames(rna_counts) <- "gene1"
    colnames(rna_counts) <- c("c1", "c2")
    obj[["RNA"]] <- Seurat::CreateAssayObject(counts = rna_counts)
    expect_error(
      MOCHA::seuratToMOCHAInputs(obj, assay = "RNA", cellPopLabel = "cellPop"),
      regexp = "not a ChromatinAssay"
    )
  })

  test_that("seuratToMOCHAInputs errors on barcode mismatch", {
    meta <- data.frame(
      Sample = rep("sample1", 2),
      cellPop = rep("TypeA", 2),
      row.names = c("c1", "c2"),
      stringsAsFactors = FALSE
    )
    obj <- make_test_seurat_chromatin(
      fragment_lines = c(
        "chr1\t760101\t760110\tAAAC\t1",
        "chr1\t760111\t760120\tAAAG\t1"
      ),
      meta = meta
    )
    skip_if(is.null(obj), "Could not build test Seurat object (indexed fragments required)")

    rownames(obj@meta.data) <- c("wrong_cell", "wrong_cell2")

    expect_error(
      MOCHA::seuratToMOCHAInputs(obj, cellPopLabel = "cellPop"),
      regexp = "Barcode prefix mismatch"
    )
  })

  test_that("fragmentPathColumn mode resolves per-sample paths", {
    frag1 <- tempfile(fileext = ".tsv")
    frag2 <- tempfile(fileext = ".tsv")
    writeLines("chr1\t760101\t760110\tAAAC\t1", frag1)
    writeLines("chr1\t760111\t760120\tAAAG\t1", frag2)

    meta <- data.frame(
      Sample = c("sample1", "sample2"),
      cellPop = c("TypeA", "TypeA"),
      frag_path = c(frag1, frag2),
      row.names = c("AAAC", "AAAG"),
      stringsAsFactors = FALSE
    )

    obj <- make_test_seurat_chromatin(
      fragment_lines = c(
        "chr1\t760101\t760110\tAAAC\t1",
        "chr1\t760111\t760120\tAAAG\t1"
      ),
      meta = meta
    )
    skip_if(is.null(obj), "Could not build test Seurat object (indexed fragments required)")

    parsed <- MOCHA::seuratToMOCHAInputs(
      seuratObj = obj,
      cellPopLabel = "cellPop",
      fragmentPathColumn = "frag_path",
      verbose = FALSE
    )

    expect_equal(length(parsed$ATACFragments), 2L)
    expect_gt(length(parsed$ATACFragments[["sample1"]]), 0L)

    unlink(c(frag1, frag2))
  })

  test_that("callOpenTiles Seurat method is registered when Seurat is installed", {
    meth <- tryCatch(
      methods::getMethod("callOpenTiles", signature = c(ATACFragments = "Seurat")),
      error = function(e) NULL
    )
    expect_false(is.null(meth))
  })
}
