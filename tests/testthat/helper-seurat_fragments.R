#' Build a minimal Seurat + Signac object for MOCHA tests.
#' @noRd
make_test_seurat_chromatin <- function(
    fragment_lines,
    meta,
    assay_name = "ATAC",
    genome = "hg19") {
  if (!requireNamespace("Seurat", quietly = TRUE) ||
      !requireNamespace("Signac", quietly = TRUE)) {
    return(NULL)
  }

  frag_path <- tempfile(fileext = ".tsv")
  writeLines(fragment_lines, con = frag_path)
  on.exit(unlink(frag_path), add = TRUE)

  peak_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1, end = 500),
    strand = "*"
  )
  counts <- Matrix::Matrix(
    1,
    nrow = 1,
    ncol = nrow(meta),
    sparse = TRUE
  )
  rownames(counts) <- "peak1"
  colnames(counts) <- rownames(meta)

  chrom_assay <- Signac::CreateChromatinAssay(
    counts = counts,
    ranges = peak_gr,
    fragments = frag_path,
    min.cells = 0,
    min.features = 0
  )

  obj <- Seurat::CreateSeuratObject(
    counts = chrom_assay,
    assay = assay_name,
    min.cells = 0,
    min.features = 0
  )
  obj@meta.data <- cbind(obj@meta.data, meta)
  rownames(obj@meta.data) <- rownames(meta)

  obj
}
