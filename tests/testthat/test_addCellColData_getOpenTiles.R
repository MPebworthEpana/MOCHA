test_that("addCellColData attaches a sample-aligned column on a SampleTileMatrix", {
  stm <- make_synthetic_sample_tile_matrix()
  samples <- rownames(SummarizedExperiment::colData(stm))

  out <- MOCHA::addCellColData(
    stm,
    name = "BatchID",
    value = paste0("batch_", seq_along(samples))
  )
  cd <- SummarizedExperiment::colData(out)
  expect_true("BatchID" %in% colnames(cd))
  expect_identical(cd$BatchID, paste0("batch_", seq_along(samples)))
})

test_that("addCellColData supports the `samples =` mapping path", {
  stm <- make_synthetic_sample_tile_matrix()
  samples <- rownames(SummarizedExperiment::colData(stm))
  subset_samples <- samples[c(1L, 3L)]

  out <- MOCHA::addCellColData(
    stm,
    name = "Subset",
    value = c("a", "b"),
    samples = subset_samples
  )
  cd <- SummarizedExperiment::colData(out)
  expect_equal(cd$Subset[match(subset_samples, samples)], c("a", "b"))
  expect_true(all(is.na(cd$Subset[-match(subset_samples, samples)])))
})

test_that("addCellColData refuses to overwrite without force = TRUE", {
  stm <- make_synthetic_sample_tile_matrix()
  expect_error(
    MOCHA::addCellColData(stm, name = "GroupA", value = rep("X", ncol(stm))),
    "already exists"
  )
  out <- MOCHA::addCellColData(
    stm,
    name = "GroupA",
    value = rep("X", ncol(stm)),
    force = TRUE
  )
  expect_true(all(SummarizedExperiment::colData(out)$GroupA == "X"))
})

make_synthetic_tile_results <- function() {
  build_sample_gr <- function(peaks, sample_offset = 0) {
    n <- length(peaks)
    gr <- GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges = IRanges::IRanges(
        start = sample_offset + seq(1, by = 100, length.out = n),
        end = sample_offset + seq(50, by = 100, length.out = n)
      )
    )
    S4Vectors::mcols(gr)$peak <- peaks
    gr
  }
  build_re <- function(peak_list) {
    grl <- lapply(peak_list, build_sample_gr)
    names(grl) <- paste0("sample", seq_along(grl))
    RaggedExperiment::RaggedExperiment(grl)
  }
  re_c2 <- build_re(list(c(TRUE, FALSE, TRUE), c(TRUE, FALSE, FALSE)))
  re_c5 <- build_re(list(c(FALSE, TRUE, TRUE), c(FALSE, TRUE, FALSE)))
  colData_df <- S4Vectors::DataFrame(
    Sample = c("sample1", "sample2"),
    row.names = c("sample1", "sample2")
  )
  MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(C2 = re_c2, C5 = re_c5),
    colData = colData_df
  )
}

test_that("getOpenTiles returns a GRangesList keyed by cell population", {
  mae <- make_synthetic_tile_results()
  grl <- MOCHA::getOpenTiles(mae)
  expect_s4_class(grl, "GRangesList")
  expect_setequal(names(grl), c("C2", "C5"))
  expect_gt(length(grl[["C2"]]), 0L)
  expect_gt(length(grl[["C5"]]), 0L)

  df <- MOCHA::getOpenTiles(mae, returnType = "data.frame")
  expect_s3_class(df, "data.frame")
  expect_true(all(c("CellPopulation", "seqnames", "start", "end") %in% colnames(df)))
  expect_setequal(unique(df$CellPopulation), c("C2", "C5"))
})

test_that("getOpenTiles errors on unknown cell populations", {
  mae <- make_synthetic_tile_results()
  expect_error(MOCHA::getOpenTiles(mae, cellPopulations = "Nope"), "not found")
})
