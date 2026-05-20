test_that("estimateDropoutModel and assessDropout run on simulated dropout data", {
  stm <- make_synthetic_dropout_tsam()
  model <- MOCHA::estimateDropoutModel(stm, cellPopulation = "C2", verbose = FALSE)
  expect_s3_class(model, "DropoutModel")
  expect_true("C2" %in% names(model$models))
  expect_false(is.na(model$diagnostics$C2$auroc))
  expect_equal(model$diagnostics$C2$auroc_method, "leave_one_sample_out_train_only_features")
  expect_gte(model$diagnostics$C2$auroc, 0.5)

  stm2 <- MOCHA::assessDropout(stm, model = model)
  expect_true("DropoutProb_C2" %in% SummarizedExperiment::assayNames(stm2))
  expect_false(is.null(S4Vectors::metadata(stm2)$dropoutModels))

  probs <- MOCHA::getDropoutProb(stm2, "C2")
  expect_equal(dim(probs), dim(SummarizedExperiment::assay(stm2, "C2")))

  cls <- MOCHA::classifyZeros(stm2, "C2")
  tech_tiles <- rownames(probs)[seq(21, 40)]
  bio_tiles <- rownames(probs)[seq_len(20)]
  group_a <- paste0("sample_", 1:5)
  tech_zero_probs <- probs[tech_tiles, group_a][probs[tech_tiles, group_a] > 0]
  bio_zero_probs <- probs[bio_tiles, group_a][probs[bio_tiles, group_a] > 0]
  expect_gt(length(tech_zero_probs), 0)
  expect_gt(length(bio_zero_probs), 0)
  expect_lt(mean(tech_zero_probs), mean(bio_zero_probs))

  tech_zeros_a <- cls[tech_tiles, group_a]
  expect_gt(mean(tech_zeros_a == "technical", na.rm = TRUE), 0.5)
})

test_that("plotDropoutDiagnostics calibration requires dropoutModels", {
  stm <- make_synthetic_dropout_tsam()
  expect_error(
    MOCHA::plotDropoutDiagnostics(stm, "C2", type = "calibration"),
    "dropoutModels"
  )
})

test_that("plotDropoutDiagnostics returns ggplot objects", {
  stm <- MOCHA::assessDropout(make_synthetic_dropout_tsam())
  p1 <- MOCHA::plotDropoutDiagnostics(stm, "C2", type = "probHist")
  expect_s3_class(p1, "ggplot")
  p2 <- MOCHA::plotDropoutDiagnostics(stm, "C2", type = "calibration")
  expect_s3_class(p2, "ggplot")
})

test_that("dropoutAdjustment none ignores DropoutProb assays for filtering", {
  stm <- MOCHA::assessDropout(make_synthetic_dropout_tsam())

  none_adj <- suppressWarnings(
    MOCHA::getDifferentialAccessibleTiles(
      stm,
      cellPopulations = "C2",
      groupColumn = "GroupA",
      foreground = "A",
      background = "B",
      signalThreshold = 8,
      qValueMethod = "experimental",
      dropoutAdjustment = "none",
      outputGRanges = FALSE,
      numCores = 1,
      verbose = FALSE
    )
  )

  stm_no_dropout <- make_synthetic_dropout_tsam()
  no_assay <- suppressWarnings(
    MOCHA::getDifferentialAccessibleTiles(
      stm_no_dropout,
      cellPopulations = "C2",
      groupColumn = "GroupA",
      foreground = "A",
      background = "B",
      signalThreshold = 8,
      qValueMethod = "experimental",
      dropoutAdjustment = "none",
      outputGRanges = FALSE,
      numCores = 1,
      verbose = FALSE
    )
  )

  expect_equal(none_adj$Tile, no_assay$Tile)
  expect_equal(none_adj$P_value, no_assay$P_value)
})

test_that(".ensure_dropout_assay fails when model cannot be fitted", {
  stm <- make_synthetic_dropout_tsam(n_bio_closed = 1L, n_tech_dropout = 1L, n_per_group = 2L)
  expect_error(
    MOCHA:::.ensure_dropout_assay(stm, "C2", verbose = FALSE),
    "DropoutProb|Insufficient data"
  )
})

test_that(".ensure_dropout_assay returns object when assay is valid", {
  stm <- MOCHA::assessDropout(make_synthetic_dropout_tsam())
  out <- MOCHA:::.ensure_dropout_assay(stm, "C2", verbose = FALSE)
  expect_true(MOCHA:::.has_valid_dropout_assay(out, "C2"))
})

test_that("dropoutAdjustment none matches legacy NA-as-zero zero rates", {
  stm <- make_synthetic_dropout_tsam(n_per_group = 5L)
  raw_mat <- MOCHA::getCellPopMatrix(stm, "C2", dropSamples = FALSE, NAtoZero = FALSE)
  group <- as.numeric(colnames(raw_mat) %in% paste0("sample_", 1:5))
  matFilled <- raw_mat
  matFilled[is.na(matFilled)] <- 0
  legacy_A <- rowMeans(matFilled[, group == 1, drop = FALSE] == 0)
  adj <- MOCHA:::.compute_adjusted_zero_rates(
    raw_mat, NULL, group, "none"
  )
  expect_equal(adj$zero_A, legacy_A)
})

test_that("dropoutAdjustment biological_only changes zero-driven tile retention", {
  stm <- make_synthetic_dropout_tsam(n_per_group = 5L)
  stm <- MOCHA::assessDropout(stm)

  raw_mat <- MOCHA::getCellPopMatrix(stm, "C2", dropSamples = FALSE, NAtoZero = FALSE)
  group <- as.numeric(colnames(raw_mat) %in% paste0("sample_", 1:5))
  prob <- MOCHA::getDropoutProb(stm, "C2")
  prob <- prob[rownames(raw_mat), colnames(raw_mat), drop = FALSE]

  raw_rates <- MOCHA:::.compute_adjusted_zero_rates(
    raw_mat, NULL, group, "none"
  )
  adj_rates <- MOCHA:::.compute_adjusted_zero_rates(
    raw_mat, prob, group, "biological_only", bioThreshold = 0.8
  )

  expect_lt(mean(adj_rates$zero_A), mean(raw_rates$zero_A))
})

test_that("DropoutProb assays survive packMOCHA and unpackMOCHA round-trip", {
  skip_if_not(requireNamespace("zip", quietly = TRUE))

  stm <- MOCHA::assessDropout(make_synthetic_dropout_tsam())
  probs_before <- MOCHA::getDropoutProb(stm, "C2")
  model_before <- S4Vectors::metadata(stm)$dropoutModels

  mdir <- tempfile("mocha_pack_dropout_")
  dir.create(mdir, recursive = TRUE)
  md <- S4Vectors::metadata(stm)
  md$Directory <- mdir
  S4Vectors::metadata(stm) <- md

  zip_path <- tempfile(fileext = ".zip")
  unpack_dir <- tempfile("mocha_unpack_dropout_")
  dir.create(unpack_dir, recursive = TRUE)

  MOCHA::packMOCHA(stm, zipfile = zip_path)
  stm_unpacked <- MOCHA::unpackMOCHA(zip_path, unpack_dir)

  probs_after <- MOCHA::getDropoutProb(stm_unpacked, "C2")
  expect_equal(probs_after, probs_before)
  expect_s3_class(S4Vectors::metadata(stm_unpacked)$dropoutModels, "DropoutModel")
  expect_equal(
    names(S4Vectors::metadata(stm_unpacked)$dropoutModels$models),
    names(model_before$models)
  )
})

test_that("DropoutProb assays and metadata survive subset and serialization", {
  stm <- MOCHA::assessDropout(make_synthetic_dropout_tsam())

  sub <- MOCHA::subsetMOCHAObject(
    stm,
    subsetBy = "celltype",
    groupList = "C2",
    verbose = FALSE
  )
  expect_true("DropoutProb_C2" %in% SummarizedExperiment::assayNames(sub))
  expect_false(is.null(S4Vectors::metadata(sub)$dropoutModels))

  sub_samples <- MOCHA::subsetMOCHAObject(
    stm,
    subsetBy = "GroupA",
    groupList = "A",
    verbose = FALSE
  )
  expect_true("DropoutProb_C2" %in% SummarizedExperiment::assayNames(sub_samples))

  tmp <- tempfile(fileext = ".rds")
  saveRDS(stm, tmp)
  loaded <- readRDS(tmp)
  unlink(tmp)
  expect_true("DropoutProb_C2" %in% SummarizedExperiment::assayNames(loaded))
  expect_s3_class(S4Vectors::metadata(loaded)$dropoutModels, "DropoutModel")
})
