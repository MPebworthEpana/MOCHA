test_that("coverage layout detection distinguishes structured, legacy, and missing", {
  legacy_dir <- make_synthetic_coverage_dir("C2", "Sample1")
  tracks_dir <- make_synthetic_tracks_dir("C3", "Sample1")

  expect_equal(MOCHA:::.detect_coverage_layout(legacy_dir, "C2"), "legacy")
  expect_equal(MOCHA:::.detect_coverage_layout(legacy_dir, "C9"), "missing")
  expect_equal(MOCHA:::.detect_coverage_layout(tracks_dir, "C3"), "structured")
  expect_true(MOCHA:::.coverage_tracks_exist(tracks_dir, "C3"))
  expect_false(MOCHA:::.coverage_tracks_exist(tracks_dir, "C9"))
})

test_that("sanitize_track_filename uses reversible percent encoding", {
  expect_equal(
    MOCHA:::.sanitize_track_filename("Sample A"),
    "Sample%20A"
  )
  expect_equal(
    MOCHA:::.sanitize_track_filename("scATAC.BMMC.R1"),
    "scATAC%2EBMMC%2ER1"
  )
  expect_equal(
    MOCHA:::.unsanitize_track_filename("Sample%20A"),
    "Sample A"
  )
  expect_false(identical(
    MOCHA:::.sanitize_track_filename("Sample A"),
    MOCHA:::.sanitize_track_filename("Sample.A")
  ))
})

test_that("validate_track_sample_ids errors on duplicate sample IDs", {
  expect_error(
    MOCHA:::.validate_track_sample_ids(c("Sample1", "Sample1")),
    "collide after filename encoding"
  )
})

test_that("coverage_tracks_exist requires complete sample sets when manifest is present", {
  tracks_dir <- make_synthetic_tracks_dir("C2", c("Sample1", "Sample2"))
  expect_true(MOCHA:::.coverage_tracks_exist(
    tracks_dir,
    "C2",
    sampleIds = c("Sample1", "Sample2")
  ))
  expect_false(MOCHA:::.coverage_tracks_exist(
    tracks_dir,
    "C2",
    sampleIds = c("Sample1", "Sample2", "Sample3")
  ))

  acc_file <- MOCHA:::.coverage_track_file_path(tracks_dir, "C2", "Accessibility", "Sample2")
  file.remove(acc_file)
  expect_false(MOCHA:::.coverage_tracks_exist(
    tracks_dir,
    "C2",
    sampleIds = c("Sample1", "Sample2")
  ))
})

test_that("mixed legacy and structured layouts warn and prefer structured", {
  legacy_dir <- make_synthetic_coverage_dir("C2", "Sample1")
  make_synthetic_tracks_dir("C2", "Sample2", dir = legacy_dir)

  expect_warning(
    bundle <- MOCHA:::.readCoverageBundle(legacy_dir, "C2"),
    "Both structured tracks and legacy RDS coverage exist"
  )
  expect_equal(names(bundle$Accessibility), "Sample2")
})

test_that("migrateCoverageToTracks works on tileResults MultiAssayExperiment", {
  skip_on_cran()

  legacy_dir <- make_synthetic_coverage_dir("C3", "scATAC_BMMC_R1")
  tiles <- MOCHA:::testTileResultsMultisample
  tiles@metadata$Directory <- legacy_dir

  migrated <- MOCHA::migrateCoverageToTracks(
    tiles,
    cellPopulations = "C3",
    removeLegacy = TRUE,
    force = TRUE,
    verbose = FALSE
  )

  expect_equal(migrated@metadata$CoverageLayout, "tracks")
  expect_s4_class(migrated, "MultiAssayExperiment")
  expect_true(MOCHA:::.coverage_tracks_exist(legacy_dir, "C3", sampleIds = "scATAC_BMMC_R1"))
})

test_that("packMOCHA archives structured tracks directory", {
  skip_on_cran()
  skip_if_not_installed("zip")

  tracks_dir <- make_synthetic_tracks_dir("C3", "scATAC_BMMC_R1")
  tiles <- MOCHA:::testTileResultsMultisample
  tiles@metadata$Directory <- tracks_dir

  zip_path <- file.path(tempdir(), "mocha_pack_tracks_test.zip")
  on.exit(unlink(c(zip_path, paste0(zip_path, "*")), recursive = TRUE), add = TRUE)

  MOCHA::packMOCHA(tiles, zipfile = zip_path)
  listing <- zip::zip_list(zip_path)
  zip_files <- if ("filename" %in% names(listing)) {
    listing$filename
  } else {
    listing[[1L]]
  }
  expect_true(any(grepl("manifest\\.rds", zip_files)))
  expect_true(any(grepl("Accessibility/", zip_files)))

  unpack_dir <- file.path(tempdir(), "mocha_unpack_tracks_test")
  on.exit(unlink(unpack_dir, recursive = TRUE), add = TRUE)
  unpacked <- MOCHA::unpackMOCHA(zip_path, unpack_dir)
  unpacked_dir <- unpacked@metadata$Directory
  expect_true(MOCHA:::.coverage_tracks_exist(unpacked_dir, "C3", sampleIds = "scATAC_BMMC_R1"))
})

test_that("structured tracks round-trip through write and load helpers", {
  tracks_dir <- make_synthetic_tracks_dir(
    cell_populations = "C2",
    sample_ids = c("Sample A", "scATAC.BMMC.R1")
  )

  manifest <- readRDS(file.path(tracks_dir, "tracks", "manifest.rds"))
  expect_equal(manifest$format, "bigwig")
  expect_true("C2" %in% manifest$cellPopulations)
  expect_equal(manifest$sample_name_map[["Sample%20A"]], "Sample A")

  loaded <- MOCHA:::.loadCoverageTracks(
    outDir = tracks_dir,
    cellPop = "C2",
    trackType = "accessibility"
  )
  expect_equal(sort(names(loaded)), c("Sample A", "scATAC.BMMC.R1"))
  expect_true(all(vapply(loaded, function(x) {
    sum(S4Vectors::mcols(x)$score) > 0
  }, logical(1L))))
})

test_that("readCoverageBundle returns legacy and structured bundles", {
  legacy_dir <- make_synthetic_coverage_dir("C2", "Sample1")
  tracks_dir <- make_synthetic_tracks_dir("C3", "Sample1")

  legacy_bundle <- MOCHA:::.readCoverageBundle(legacy_dir, "C2")
  expect_equal(names(legacy_bundle), c("Accessibility", "Insertions"))

  structured_bundle <- MOCHA:::.readCoverageBundle(tracks_dir, "C3")
  expect_equal(names(structured_bundle), c("Accessibility", "Insertions"))
  expect_equal(names(structured_bundle$Accessibility), "Sample1")
})

test_that("loadCoverageTracks supports region-aware structured reads", {
  tracks_dir <- make_synthetic_tracks_dir("C2", "Sample1", region = "chr1:100-500")
  region <- MOCHA::StringsToGRanges("chr1:150-250")

  loaded <- MOCHA:::.loadCoverageTracks(
    outDir = tracks_dir,
    cellPop = "C2",
    trackType = "accessibility",
    region = region
  )

  expect_length(loaded, 1L)
  expect_true(all(GenomicRanges::start(loaded[[1]]) >= 150))
  expect_true(all(GenomicRanges::end(loaded[[1]]) <= 250))
})

test_that("migrateCoverageToTracks converts legacy bundles to structured tracks", {
  skip_on_cran()
  skip_if_not_installed("SummarizedExperiment")

  legacy_dir <- make_synthetic_coverage_dir("C3", c("scATAC_BMMC_R1", "scATAC_CD34_BMMC_R1"))
  capture.output(
    stm <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResultsMultisample,
      cellPopulations = "all",
      reproducibilityThreshold = 0,
      numCores = 1
    )
  )
  stm@metadata$Directory <- legacy_dir

  migrated <- MOCHA::migrateCoverageToTracks(
    stm,
    cellPopulations = "C3",
    removeLegacy = TRUE,
    force = TRUE,
    verbose = FALSE
  )

  expect_equal(migrated@metadata$CoverageLayout, "tracks")
  expect_false(file.exists(file.path(legacy_dir, "C3_CoverageFiles.RDS")))
  expect_true(MOCHA:::.coverage_tracks_exist(
    legacy_dir,
    "C3",
    sampleIds = c("scATAC_BMMC_R1", "scATAC_CD34_BMMC_R1")
  ))
})

test_that("extractRegion works with structured and legacy coverage layouts", {
  skip_on_cran()
  skip_if_not_installed("SummarizedExperiment")

  region <- "chr1:18137866-18237865"
  sample_ids <- c("scATAC_BMMC_R1", "scATAC_CD34_BMMC_R1")

  capture.output(
    stm <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResultsMultisample,
      cellPopulations = "all",
      reproducibilityThreshold = 0,
      numCores = 1
    )
  )
  stm <- stm[, stm$Sample %in% sample_ids]

  legacy_dir <- make_synthetic_coverage_dir("C3", sample_ids, region = region)
  tracks_dir <- make_synthetic_tracks_dir("C3", sample_ids, region = region)

  stm_legacy <- stm
  stm_legacy@metadata$Directory <- legacy_dir
  legacy_se <- MOCHA::extractRegion(
    SampleTileObj = stm_legacy,
    cellPopulations = "C3",
    region = region,
    numCores = 1,
    sampleSpecific = FALSE
  )

  stm_tracks <- stm
  stm_tracks@metadata$Directory <- tracks_dir
  tracks_se <- MOCHA::extractRegion(
    SampleTileObj = stm_tracks,
    cellPopulations = "C3",
    region = region,
    numCores = 1,
    sampleSpecific = FALSE
  )

  expect_equal(dim(legacy_se), dim(tracks_se))
  expect_equal(
    SummarizedExperiment::assayNames(legacy_se),
    SummarizedExperiment::assayNames(tracks_se)
  )
})

test_that("validate_track_path_label rejects unsafe characters", {
  expect_error(
    MOCHA:::.validate_track_path_label("Sample/1", labelType = "sample"),
    "Invalid sample label"
  )
  expect_error(
    MOCHA:::.validate_track_path_label("C2/Bad", labelType = "cellPop"),
    "Invalid cellPop label"
  )
})

test_that("writeCoverageTracks removes orphaned bigWigs when sample set shrinks", {
  tracks_dir <- make_synthetic_tracks_dir("C2", c("Sample1", "Sample2"))
  region_gr <- MOCHA::StringsToGRanges("chr1:100-500")
  S4Vectors::mcols(region_gr)$score <- 1L
  bundle <- list(
    Accessibility = stats::setNames(list(region_gr), "Sample1"),
    Insertions = stats::setNames(list(region_gr), "Sample1")
  )

  MOCHA:::.writeCoverageTracks(
    covFiles = bundle,
    outDir = tracks_dir,
    cellPop = "C2",
    force = TRUE,
    verbose = FALSE
  )

  expect_false(file.exists(
    MOCHA:::.coverage_track_file_path(tracks_dir, "C2", "Accessibility", "Sample2")
  ))
  expect_false(file.exists(
    MOCHA:::.coverage_track_file_path(tracks_dir, "C2", "Insertions", "Sample2")
  ))
  manifest <- readRDS(file.path(tracks_dir, "tracks", "manifest.rds"))
  expect_equal(manifest$cellPopulationSamples$C2, "Sample1")
})

test_that("migrateCoverageToTracks repairs incomplete structured tracks from legacy", {
  skip_on_cran()

  legacy_dir <- make_synthetic_coverage_dir("C3", c("Sample1", "Sample2"))
  make_synthetic_tracks_dir("C3", c("Sample1"), dir = legacy_dir)

  tiles <- MOCHA:::testTileResultsMultisample
  tiles@metadata$Directory <- legacy_dir

  migrated <- MOCHA::migrateCoverageToTracks(
    tiles,
    cellPopulations = "C3",
    removeLegacy = FALSE,
    force = FALSE,
    verbose = FALSE
  )

  expect_equal(migrated@metadata$CoverageLayout, "tracks")
  expect_true(MOCHA:::.coverage_tracks_exist(
    legacy_dir,
    "C3",
    sampleIds = c("Sample1", "Sample2")
  ))
})

test_that("migrateCoverageToTracks warns when structured tracks incomplete without legacy", {
  skip_on_cran()

  tracks_dir <- make_synthetic_tracks_dir("C3", c("Sample1", "Sample2"))
  acc_file <- MOCHA:::.coverage_track_file_path(tracks_dir, "C3", "Accessibility", "Sample2")
  file.remove(acc_file)

  tiles <- MOCHA:::testTileResultsMultisample
  tiles@metadata$Directory <- tracks_dir
  tiles@metadata$CoverageLayout <- NULL

  expect_warning(
    out <- MOCHA::migrateCoverageToTracks(
      tiles,
      cellPopulations = "C3",
      verbose = FALSE
    ),
    "incomplete and no legacy RDS bundle"
  )
  expect_null(out@metadata$CoverageLayout)
})

test_that("loadCoverageTracks caches full-genome structured reads", {
  MOCHA:::.clear_coverage_cache()
  on.exit(MOCHA:::.clear_coverage_cache(), add = TRUE)

  tracks_dir <- make_synthetic_tracks_dir("C2", "Sample1")
  MOCHA:::.loadCoverageTracks(
    outDir = tracks_dir,
    cellPop = "C2",
    trackType = "accessibility"
  )
  cache_keys <- ls(envir = MOCHA:::.mocha_coverage_cache)
  expect_length(cache_keys, 1L)

  bundle <- MOCHA:::.readCoverageBundle(tracks_dir, "C2")
  expect_equal(names(bundle), c("Accessibility", "Insertions"))
  expect_length(ls(envir = MOCHA:::.mocha_coverage_cache), 2L)
})

test_that("exportCoverage copies structured sample-specific tracks when available", {
  skip_on_cran()

  tracks_dir <- make_synthetic_tracks_dir("C3", "scATAC_BMMC_R1")
  capture.output(
    stm <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResultsMultisample,
      cellPopulations = "all",
      reproducibilityThreshold = 0,
      numCores = 1
    )
  )
  stm@metadata$Directory <- tracks_dir
  stm <- stm[, stm$Sample == "scATAC_BMMC_R1"]

  export_dir <- tempfile(pattern = "mocha_export_cov_")
  on.exit(unlink(export_dir, recursive = TRUE), add = TRUE)

  expect_message(
    out <- MOCHA::exportCoverage(
      SampleTileObject = stm,
      dir = export_dir,
      cellPopulations = "C3",
      sampleSpecific = TRUE,
      saveFile = TRUE,
      numCores = 1,
      verbose = TRUE
    ),
    "Using existing structured coverage tracks"
  )

  expected_file <- file.path(export_dir, "C3__scATAC_BMMC_R1_Coverage.bw")
  expect_true(file.exists(expected_file))
  expect_true(length(out) >= 1L)
})
