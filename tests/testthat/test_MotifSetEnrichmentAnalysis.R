test_that("MotifSetEnrichmentAnalysis runs on bundled fixtures", {
  motifEnrichmentDF <- read.csv(mocha_fixture_path("msea_motif_enrichment.csv"))
  ligandTFMatrix <- matrix(
    c(
      1, 1, 0,
      1, 0, 1,
      0, 1, 1,
      1, 0, 0,
      0, 0, 0
    ),
    nrow = 5,
    ncol = 3,
    byrow = TRUE,
    dimnames = list(
      c("SPI1", "GATA1", "RUNX1", "CEBPA", "MEF2C"),
      c("TGFB1", "VEGFA", "IL6")
    )
  )

  results <- MOCHA::MotifSetEnrichmentAnalysis(
    ligandTFMatrix,
    motifEnrichmentDF,
    motifColumn = "TranscriptionFactor",
    ligands = c("TGFB1", "VEGFA"),
    statColumn = "mlog10Padj",
    statThreshold = 2,
    verbose = FALSE,
    numCores = 1
  )

  expect_true(is.data.frame(results))
  expect_equal(sort(results$ligand), c("TGFB1", "VEGFA"))
  expect_true(all(c("adjp_val", "p_val", "PercentSigTF", "PercInNicheNet") %in% colnames(results)))
  expect_true(all(results$p_val >= 0 & results$p_val <= 1, na.rm = TRUE))
})

test_that("MotifSetEnrichmentAnalysis errors for invalid ligands", {
  motifEnrichmentDF <- read.csv(mocha_fixture_path("msea_motif_enrichment.csv"))
  ligandTFMatrix <- matrix(1, nrow = 2, ncol = 1, dimnames = list(c("SPI1", "GATA1"), "TGFB1"))

  expect_error(
    MOCHA::MotifSetEnrichmentAnalysis(
      ligandTFMatrix,
      motifEnrichmentDF,
      motifColumn = "TranscriptionFactor",
      ligands = "NOT_IN_MATRIX",
      statColumn = "mlog10Padj",
      statThreshold = 2
    ),
    "does not appear in NicheNet matrix"
  )
})

test_that("MotifSetEnrichmentAnalysis heavy regression matches fixture expectations", {
  skip_unless_mocha_heavy()
  motifEnrichmentTestDataFP <- Sys.getenv(
    "MOCHA_MSEA_INPUT",
    unset = "/home/jupyter/MOCHA/input_motifenrichment.csv"
  )
  expectedResultsFP <- Sys.getenv(
    "MOCHA_MSEA_EXPECTED",
    unset = "/home/jupyter/MOCHA/results_MSEA.csv"
  )
  skip_if_not(
    file.exists(motifEnrichmentTestDataFP) && file.exists(expectedResultsFP),
    "Heavy MSEA reference files not available"
  )

  ligandTFMatrix <- readRDS(
    "https://zenodo.org/record/3260758/files/ligand_tf_matrix.rds"
  )
  motifEnrichmentDF <- read.csv(motifEnrichmentTestDataFP)
  expectedResults <- read.csv(expectedResultsFP)

  filteredligandTFMatrix <- ligandTFMatrix[
    rownames(ligandTFMatrix) %in% unique(motifEnrichmentDF$TranscriptionFactor),
  ]

  results <- MOCHA::MotifSetEnrichmentAnalysis(
    filteredligandTFMatrix,
    motifEnrichmentDF,
    motifColumn = "TranscriptionFactor",
    ligands = expectedResults$ligand,
    statColumn = "mlog10Padj",
    statThreshold = 2,
    verbose = FALSE
  )

  results <- dplyr::arrange(results, ligand)
  expectedResults <- dplyr::arrange(expectedResults, ligand)
  expect_equal(round(results$adjp_val, 4), round(expectedResults$adjp_val, 4))
})
