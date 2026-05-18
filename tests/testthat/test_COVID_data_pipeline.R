skip_on_cran()
skip_unless_mocha_heavy()

if (
  requireNamespace("TxDb.Hsapiens.UCSC.hg38.refGene", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
    requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE) &&
    requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE) &&
    requireNamespace("ArchR", quietly = TRUE)
) {
  ArchRProjDir <- mocha_archr_project_dir("FullCovid")
  skip_if(is.na(ArchRProjDir), "FullCovid ArchR project not available")

  oldw <- getOption("warn")
  options(warn = -1)
  test_that("We reproduce the COVID CD16 Monocyte analysis", {
    ArchRProj <- ArchR::loadArchRProject(ArchRProjDir)
    metadata <- data.table::as.data.table(ArchR::getCellColData(ArchRProj))
    studySignal <- median(metadata$nFrags)

    expect_equal(studySignal, 3628)

    lookup_table <- unique(
      metadata[, c(
        "Sample",
        "COVID_status",
        "Visit",
        "days_since_symptoms"
      ),
      with = FALSE
      ]
    )

    samplesToKeep <- lookup_table$Sample[
      lookup_table$Visit == "FH3 COVID-19 Visit 1" &
        lookup_table$days_since_symptoms <= 15 |
        is.na(lookup_table$days_since_symptoms)
    ]

    idxSample <- BiocGenerics::which(ArchRProj$Sample %in% samplesToKeep)
    cellsSample <- ArchRProj$cellNames[idxSample]
    ArchRProj <- ArchRProj[cellsSample, ]

    capture.output(tileResults <- MOCHA::callOpenTiles(
      ArchRProj,
      cellPopLabel = "CellSubsets",
      cellPopulations = "CD16 Mono",
      TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene",
      Org = "org.Hs.eg.db",
      numCores = 2,
      studySignal = studySignal,
      outDir = tempdir()
    ))

    tileResultAssay <- RaggedExperiment::compactAssay(
      MultiAssayExperiment::experiments(tileResults)[[1]],
      i = "TotalIntensity"
    )
    expect_gt(nrow(tileResultAssay), 0L)
    expect_gt(ncol(tileResultAssay), 0L)

    capture.output(SampleTileMatrices <- MOCHA::getSampleTileMatrix(
      tileResults,
      cellPopulations = "CD16 Mono",
      groupColumn = "COVID_status",
      threshold = 0.2,
      verbose = FALSE
    ))

    SampleTileMatrixAssay <- SummarizedExperiment::assays(SampleTileMatrices)[[1]]
    expect_gt(nrow(SampleTileMatrixAssay), 0L)

    capture.output(cd16matrix <- MOCHA::getCellPopMatrix(
      SampleTileObj = SampleTileMatrices,
      cellPopulation = "CD16 Mono",
      dropSamples = TRUE,
      NAtoZero = TRUE
    ))

    expect_gt(nrow(cd16matrix), 0L)

    capture.output(
      differentials <- MOCHA::getDifferentialAccessibleTiles(
        SampleTileObj = SampleTileMatrices,
        cellPopulation = "CD16 Mono",
        groupColumn = "COVID_status",
        foreground = "Positive",
        background = "Negative",
        outputGRanges = TRUE,
        numCores = 2
      )
    )
    cd16_differentials <- plyranges::filter(differentials, FDR <= 0.2)

    expect_s4_class(cd16_differentials, "GRanges")
    expect_equal(length(cd16_differentials), 6211)

    FDRSummary <- summary(cd16_differentials$FDR)
    expect_equal(round(FDRSummary[["Min."]], 11), 0.01887700776)
    expect_equal(round(FDRSummary[["Max."]], 7), 0.1983282)
  })
  options(warn = oldw)
}
