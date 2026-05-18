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
  test_that("plotRegion runs end-to-end on COVID subset", {
    ArchRProj <- ArchR::loadArchRProject(ArchRProjDir)
    metadata <- data.table::as.data.table(ArchR::getCellColData(ArchRProj))
    studySignal <- median(metadata$nFrags)
    expect_equal(studySignal, 3628)

    lookup_table <- unique(
      metadata[, c("Sample", "COVID_status", "Visit", "days_since_symptoms"), with = FALSE]
    )
    samplesToKeep <- lookup_table$Sample[
      lookup_table$Visit == "FH3 COVID-19 Visit 1" &
        lookup_table$days_since_symptoms <= 15 |
        is.na(lookup_table$days_since_symptoms)
    ]
    idxSample <- BiocGenerics::which(ArchRProj$Sample %in% samplesToKeep)
    ArchRProj <- ArchRProj[ArchRProj$cellNames[idxSample], ]

    mytempdir <- tempdir()
    capture.output(tileResults <- MOCHA::callOpenTiles(
      ArchRProj,
      cellPopLabel = "CellSubsets",
      cellPopulations = "CD16 Mono",
      TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene",
      Org = "org.Hs.eg.db",
      numCores = 2,
      studySignal = studySignal,
      outDir = mytempdir,
      verbose = FALSE
    ))

    capture.output(SampleTileMatrices <- MOCHA::getSampleTileMatrix(
      tileResults,
      cellPopulations = "CD16 Mono",
      groupColumn = "COVID_status",
      threshold = 0.2,
      verbose = FALSE
    ))

    region <- "chr4:122610000-122625000"
    capture.output(
      countSE <- MOCHA::extractRegion(
        SampleTileObj = SampleTileMatrices,
        region = region,
        cellPopulations = "CD16 Mono",
        groupColumn = "COVID_status",
        sampleSpecific = FALSE,
        binSize = 250,
        numCores = 2,
        verbose = FALSE
      )
    )

    pdf(file.path(mytempdir, "testplotregion.pdf"))
    MOCHA::plotRegion(countSE = countSE)
    dev.off()
    expect_true(file.exists(file.path(mytempdir, "testplotregion.pdf")))
  })
  options(warn = oldw)
}
