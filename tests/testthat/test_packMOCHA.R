skip_on_cran()
if (
  require("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE) &&
  require("org.Hs.eg.db", quietly = TRUE) &&
  require("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE) &&
  require("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
) {
  if (mocha_heavy_tests_enabled() && requireNamespace("ArchR", quietly = TRUE)) {
    ArchRProjDir <- mocha_archr_project_dir("PBMCSmall")
    if (!is.na(ArchRProjDir)) {
    test_that("We can pack/unpack a MOCHA object", {
      capture.output(
        testProj <- ArchR::loadArchRProject(ArchRProjDir),
        type = "message"
      )

      TxDb <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
      OrgDb <- "org.Hs.eg.db"
      
      mytempdir <- tempdir()
      
      capture.output(
        tiles <- MOCHA::callOpenTiles(
          ATACFragments = testProj,
          TxDb = TxDb,
          OrgDb = OrgDb,
          cellPopLabel = "Clusters",
          cellPopulations = c("C2", "C5"),
          numCores = 1,
          outDir = mytempdir
        ),
        type = "message"
      )
      
      zipPath <- MOCHA::packMOCHA(tiles, zipfile = file.path(mytempdir, "testzip.zip"))
      expect_true(file.exists(zipPath))
      expect_error(
        MOCHA::packMOCHA(tiles, zipfile = file.path(mytempdir, "testzip.zap"))
      )
      
      if (requireNamespace("waldo")) {
        unpackedmochaObj <- MOCHA::unpackMOCHA(zipPath, mytempdir)
        diff <- waldo::compare(unpackedmochaObj, tiles)
        expect_length(diff, 1)
        expect_true(grepl('metadata$Directory', diff[1], fixed = TRUE))
      }
    })
    }
  }
}
