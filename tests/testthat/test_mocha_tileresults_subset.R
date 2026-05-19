test_that("subset() on MochaTileResults matches subsetMOCHAObject", {
  mtr <- MOCHA::asMochaTileResults(MOCHA:::testTileResults)
  legacy <- MOCHA:::testTileResults

  by_fn <- MOCHA::subsetMOCHAObject(
    legacy,
    subsetBy = "celltype",
    groupList = c("C2", "C5")
  )
  by_method <- subset(mtr, cells = c("C2", "C5"))

  expect_equal(names(by_method), names(by_fn))
  expect_s4_class(by_method, "MochaTileResults")
  expect_equal(
    nrow(by_method@metadata$summarizedData),
    nrow(by_fn@metadata$summarizedData)
  )
})
