# MOCHA exported API test coverage matrix

Generated for test-suite cleanup. **Reference** means the symbol appears in `tests/testthat/test_*.R`.
**Dedicated** means a test file primarily exercises that API.

| Export | Risk | Reference | Dedicated / notes |
|--------|------|-----------|-------------------|
| `callOpenTiles` | High | Yes | `test_callOpenTiles.R` |
| `getSampleTileMatrix` | High | Yes | `test_getSampleTileMatrix.R` |
| `getDifferentialAccessibleTiles` | High | Partial | `test_getDifferentialAccessibleTiles.R` (restored) |
| `getCoAccessibleLinks` | High | Yes | `test_getCoAccessibleLinks.R` |
| `testCoAccessibility` | High | Yes | `test_testCoAccessibility.R` |
| `runLMEM` / `pilotLMEM` | High | Yes | `test_runLMEM.R` (errors only) |
| `runZIGLMM` / `pilotZIGLMM` / `varZIGLMM` | High | Yes | `test_modeling_ziglmm.R` |
| `linearModeling` | High | Yes | `test_modeling_ziglmm.R` |
| `getModelValues` | Medium | Yes | `test_modeling_ziglmm.R` |
| `combineSampleTileMatrix` | Medium | Yes | `test_combineSampleTileMatrix.R` |
| `subsetMOCHAObject` | Medium | Yes | `test_subsetMOCHAObject.R` |
| `annotateTiles` | Medium | Yes | `test_annotateTiles.R` |
| `MotifEnrichment` | Medium | Yes | `test_MotifEnrichment.R` |
| `MotifSetEnrichmentAnalysis` | Medium | Yes | `test_MotifSetEnrichmentAnalysis.R` (fixtures) |
| `filterCoAccessibleLinks` | Medium | Yes | `test_filterCoAccessibleLinks.R` |
| `plotConsensus` | Medium | Yes | `test_plotConsensus.R` |
| `getSampleCellTypeMetadata` | Medium | Yes | `test_utils_metadata.R` |
| `StringsToGRanges` / `GRangesToString` | Low | Yes | `test_utils_strings.R` |
| `isMOCHAObject` | Low | Yes | `test_utils_strings.R` |
| `differentialsToGRanges` | Low | Yes | `test_MotifEnrichment.R` |
| `packMOCHA` / `unpackMOCHA` | Medium | Yes | `test_packMOCHA.R` |
| `exportCoverage` / `exportOpenTiles` / `exportMotifs` | Medium | Yes | export tests |
| `getCoverage` | Medium | Yes | `test_getCoverage.R` |
| `extractRegion` | Medium | Partial | `test_extractRegion.R` (heavy data optional) |
| `getPopFrags` | Medium | Partial | `test_getPopFrags.R` (ArchR optional) |
| `import_snap_atac` | Medium | Yes | `test_import_snap_atac.R` (reticulate/anndata optional) |
| `bulkDimReduction` | Medium | Yes | `test_dimensionalityReduction.R` |
| `bulkUMAP` | Low | — | allowlisted (optional uwot) |
| `plotRegion` | Medium | Heavy | `test_plotRegion.R` (`MOCHA_HEAVY_TESTS`) |
| `addMotifSet` | Medium | Yes | `test_addMotifSet.R` |
| `getCellPopMatrix` | Medium | Yes | `test_getCellPopMatrix.R` |
| `testCoAccessibilityChromVar` | Medium | — | allowlisted (chromVAR optional) |
| `testCoAccessibilityRandom` | Medium | — | allowlisted (optional) |
| `addAccessibilityShift` / `addInsertionBias` | Low | — | allowlisted |
| `correctGenome` | Low | — | allowlisted |
| `exportDifferentials` / `exportLocalFootprints` | Low | — | allowlisted |
| `getAltTSS` / `getPromoterGenes` | Low | — | allowlisted |
| `getCellTypeMotifs` / `getCellTypeTiles` / `getCellTypes` | Low | — | allowlisted |
| `mergeTileResults` | Low | — | allowlisted |
| `motifFootprint` | Medium | — | allowlisted (heavy) |
| `plotMotifs` / `plotIntensityDistribution` | Low | — | allowlisted |
| `renameCellTypes` | Low | — | allowlisted |
| `updateDirectoryPath` | Low | Yes | `test_updateDirectoryPath.R` |

## Test profiles

- **Default (CI):** `R CMD check` / `testthat::test_local()` — uses bundled `testTileResults*` and fixtures; skips `MOCHA_HEAVY_TESTS`.
- **Heavy:** set `MOCHA_HEAVY_TESTS=true` and provide external ArchR projects (`PBMCSmall`, `FullCovid`, `HemeTutorial/MOCHA`) as documented in `tests/testthat/helper-mocha-test-profiles.R`.
