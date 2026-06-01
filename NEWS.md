# MOCHA 2.0.0

* Prepare package for initial Bioconductor submission on branch `bioc/submission`.
* Bump package version to `2.0.0`; populate `biocViews`, `BiocType`, and remove CRAN-only `Additional_repositories`.
* Replace `TxDb.Hsapiens.UCSC.hg38.refGene` with Bioconductor `TxDb.Hsapiens.UCSC.hg38.knownGene`.
* Add runnable `COVID-walkthrough` vignette (BiocStyle) using bundled example data.
* Add `inst/CITATION`, `Bioc-check` CI workflow, and `requireNamespace()` guards for optional ArchR usage.
* Retire CRAN submission artifacts (`cran-comments.md`, `CRAN-SUBMISSION`).

# MOCHA (development version)

## Documentation

* Renamed vignettes: `MOCHA-workflow-tutorial`, `Data-Import-Tutorial`, `Alternative-TSS-TF-regulation`.
* Merged `ImportingFromOtherSources` into `Data-Import-Tutorial` (ArchR, Signac, SnapATAC2, legacy SnapATAC v1, and manual fragment workflows).
* pkgdown redirects from former article URLs to the new vignette pages.
* Expanded tutorial library: tuning parameters and differential options in the workflow vignette;
  new articles `MOCHA-downstream-workflows`, `MOCHA-export-and-sharing`, and
  `MOCHA-advanced-modeling`; expanded alt-TSS/motif and import-utility sections.
* Bioconductor build policy: only bundled-data chunks execute on builders; other vignette
  chunks are reference-only (`eval = FALSE`).
* Tutorial code lives in `inst/tutorials/*.R` (single source of truth); vignettes include
  scripts via `knitr::read_chunk()`. Verified by `tests/scripts/run_tutorials.R` and
  `.github/workflows/tutorial-check.yml`.

## Coverage storage

* Coverage tracks from `callOpenTiles()` are now saved as per-sample bigWig
  files under `tracks/{cellPop}/{Accessibility|Insertions}/`, with legacy
  `{cellPop}_CoverageFiles.RDS` bundles still supported for reading.
* Added `migrateCoverageToTracks()` to convert legacy RDS bundles to the
  structured layout.
* Track filenames use reversible percent-encoding to avoid sample-name
  collisions; `exportCoverage()` copies existing structured tracks when
  exporting sample-specific bigWigs without regrouping.
* Coverage I/O validates unsafe path characters, prunes orphaned bigWigs
  when sample sets change, and repairs incomplete migrations from legacy RDS
  bundles when available.

## Bug fixes

* Fixed `getSampleCellTypeMetadata` export (`@export` was incorrectly on internal `.get_sample_celltype_count_tables()`).
* Fixed `runZIGLMM()` / `pilotZIGLMM()` when `cellPopulation = "counts"` on combined sample-tile matrices (aligned with `varZIGLMM()`).
* Fixed `varZIGLMM()` undefined `cl` when `numCores == 1`.
* Fixed matrix subsetting (`drop = FALSE`) in `linearModeling()`, `varZIGLMM()`, and `runZIGLMM()` preventing empty modeling inputs.
* Fixed `plotConsensus()` / `cellTypeDF()` when a single sample remains after filtering (`drop = FALSE`, empty-column guard).
* Fixed `filterCoAccessibleLinks()` when all correlations are `NA` (`na.rm = TRUE` on threshold check).
* Fixed `StringsToGRanges()` to error on unparseable region strings instead of returning invalid ranges.

* MOCHA S4 subclasses (opt-in, backward compatible):
  - `MochaTileResults` and `MochaSampleTileMatrix` extend the existing
    Bioconductor output types with validity checks and optional `returnClass =
    "mocha"` on constructors.
  - `asMochaTileResults()`, `asMochaSTM()`, and subclass methods `cellTypes()`
    and `openTiles()` (tileResults only).
  - `subset()` on subclass objects mirrors `subsetMOCHAObject()`; legacy defaults
    and `subsetMOCHAObject()` are unchanged.
* New Functions:
  - `trainPeakModel`: retrain the MOCHA peak-calling LRM and Youden threshold
    at a user-chosen tile size, following the paper's Methods. Optional
    dependencies: `cutpointr` for threshold selection and a MACS2 binary (or
    user-supplied `groundTruthPeaks`) for training labels.
  - `callOpenTiles` gains a `peakModel` argument; default behaviour is unchanged.
* New dropout assessment module for distinguishing technical vs biological zeros
  in sample-tile matrices:
  - `estimateDropoutModel()`
  - `assessDropout()`
  - `getDropoutProb()`
  - `classifyZeros()`
  - `plotDropoutDiagnostics()`
* New getters / setters:
  - `addCellColData()` — append a sample-level column to a MOCHA tileResults
    or SampleTileMatrix object's colData (issue #84).
  - `getOpenTiles()` — extract per-cell-population called peaks from a
    tileResults `MultiAssayExperiment` as a `GRangesList` or flat
    `data.frame` (issue #109).
* `getDifferentialAccessibleTiles()` gains optional `dropoutAdjustment` and
  `bioThreshold` arguments (appended after `verbose` to preserve positional
  compatibility). Deprecated aliases: `techThreshold` → `bioThreshold`,
  `fdrToDisplay` → `qValueThreshold`. Dropout-adjusted modes fail fast with a
  clear error when `DropoutProb_*` assays cannot be fitted.
* `subsetMOCHAObject()` retains matching `DropoutProb_*` assays when subsetting
  by cell type.
* **Breaking change:** `DropoutProb_*` assays now store fitted `P(zero)` (not
  `1 - P(zero)`). Lower values indicate technical (unexpected) zeros; higher
  values indicate biological closure. `classifyZeros()` defaults are now
  `techThreshold = 0.2` and `bioThreshold = 0.8`. Prior releases stored a
  technical score where higher meant more technical.
* Dropout diagnostics report leave-one-sample-out AUROC with tile means
  (`log_mu`) recomputed from training samples only within each fold.
* `plotDropoutDiagnostics(type = "calibration")` requires
  `metadata(TSAM)$dropoutModels` (no zero-only fallback plot).
* `runLMEM()`, `runZIGLMM()`, `varZIGLMM()`, and their companion functions
  now emit `lifecycle::deprecate_warn()` at call time, matching their
  documented deprecated status (issue #29).
* Dependency migration (issue #69):
  - `AnnotationDbi` moved from `Imports` to `Suggests`; `biomaRt` added to
    `Suggests`; `RMariaDB` removed.
  - Internal helper `.map_gene_ids()` routes gene-id mapping through either
    `AnnotationDbi::mapIds()` (offline OrgDb) or `biomaRt::getBM()` when a
    `Mart` connection is passed. Existing `annotateTiles()`, `getAltTSS()`,
    and `plotRegion()` calls work unchanged for users with an OrgDb
    installed.
* Consensus thresholding (issue #85):
  - `suggestConsensusThreshold()` returns an automated reproducibility
    threshold per cell population using either a kneedle elbow or a
    second-derivative inflection on the `log10(PeakNumber)` curve.
  - `plotConsensus()` gains `showSuggested = TRUE` to overlay the
    recommendation as a dashed vertical line.
* Non-parametric extensions (issue #30):
  - New paired two-part test (`TwoPartPaired()`) combining McNemar on the
    binary component with a paired Wilcoxon (or paired t) on the
    non-zero pairs.
  - `getDifferentialAccessibleTiles()` gains `method =` and `pairColumn =`
    arguments. Choose `"wilcoxon"` (default, unpaired), `"paired_wilcoxon"`
    (paired via `pairColumn`), or `"polr"` (proportional-odds cumulative
    logit; requires the `MASS` package, added to `Suggests`).
* Parameter-naming consistency sweep (issue #61):
  - `getDifferentialAccessibleTiles(cellPopulations = )` is now the
    canonical kwarg; `cellPopulation =` is a deprecated alias that emits
    `lifecycle::deprecate_warn("2.0.0", ...)`.
  - `getSampleTileMatrix(reproducibilityThreshold = )` is now canonical;
    `threshold =` is a deprecated alias.
  - All in-tree tests and the `COVID-walkthrough` vignette were migrated
    to the new canonical names. Existing user code continues to work
    until the next major version.
* New `MotifAndAltTSS-walkthrough` vignette (issue #27) — motif
  enrichment, alternative TSS regulation, and motif footprinting
  reference using bundled example data shape.
* `COVID-walkthrough` vignette gains a "Tuning parameters" section
  threading `plotConsensus()` / `suggestConsensusThreshold()` →
  `plotIntensityDistribution()` → `signalThreshold` → method /
  `pairColumn` → dropout adjustment, with a cheat-sheet table mapping
  each knob to the helper that inspects it (issue #61).
* Internal refactor to reduce duplicated helper logic across coverage extraction,
  export, motif footprinting, and co-accessibility workflows. No intended
  user-facing behavior changes.

# MOCHA 1.1.0
* New Functions:
  - Sharing MOCHA objects between file systems
    - packMOCHA
    - unpackMOCHA
  - Exporting for genome browsers (bigwig, bigbed)
    - exportCoverage
    - exportDifferentials
    - exportMotifs
    - exportOpenTiles
    - exportSmoothedInsertions
  - Getters 
    - getCellTypeTiles
    - getCellTypes
    - getPromoterGenes
    - getSampleCellTypeMetadata
  - Other
    - mergeTileResults
    - plotIntensityDistribution
    - renameCellTypes

* Deprecated Functions:
  - runLMEM, pilotLMEM, runZIGLMM, pilotZIGLMM, IndividualZIGLMM, getModelValues, varZIGLMM, processModelOutputs
  
* Updates test data to consistently use TxDb hg19 references.

* Various minor bug and documentation fixes.

# MOCHA 1.0.2
* Addressed check errors in "donttest" examples.

# MOCHA 1.0.1

* Deprecates `testCoAccessibilityChromVAR()` and `testCoAccessibilityRandom()` in favor of `testCoAccessibility()`
* Updates maintainer email

# MOCHA 1.0.0

* Adopting semantic versioning starting with this version, versioning reflects breaking changes compared to previous CRAN release.
* Includes test improvements
* New functions bulkDimReduction, bulkUMAP, MotifEnrichment, MotifSetEnrichmentAnalysis, pilotLMEM, runLMEM, pilotZIGLMM, runZIGLMM, combineSampleTileMatrix, getCoverage
* Improvements to metadata carried by output objects

# MOCHA 0.2.5

* patches a bug (rounding error) found in getDifferentialAccessibleTiles (#125), and reverts to using mclapply parallelization forgetDifferentialAccessibleTiles.
* adds conditional tests (and snapshot tests) on the COVID dataset to ensure reproducibility with results in the MOCHA manuscript
* updates the COVID vignette through differentials to reflect the latest usage.

# MOCHA 0.2.4

* Fixes bug in callOpenTiles where "Clusters" was hardcoded in the step computing fragment counts table.
* Parallelization overhaul to address memory leaks when using parLapply. ParLapply is now passed a helper function directly and a single object input to that function, where the input object is a list containing all variables needed in the helper function.

# MOCHA 0.2.3

* MOCHA MultiAssayExperiment and MOCHA SummarizedExperiment objects now contain new metadata 
* CallOpenTiles now only only accepts the database package names (strings) for Genome, OrgDb, and TxDb, and not the in-memory R objects for those database packages
* Downstream functions where Genome, OrgDb, and TxDb were previously inputs now check the input MOCHA object's metadata to load the relevant databases.

# MOCHA 0.2.2

* getCoAccessibleLinks
* testCoAccessibilityChromVar and testCoAccessibilityRandom
* combineSampleTileMatrix

# MOCHA 0.2.1

* This release adds additional test coverage with new test data, covering edge cases in callOpenTiles.

* Removes option log2Intensity from getSampleTileMatrix (done by default in getDifferentialAccessibleTiles)

# MOCHA 0.2.0

* This includes the initial release of package on CRAN, adds updated requirements for R >= 4.1.0 and plyranges >1.14.0

# MOCHA 0.1.0

* Added a `NEWS.md` file to track changes to the package.
* MOCHA is submitted to CRAN as an initial release.
