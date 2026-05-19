# MOCHA Open Issues: Remaining Work (Codebase Review)

This document lists open GitHub issues that are **not yet fully addressed** in the current codebase. Issues implemented in code (even if still open on GitHub) are omitted.

## Legend

- **Partially Addressed**: Some related implementation exists, but not fully aligned with the issue.
- **Open**: Not implemented or not sufficiently addressed in current code.

---

## #109 MOCHA needs getters and subsetters

- **Status:** Addressed (Phase 1)
- **Summary:** `getOpenTiles()` added in `R/utils.R`, alongside the existing
  `subsetMOCHAObject()`, `getCellTypes()`, `getCellTypeTiles()`, and
  `getSampleCellTypeMetadata()` helpers. Returns a `GRangesList` (default) or
  flat `data.frame` of called peaks per cell population.

---

## #85 adding automated thresholding for consensus tiles

- **Status:** Open
- **Summary:** Consensus plotting exists, but no automated derivative/changepoint threshold recommendation found.
- **Code evidence:** `plotConsensus()` exists.
- **Gap:** No automated threshold-selection algorithm matching issue request.

---

## #84 addCellColData for MOCHA TileResults

- **Status:** Addressed (Phase 1)
- **Summary:** `addCellColData()` added in `R/utils.R`. Works on both MOCHA
  tileResults (`MultiAssayExperiment`) and SampleTileMatrix
  (`SummarizedExperiment`) objects. Supports positional value alignment or a
  `samples =` mapping, with `force = TRUE` to overwrite.

---

## #69 Changing dependency from AnnotationDbi/RMariaDb to biomaRt

- **Status:** Addressed (Phase 2)
- **Summary:** `AnnotationDbi` moved from `Imports` to `Suggests`; `biomaRt`
  added to `Suggests`; `RMariaDB` removed entirely. A new internal helper
  `.map_gene_ids()` in `R/utils.R` routes through `biomaRt::getBM()` when a
  `Mart` object is supplied and falls back to `AnnotationDbi::mapIds()`
  otherwise, preserving offline workflows for users with an `OrgDb`
  installed.
- **Code evidence:** `R/utils.R:6` (`.map_gene_ids`), call sites at
  `R/annotateTiles.R:58`, `R/getAltTSS.R:79`,
  `R/plottingUtils.R:299,309`.

---

## #61 Handling optional parameters for downstream steps

- **Status:** Partially Addressed
- **Summary:** Some optional tuning parameters and plotting utilities exist.
- **Code evidence:** `minZeroDiff`, `signalThreshold`, `plotConsensus()`, `plotIntensityDistribution()`.
- **Gap:** Broader consistency and workflow guidance requested in issue appears only partially covered.

---

## #30 Extend non-parametric functionality

- **Status:** Open
- **Summary:** Requested paired 2-part tests / proportional-odds extensions were not found.
- **Gap:** No direct implementation of the issue's specific proposed methods.

---

## #29 Longitudinal, LMM, and Variance decomposition integration

- **Status:** Addressed (Phase 1)
- **Summary:** Modeling functions remain exported but now emit
  `lifecycle::deprecate_warn()` / `deprecate_soft()` at call time (previously
  commented out). The `browser()` breakpoint in `varZIGLMM()` is no longer
  present. Users are pointed to the `ChAI` package for improved modeling.
- **Code evidence:** `runLMEM()`, `runZIGLMM()`, `varZIGLMM()`,
  `pilotLMEM()`, `pilotZIGLMM()`, `individualZIGLMM()`, `getModelValues()`,
  `processModelOutputs()`.

---

## #27 Draft Vignettes for Figure 4-5 Analyses

- **Status:** Partially Addressed
- **Summary:** Related analysis functions exist (motif enrichment, alt TSS, etc.), but dedicated figure-4/5-style vignette coverage appears incomplete.
- **Code evidence:** Functions like `getAltTSS()`, `MotifEnrichment()`, `MotifSetEnrichmentAnalysis()`.
- **Gap:** Vignette scope does not clearly map to issue's requested downstream analysis drafts.

