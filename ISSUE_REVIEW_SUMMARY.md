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

- **Status:** Addressed (Phase 3)
- **Summary:** `suggestConsensusThreshold()` added in `R/plotConsensus.R`.
  Supports a kneedle elbow method (default) and a second-derivative
  inflection method on the `log10(PeakNumber)` vs reproducibility curve.
  `plotConsensus()` gains a `showSuggested = TRUE` flag that overlays the
  recommendation as a dashed vertical line.
- **Code evidence:** `R/plotConsensus.R` —
  `suggestConsensusThreshold()`, `.suggest_threshold_from_curve()`,
  `.kneedle_index()`, `.second_derivative_index()`.

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

- **Status:** Addressed (Phase 4)
- **Summary:** Two halves addressed.
  1. **Consistency.** Two clear naming drifts unified with backwards-
     compatible deprecation aliases:
     - `getDifferentialAccessibleTiles(cellPopulations =)` now canonical;
       `cellPopulation =` deprecated alias.
     - `getSampleTileMatrix(reproducibilityThreshold =)` now canonical;
       `threshold =` deprecated alias.
     All in-tree callers (tests + vignettes) migrated.
  2. **Workflow guidance.** A new "Tuning parameters" section in the
     `COVID-walkthrough` vignette threads the major decision points
     (`reproducibilityThreshold`, `signalThreshold`, `minZeroDiff`,
     `qValueThreshold`, `method`/`pairColumn`, `dropoutAdjustment`) and
     ends with a cheat-sheet table mapping each knob to the helper that
     inspects it.
- **Code evidence:** `R/getDifferentialAccessibleTiles.R`,
  `R/getSampleTileMatrix.R`, `vignettes/COVID-walkthrough.Rmd`,
  `tests/testthat/test_getDifferentialAccessibleTiles.R` (two new
  alias-coverage tests).

---

## #30 Extend non-parametric functionality

- **Status:** Addressed (Phase 3)
- **Summary:** Two new statistical paths added.
  - `TwoPartPaired()` (in `R/two_part.R`) combines McNemar on the
    zero / non-zero component with a paired Wilcoxon (or paired t-test)
    on the non-zero pairs, combined into a chi-square in the same shape
    as the existing `TwoPart()`.
  - `.twoPart_polr()` (in `R/estimate_differential_accessibility.R`)
    fits a proportional-odds cumulative-logit model on zero / low /
    high-binned values via `MASS::polr`, returning the Wald χ² on the
    group coefficient.
  - Both are wired into `getDifferentialAccessibleTiles()` via new
    `method =` and `pairColumn =` arguments (default behaviour
    unchanged).
- **Code evidence:** `R/two_part.R` (TwoPartPaired),
  `R/estimate_differential_accessibility.R` (method dispatcher +
  `.twoPart_polr`), `R/getDifferentialAccessibleTiles.R` (argument
  plumbing). `MASS` is now in `Suggests`.

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

- **Status:** Addressed (Phase 4)
- **Summary:** New `vignettes/MotifAndAltTSS-walkthrough.Rmd` covering
  the three Figure 4–5 modules:
  1. Motif enrichment (`MotifEnrichment()` +
     `MotifSetEnrichmentAnalysis()`).
  2. Alternative TSS regulation (`getAltTSS()` + `plotRegion()`).
  3. Motif footprinting (`motifFootprint()` /
     `analyzeFootprints()`).
  Chunks use `eval = FALSE` so the vignette builds quickly under
  BiocStyle without requiring a full motif/TF fixture in the package.
- **Code evidence:** `vignettes/MotifAndAltTSS-walkthrough.Rmd`.

