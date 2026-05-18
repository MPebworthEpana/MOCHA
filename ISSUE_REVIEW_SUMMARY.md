# MOCHA Open Issues: Codebase Review Summary

This document summarizes currently open GitHub issues against the current codebase state.

## Legend

- **Addressed**: Implemented in codebase (issue may still be open on GitHub).
- **Partially Addressed**: Some related implementation exists, but not fully aligned with issue.
- **Open**: Not implemented or not sufficiently addressed in current code.

---

## #166 Common Issues re: plotting, exporting, footprinting, and sharing MOCHA objects

- **Status:** Addressed (kept open intentionally as reference)
- **Summary:** Portability and path-related problems are documented.
- **Code evidence:** `packMOCHA()`, `unpackMOCHA()`, and `updateDirectoryPath()` exist and are exported.
- **Notes:** Issue body itself indicates it is intended as a persistent reference.

---

## #119 Enhancement: Clean up data input

- **Status:** Addressed (backward-compatible)
- **Summary:** `callOpenTiles()` now accepts sample-level `GRangesList` input (one element per sample) and derives `CellPopulation#Sample` fragments from `cellColData`.
- **Code evidence:** `.detectATACFragmentsInputMode()`, `.normalizeSampleLevelFragments()` in `R/callOpenTiles.R`; tests in `tests/testthat/test_callOpenTiles.R`.
- **Notes:** Legacy `CellPopulation#Sample` naming remains supported.

---

## #109 MOCHA needs getters and subsetters

- **Status:** Partially Addressed
- **Summary:** Subsetting and multiple getters are implemented.
- **Code evidence:** `subsetMOCHAObject()`, `getCellTypes()`, `getCellTypeTiles()`, `getSampleCellTypeMetadata()`.
- **Gap:** The specific getter named in issue (`getOpenTiles`) is not present.

---

## #85 adding automated thresholding for consensus tiles

- **Status:** Open
- **Summary:** Consensus plotting exists, but no automated derivative/changepoint threshold recommendation found.
- **Code evidence:** `plotConsensus()` exists.
- **Gap:** No automated threshold-selection algorithm matching issue request.

---

## #84 addCellColData for MOCHA TileResults

- **Status:** Open
- **Summary:** No `addCellColData` helper found in exported API.
- **Gap:** Functionality requested by issue appears absent.

---

## #69 Changing dependency from AnnotationDbi/RMariaDb to biomaRt

- **Status:** Open
- **Summary:** Current package still uses `AnnotationDbi`; `RMariaDB` still appears in `Suggests`.
- **Code evidence:** `DESCRIPTION` includes `AnnotationDbi` import and `RMariaDB` suggestion.
- **Gap:** Requested dependency transition not completed.

---

## #61 Handling optional parameters for downstream steps

- **Status:** Partially Addressed
- **Summary:** Some optional tuning parameters and plotting utilities exist.
- **Code evidence:** `minZeroDiff`, `signalThreshold`, `plotConsensus()`, `plotIntensityDistribution()`.
- **Gap:** Broader consistency and workflow guidance requested in issue appears only partially covered.

---

## #36 Plot consensus tiles at different thresholds

- **Status:** Addressed
- **Summary:** Functionality for reproducibility-threshold exploration is implemented.
- **Code evidence:** `plotConsensus()`.

---

## #31 Assessing Dropout: Technical vs. Biological 0s

- **Status:** Open
- **Summary:** Basic dropout filtering exists, but no broad expanded module specifically separating technical vs biological zero inflation was identified.
- **Code evidence:** `minZeroDiff` in differential testing.
- **Gap:** Requested expanded dropout modeling functionality appears unresolved.

---

## #30 Extend non-parametric functionality

- **Status:** Open
- **Summary:** Requested paired 2-part tests / proportional-odds extensions were not found.
- **Gap:** No direct implementation of the issue's specific proposed methods.

---

## #29 Longitudinal, LMM, and Variance decomposition integration

- **Status:** Partially Addressed
- **Summary:** Modeling functions exist, including LMEM/ZIGLMM/variance decomposition.
- **Code evidence:** `runLMEM()`, `runZIGLMM()`, `varZIGLMM()`.
- **Gap / Risk:** These functions are marked deprecated in `NEWS.md`; `varZIGLMM()` currently contains a `browser()` debug breakpoint.

---

## #27 Draft Vignettes for Figure 4-5 Analyses

- **Status:** Partially Addressed
- **Summary:** Related analysis functions exist (motif enrichment, alt TSS, etc.), but dedicated figure-4/5-style vignette coverage appears incomplete.
- **Code evidence:** Functions like `getAltTSS()`, `MotifEnrichment()`, `MotifSetEnrichmentAnalysis()`.
- **Gap:** Vignette scope does not clearly map to issue's requested downstream analysis drafts.

---

## Notes

- Open PRs can appear in issue listings via GitHub API and should be excluded from issue triage summaries.
- This summary reflects codebase state, not project-management intent (some issues may intentionally remain open).
