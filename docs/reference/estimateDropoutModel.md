# Estimate a logistic dropout model for a cell population

`estimateDropoutModel` fits a per-cell-population logistic model
relating the probability of an observed zero to tile signal and sample
depth, to distinguish technical dropout from biological closure.

## Usage

``` r
estimateDropoutModel(
  TSAM_Object,
  cellPopulation = "all",
  additionalCovariates = NULL,
  depthCols = c("FragNumber", "CellCounts"),
  minNonzero = 5,
  verbose = FALSE
)
```

## Arguments

- TSAM_Object:

  A SummarizedExperiment from
  [`getSampleTileMatrix`](https://aifimmunology.github.io/MOCHA/reference/getSampleTileMatrix.md).

- cellPopulation:

  Character vector of assay (cell population) names, or `"all"` for
  every intensity assay (excluding existing `DropoutProb_*` assays).

- additionalCovariates:

  Optional character vector of `colData` column names to include in the
  model.

- depthCols:

  Character vector of depth covariate names. Values are read from
  `colData` when present; otherwise from `metadata(TSAM)$summarizedData`
  via fragment and cell count tables. `FragmentCounts` is aliased to
  `FragNumber`.

- minNonzero:

  Minimum number of non-zero observations required to include a tile in
  model fitting. Default is 5.

- verbose:

  Logical. Default is `FALSE`.

## Value

An S3 object of class `DropoutModel` containing fitted models, formulas,
depth column names, and per-cell-population diagnostics (including
leave-one-sample-out AUROC with train-only tile means per fold).

## Examples

``` r
if (FALSE) { # \dontrun{
model <- MOCHA::estimateDropoutModel(SampleTileMatrices, cellPopulation = "CD16 Mono")
} # }
```
