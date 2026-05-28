# Assess technical vs biological zeros on a sample-tile matrix

`assessDropout` fits a dropout model and stores predicted `P(zero)` per
tile and sample in `DropoutProb_<CellPop>` assays. \*\*Lower\*\* values
indicate a zero is more likely technical (unexpected, given tile signal
and depth). Model objects are saved in `metadata(TSAM)$dropoutModels`.

## Usage

``` r
assessDropout(
  TSAM_Object,
  model = NULL,
  cellPopulation = "all",
  techThreshold = 0.2,
  bioThreshold = 0.8,
  verbose = FALSE
)
```

## Arguments

- TSAM_Object:

  A SummarizedExperiment from
  [`getSampleTileMatrix`](https://aifimmunology.github.io/MOCHA/reference/getSampleTileMatrix.md).

- model:

  Optional `DropoutModel` from `estimateDropoutModel`. If `NULL`, a
  model is estimated with default settings.

- cellPopulation:

  Cell population name(s) to assess, or `"all"` for every population in
  `model`. When estimating a new model, only these populations are
  fitted.

- techThreshold:

  Maximum `P(zero)` for classifying an observed zero as technical
  (unexpected). Default is 0.2.

- bioThreshold:

  Minimum `P(zero)` for classifying an observed zero as biological
  (expected closure). Default is 0.8.

- verbose:

  Logical. Default is `FALSE`.

## Value

The input `TSAM_Object` with `DropoutProb_*` assays added.

## Examples

``` r
if (FALSE) { # \dontrun{
TSAM_with_dropout <- MOCHA::assessDropout(SampleTileMatrices)
probs <- MOCHA::getDropoutProb(TSAM_with_dropout, "CD16 Mono")
} # }
```
