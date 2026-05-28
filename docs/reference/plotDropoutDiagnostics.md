# Plot dropout model diagnostics

Plot dropout model diagnostics

## Usage

``` r
plotDropoutDiagnostics(
  TSAM_Object,
  cellPopulation,
  type = c("calibration", "probHist", "perSample")
)
```

## Arguments

- TSAM_Object:

  A SummarizedExperiment with dropout results from
  [`assessDropout`](https://aifimmunology.github.io/MOCHA/reference/assessDropout.md).

- cellPopulation:

  Cell population to plot.

- type:

  One of `"calibration"`, `"probHist"`, or `"perSample"`.

## Value

A `ggplot2` object.

## Examples

``` r
if (FALSE) { # \dontrun{
MOCHA::plotDropoutDiagnostics(TSAM_with_dropout, "CD16 Mono", type = "calibration")
} # }
```
