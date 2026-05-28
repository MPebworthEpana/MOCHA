# Suggest a reproducibility threshold for consensus tile selection

`suggestConsensusThreshold` examines the reproducibility-vs-peak-number
curve produced internally by
[`plotConsensus()`](https://aifimmunology.github.io/MOCHA/reference/plotConsensus.md)
and returns a recommended threshold per cell population (and optionally
per group within population). Two automated methods are available:
`"kneedle"` finds the point farthest from the secant joining the
endpoints (the classic knee detection algorithm); `"second_derivative"`
returns the reproducibility value at which the numerical second
difference of `log10(PeakNumber)` is maximised in magnitude (the
inflection in steepness).

## Usage

``` r
suggestConsensusThreshold(
  tileObject,
  cellPopulations = "all",
  groupColumn = NULL,
  method = c("kneedle", "second_derivative"),
  numCores = 1
)
```

## Arguments

- tileObject:

  A MultiAssayExperiment object from
  [`callOpenTiles()`](https://aifimmunology.github.io/MOCHA/reference/callOpenTiles-methods.md).

- cellPopulations:

  Cell populations to evaluate (default `"all"`).

- groupColumn:

  Optional grouping column from `colData(tileObject)`. If supplied, one
  threshold is returned per (population, group).

- method:

  One of `"kneedle"` (default) or `"second_derivative"`.

- numCores:

  Cores passed through to
  [`plotConsensus()`](https://aifimmunology.github.io/MOCHA/reference/plotConsensus.md).

## Value

A `data.frame` with columns `CellPopulation`, `Reproducibility`
(suggested threshold), `PeakNumber` (peaks retained at that threshold),
`Method`, and `GroupName` when a grouping is requested.
