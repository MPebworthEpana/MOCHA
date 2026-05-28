# Get per-cell-population open tiles from a MOCHA tileResults object

`getOpenTiles` extracts the called open tiles (peaks) for one or more
cell populations from a `MultiAssayExperiment` returned by
`callOpenTiles`. By default tiles are returned as a `GRangesList` keyed
by cell population; with `returnType = "data.frame"` they are flattened
into a single data frame with a `CellPopulation` column.

## Usage

``` r
getOpenTiles(
  tileResults,
  cellPopulations = "all",
  returnType = c("GRangesList", "data.frame")
)
```

## Arguments

- tileResults:

  A `MultiAssayExperiment` from `callOpenTiles`.

- cellPopulations:

  Character vector of cell population names, or `"all"` (default) to
  return all populations.

- returnType:

  One of `"GRangesList"` (default) or `"data.frame"`.

## Value

A `GRangesList` or `data.frame` of open tiles per cell population.
