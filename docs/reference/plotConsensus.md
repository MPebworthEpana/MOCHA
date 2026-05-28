# Plot to determine the reproducibility threshold

`plotConsensus` Extracts the peak reproducibility and generates a
heuristic plots that can be used to determine the reproducibility
threshold used within getSampleTileMatrix.

## Usage

``` r
plotConsensus(
  tileObject,
  cellPopulations = "All",
  groupColumn = NULL,
  returnPlotList = FALSE,
  returnDFs = FALSE,
  showSuggested = FALSE,
  suggestMethod = c("kneedle", "second_derivative"),
  numCores = 1
)
```

## Arguments

- tileObject:

  A MultiAssayExperiment object from callOpenTiles,

- cellPopulations:

  the cell populations you want to visualize.

- groupColumn:

  Optional parameter, same as in getSampleTileMatrix, which defines
  whether you want to plot reproducibility within each

- returnPlotList:

  Instead of one plot with all celltypes/conditions, it returns a list
  of plots for each cell types

- returnDFs:

  Instead of a plot, returns a data.frame of the reproducibility across
  samples. If set to false, then it plots the data.frame instead of
  returning it.

- showSuggested:

  Logical. If `TRUE`, overlay a dashed vertical line at the threshold
  recommended by
  [`suggestConsensusThreshold()`](https://aifimmunology.github.io/MOCHA/reference/suggestConsensusThreshold.md)
  for each cell population. Default `FALSE`.

- suggestMethod:

  Method passed to
  [`suggestConsensusThreshold()`](https://aifimmunology.github.io/MOCHA/reference/suggestConsensusThreshold.md)
  when `showSuggested = TRUE`. One of `"kneedle"` (default) or
  `"second_derivative"`.

- numCores:

  Number of cores to multithread over.

## Value

A data.frame of reproducibility, or plots
