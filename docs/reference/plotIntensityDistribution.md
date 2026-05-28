# Plots the distribution of sample-tile intensities for a given cell population

`plotIntensityDistribution` Plots the distribution of sample-tile
intensities for a give cell population.

## Usage

``` r
plotIntensityDistribution(
  TSAM_object,
  cellPopulation,
  returnDF = FALSE,
  density = TRUE
)
```

## Arguments

- TSAM_object:

  SummarizedExperiment from getSampleTileMatrix

- cellPopulation:

  Cell type names (assay name) within the TSAM_object

- returnDF:

  If TRUE, return the data frame without plotting. Default is FALSE.

- density:

  Boolean to determine whether to plot density or histogram. Default is
  TRUE (plots density).

## Value

data.frame or ggplot histogram.
