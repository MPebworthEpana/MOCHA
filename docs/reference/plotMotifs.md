# plotMotifs - plots motif footprints, exported from motifFootprint

`plotMotifs` Generates plot of motif footprints

## Usage

``` r
plotMotifs(
  motifSE,
  footprint = NULL,
  groupColumn = NULL,
  returnDF = FALSE,
  returnPlotList = FALSE,
  topPercentage = 1,
  plotIndividualRegions = TRUE,
  relHeights = c(0.3, 0.7)
)
```

## Arguments

- motifSE:

  A SummarizedExperiment with motif footprinting information, from
  motifFootprint

- footprint:

  Optional string, to describe which footprint within the motifSE to
  analyze. Default is NULL, at which point it will pull from all
  footprint present.

- groupColumn:

  An optional string, that will contain the group-level labels for
  either calling footprints as significant, or for comparing motif
  footprintings stats betweem groups. Default is null, at which point no
  wilcoxon tests will be conducted. groupColumn can be a sample-level
  group from the MOCHA's object colData slot, 'CellType' (if you want to
  compare across cell types), or 'Motif' if you want to compare across
  motifs.

- returnDF:

  Boolean, default is FALSE, determines whether or not to return a
  data.frame, rather than plotting.

- returnPlotList:

  A boolean, default is false, determines whether to return the full
  plot, or a list of subplots (ggplot2-based) for custom arrangements.

- topPercentage:

  A number 1 or less, that describes the top percent of motif locations
  to use for plotting. If 0.9, then the top 90 percent of regions, by
  total insertions, will be used in the plot. If 0.1, then the top 10
  percent of regions will be used. Default of 1 uses all regions.

- plotIndividualRegions:

  A Boolean that determines whether to plot the individual motif
  regions. If FALSE, only the motif footprint will be returned. Default
  is TRUE.

- relHeights:

  A vector of two numbers, describing the relative space in the plot
  given to the motif summary and individual location plot.

## Value

a data.frame, containing motif footprint stats.

## Examples

``` r
if (FALSE) { # \dontrun{
p1 <- MOCHA::plotMotifs(
  motifSE,
  footprint = "CD4_Naive_ARID5A", groupColumn = 'COVID_status', topPercentages = 0.1)
} # }
```
