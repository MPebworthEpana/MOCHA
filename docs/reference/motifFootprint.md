# Generate motif footprints

Generate a plot of average normalized insertions around a motif center

## Usage

``` r
motifFootprint(
  SampleTileObj,
  motifName = "Motifs",
  specMotif = NULL,
  regions = NULL,
  cellPopulations = "ALL",
  windowSize = 500,
  normTn5 = TRUE,
  smoothTn5 = 10,
  groupColumn = NULL,
  subGroups = NULL,
  sampleSpecific = FALSE,
  numCores = 1,
  force = FALSE,
  verbose = FALSE
)
```

## Arguments

- SampleTileObj:

  A Sample-tile object from MOCHA's getSampleTileMatrix()

- motifName:

  The name of metadata entry with Motif location, added via
  addMotifSet()

- specMotif:

  An optional string specifying which motif to analyze. If blank, it
  will analyze all motifs.

- regions:

  An optional GRanges object or list of strings in the format
  chr1:100-200, specifiying which specific regions to look at when
  conducting motif footprinting.

- cellPopulations:

  A list of cell populations to conduct motif footprinting on.

- windowSize:

  A number, representing the window to analyze around each motif
  location. Default is 500 bp.

- normTn5:

  A boolean for whether to normalize by Tn5 insertion bias.

- smoothTn5:

  The window size for smoothen Tn5 insertions at each location. Ideal
  when looking at motif footprints over a smaller number of regions,
  rarer cell types, or sparse regions, where local noise can make it
  harder to see the overall pattern. Can be set to 0.

- groupColumn:

  A string, corresponding to a metadata column within the SampleTileObj,
  that describes the groups by which you want to summarize motif
  footprints. If sampleSpecific = FALSE, then motifFootprint will
  average insertions across samples within each group.

- subGroups:

  A list of subgroups, if you want to only look at specific groups
  within the groupColumn.

- sampleSpecific:

  A boolean for whether to generate average motif footprints within each
  group, or to return a data.frame for all samples.

- numCores:

  Number of cores to parallelize over

- force:

  Boolean. If FALSE, it will through an error if there's an empty sample
  or not overlap with regions and motifs. If TRUE will ignore these
  issues and continue (or return NULL)

- verbose:

  Boolean. Default is FALSE. Will print more messages if TRUE.

## Value

A SummarizedExperiment containing motif footprinting data
