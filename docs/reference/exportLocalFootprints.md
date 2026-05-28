# Export insertion counts to per-sample BigWig files after applying a rolling sum and rolling median smoothing filter.

`exportLocalFootprints` Takes a SampleTileMatrix with linked insertion
files, corrects Tn5 insertion bias, and applies a smoothing filter (a
rolling sum then rolling median) to the insertions. This generates a
local footprint track for visually identifying evidence of potential
binding sites. These files are then written to bigwig format.

## Usage

``` r
exportLocalFootprints(
  SampleTileObj,
  cellPopulation,
  outDir = NULL,
  windowSize = 10,
  groupColumn = NULL,
  subGroups = NULL,
  sampleSpecific = FALSE,
  normTn5 = TRUE,
  force = FALSE,
  slow = FALSE,
  verbose = FALSE,
  numCores = 5
)
```

## Arguments

- SampleTileObj:

  A MultiAssayExperiment or RangedSummarizedExperiment from MOCHA

- cellPopulation:

  A string denoting the cell population of interest

- outDir:

  Directory to write output bigwig files. Default is NULL, where the
  directory in \`SampleTileObj@metadata\$Directory\` will be used.

- windowSize:

  Window size for rolling sum & median over basepairs. Default is 10.

- groupColumn:

  A string, corresponding to a metadata column within the SampleTileObj,
  that describes the groups by which you want to summarize motif
  footprints. If sampleSpecific = FALSE, then the function will average
  insertions across samples within each group.

- subGroups:

  A list of subgroups, if you want to only look at specific groups
  within the groupColumn.

- sampleSpecific:

  A boolean for whether to generate average motif footprints within each
  group, or to return a data.frame for all samples.

- normTn5:

  A boolean for whether to normalize by Tn5 insertion bias.

- force:

  Set TRUE to overwrite existing files. Default is FALSE.

- slow:

  Set TRUE to bypass optimisations and compute smoothing filter directly
  on the whole genome. May run slower and consume more RAM. Default is
  FALSE.

- verbose:

  Set TRUE to display additional messages. Default is FALSE.

- numCores:

  Number of cores to parallelize over

## Value

outPaths List of paths of exported insertion files

## Examples

``` r
if (FALSE) { # \dontrun{
# Depends on and manipulates files on filesystem
outPath <- MOCHA::exportSmoothedInsertions(
  SampleTileObj,
  cellPopulation = "CD4 Naive", sumWidth = 10, medianWidth = 11, verbose = FALSE
)
} # }
```
