# Export normalized coverage files to per-sample or sample-averaged (per cell population) BigWig files.

`exportCoverage` will export normalized coverage files to BigWig files,
either as sample-specific or sample-averaged files, for visualization in
genome browsers.

## Usage

``` r
exportCoverage(
  SampleTileObject,
  dir = getwd(),
  cellPopulations = "ALL",
  groupColumn = NULL,
  subGroups = NULL,
  sampleSpecific = FALSE,
  saveFile = TRUE,
  numCores = 1,
  verbose = FALSE
)
```

## Arguments

- SampleTileObject:

  The SummarizedExperiment object output from getSampleTileMatrix

- dir:

  string. Directory to save files to.

- cellPopulations:

  vector of strings. Cell subsets for which to call peaks. This list of
  group names must be identical to names that appear in the
  SampleTileObject. Optional, if cellPopulations='ALL', then peak
  calling is done on all cell populations. Default is 'ALL'.

- groupColumn:

  Optional, the column containing sample group labels for returning
  coverage within sample groups. Default is NULL, all samples will be
  used.

- subGroups:

  a list of subgroup(s) within the groupColumn from the metadata.
  Optional, default is NULL, all labels within groupColumn will be used.

- sampleSpecific:

  If TRUE, a BigWig will export for each sample-cell type combination.

- saveFile:

  Boolean. If TRUE, it will save to a BigWig. If FALSE, it will return
  the GRangesList without writing a BigWig.

- numCores:

  integer. Number of cores to parallelize over

- verbose:

  Set TRUE to display additional messages. Default is FALSE.

## Value

A list of GRanges objects, or a list of paths to where the BigWigs are
saved

## Examples

``` r
if (FALSE) { # \dontrun{
MOCHA::exportCoverage(
  SampleTileObject = SampleTileMatrices,
  cellPopulations = "ALL",
  numCores = 30,
  sampleSpecific = FALSE
)
} # }
```
