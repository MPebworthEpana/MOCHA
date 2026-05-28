# Extract accessibility coverage for a given region

`extractRegion` will extract the coverage files created by callOpenTiles
and return a specific region's coverage

## Usage

``` r
extractRegion(
  SampleTileObj,
  type = TRUE,
  region,
  cellPopulations = "ALL",
  groupColumn = NULL,
  subGroups = NULL,
  sampleSpecific = FALSE,
  approxLimit = 1e+05,
  skipEmpty = TRUE,
  binSize = 250,
  sliding = NULL,
  numCores = 1,
  verbose = FALSE
)
```

## Arguments

- SampleTileObj:

  The SummarizedExperiment object output from getSampleTileMatrix

- type:

  Boolean. Default is true, and exports Coverage. If set to FALSE,
  exports Insertions.

- region:

  a GRanges object or vector or strings containing the regions of
  interest. Strings must be in the format "chr:start-end", e.g.
  "chr4:1300-2222".

- cellPopulations:

  vector of strings. Cell subsets for which to call peaks. This list of
  group names must be identical to names that appear in the
  SampleTileObj. Optional, if cellPopulations='ALL', then peak calling
  is done on all cell populations. Default is 'ALL'.

- groupColumn:

  Optional, the column containing sample group labels for returning
  coverage within sample groups. Default is NULL, all samples will be
  used.

- subGroups:

  a list of subgroup(s) within the groupColumn from the metadata.
  Optional, default is NULL, all labels within groupColumn will be used.

- sampleSpecific:

  If TRUE, get a sample-specific count dataframe out. Default is FALSE,
  average across samples and get a dataframe out.

- approxLimit:

  Optional limit to region size, where if region is larger than
  approxLimit basepairs, binning will be used. Default is 100000.

- binSize:

  Optional numeric, size of bins in basepairs when binning is used.
  Default is 250.

- sliding:

  Optional numeric. Default is NULL. This number is the size of the
  sliding window for generating average intensities.

- numCores:

  integer. Number of cores to parallelize peak-calling across multiple
  cell populations

- verbose:

  Set TRUE to display additional messages. Default is FALSE.

## Value

countSE a SummarizedExperiment containing coverage for the given input
cell populations.

## Examples

``` r
if (FALSE) { # \dontrun{
countSE <- MOCHA::extractRegion(
  SampleTileObj = SampleTileMatrices,
  cellPopulations = "ALL",
  region = "chr1:18137866-38139912",
  numCores = 30,
  sampleSpecific = FALSE
)
} # }
```
