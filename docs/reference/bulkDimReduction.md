# Run PCA or LSI dimensionality reduction on tiles

`bulkDimReduction` runs dimensionality reduction (either PCA or LSI)

## Usage

``` r
bulkDimReduction(
  SampleTileObj,
  cellType = "All",
  componentNumber = 30,
  method = "LSI",
  verbose = FALSE
)
```

## Arguments

- SampleTileObj:

  The SummarizedExperiment object output from getSampleTileMatrix

- cellType:

  vector of strings. Cell subsets for which to call peaks. This list of
  group names must be identical to names that appear in the
  SampleTileObj. Optional, if cellPopulations='ALL', then peak calling
  is done on all cell populations. Default is 'ALL'.

- componentNumber:

  integer. Number of components to include in LSI, or PCA This must be
  strictly less than

- method:

  a string. Represents the method to use. Includes LSI or PCA, but we do
  not recommend PCA for scATAC pseudobulk.

- verbose:

  Set TRUE to display additional messages. Default is FALSE.

## Value

SEObj a SummarizedExperiment containing PC components from
dimensionality reduction and metadata from the SampleTileObj

## References

LSI method adapted from Andrew Hill:
http://andrewjohnhill.com/blog/2019/05/06/dimensionality-reduction-for-scatac-data/

## Examples

``` r
if (FALSE) { # \dontrun{
LSIObj <- MOCHA::bulkDimReduction(SampleTileObj, cellType = "CD16_Mono")
} # }
```
