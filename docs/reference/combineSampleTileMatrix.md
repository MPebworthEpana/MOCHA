# Merge the TSAM from multiple cell populations into a single matrix

`combineSampleTileMatrix` combines all celltypes in a SampleTileMatrix
object into a SummarizedExperiment with one single matrix across all
cell types and samples,

## Usage

``` r
combineSampleTileMatrix(
  SampleTileObj,
  NAtoZero = TRUE,
  verbose = FALSE,
  returnClass = c("legacy", "mocha")
)
```

## Arguments

- SampleTileObj:

  The SummarizedExperiment object output from getSampleTileMatrix
  containing your sample-tile matrices

- NAtoZero:

  Set NA values in the sample-tile matrix to zero

- verbose:

  Set TRUE to display additional messages. Default is FALSE.

- returnClass:

  Return type: `"legacy"` (default) or `"mocha"`
  (`MochaSampleTileMatrix` when the input is a MOCHA sample-tile
  matrix).

## Value

A `RangedSummarizedExperiment` with one matrix across cell types.

## Examples

``` r
# \donttest{
if (
  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE) &&
    requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)
) {
  tiles <- MOCHA::callOpenTiles(
    ATACFragments = MOCHA::exampleFragments,
    cellColData = MOCHA::exampleCellColData,
    blackList = MOCHA::exampleBlackList,
    genome = "BSgenome.Hsapiens.UCSC.hg19",
    TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    OrgDb = "org.Hs.eg.db",
    outDir = tempdir(),
    cellPopLabel = "Clusters",
    cellPopulations = c("C2", "C5"),
    numCores = 1
  )
  stm <- MOCHA::getSampleTileMatrix(
    tiles,
    cellPopulations = c("C2", "C5"),
    threshold = 0
  )
  combined <- MOCHA::combineSampleTileMatrix(stm)
}
#> Error in MOCHA::combineSampleTileMatrix(stm): Could not map all combined sample keys to CellCounts metadata. Missing: C2__, C5__
# }
```
