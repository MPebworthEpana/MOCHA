# Subset a tileResults object by metadata

`subsetMOCHAObject` subsets a tileResults-type object (from
callOpenTiles), or a SummarizedExperiment-type object (from
getSampleTileMatrix), either by cell type or sample metadata.

## Usage

``` r
subsetMOCHAObject(
  Object,
  subsetBy,
  groupList,
  removeNA = TRUE,
  subsetPeaks = TRUE,
  verbose = FALSE
)
```

## Arguments

- Object:

  A MultiAssayExperiment or RangedSummarizedExperiment,

- subsetBy:

  The variable to subset by. Can either be 'celltype', or a column from
  the sample metadata (see \`colData(Object)\`).

- groupList:

  the list of cell type names or sample-associated data that should be
  used to subset the Object

- removeNA:

  If TRUE, removes groups in groupList that are NA. If FALSE, keep
  groups that are NA.

- subsetPeaks:

  If \`subsetBy\` = 'celltype', subset the tile set to tiles only called
  in those cell types. Default is TRUE.

- verbose:

  Set TRUE to display additional messages. Default is FALSE.

## Value

Object the input Object, filtered down to either the cell type or
samples desired.

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
  subsetTiles <- MOCHA::subsetMOCHAObject(
    tiles,
    subsetBy = "celltype",
    groupList = "C2"
  )
}
#> Warning: 'experiments' dropped; see 'drops()'
#> harmonizing input:
#>   removing 1 sampleMap rows not in names(experiments)
# }
```
