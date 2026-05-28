# Get consensus sample-tile matrices containing the signal intensity at each tile

`getSampleTileMatrix` takes the output of peak calling with
callOpenTiles and creates sample-tile matrices containing the signal
intensity at each tile.

## Usage

``` r
getSampleTileMatrix(
  tileResults,
  cellPopulations = "ALL",
  groupColumn = NULL,
  reproducibilityThreshold = 0.2,
  numCores = 1,
  verbose = FALSE,
  returnClass = c("legacy", "mocha"),
  threshold = NULL
)
```

## Arguments

- tileResults:

  a MultiAssayExperiment returned by callOpenTiles containing containing
  peak calling results.

- cellPopulations:

  vector of strings. Cell subsets in TileResults for which to generate
  sample-tile matrices. This list of group names must be identical to
  names that appear in the ArchRProject metadata. If
  cellPopulations='ALL', then peak calling is done on all cell
  populations in the ArchR project metadata. Default is 'ALL'.

- groupColumn:

  Optional, the column containing sample group labels for determining
  consensus tiles within sample groups. Default is NULL, all samples
  will be used for determining consensus tiles.

- reproducibilityThreshold:

  Threshold for consensus tiles, the minimum % of samples (within a
  sample group, if groupColumn is set) that a peak must be called in to
  be retained. If set to 0, retain the union of all samples' peaks (this
  is equivalent to a threshold of 1/numSamples). Tune to omit
  potentially spurious peaks;
  [`suggestConsensusThreshold`](https://aifimmunology.github.io/MOCHA/reference/suggestConsensusThreshold.md)
  recommends a value automatically. Default 0.2.

- numCores:

  Optional, the number of cores to use with multiprocessing. Default is
  1.

- verbose:

  Set TRUE to display additional messages. Default is FALSE.

- returnClass:

  Return type: `"legacy"` (default, `RangedSummarizedExperiment`) or
  `"mocha"` (`MochaSampleTileMatrix`).

- threshold:

  Deprecated alias for `reproducibilityThreshold`. Use
  `reproducibilityThreshold` instead.

## Value

SampleTileMatrices a MultiAssayExperiment containing a sample-tile
intensity matrix for each cell population

## Examples

``` r
# \donttest{
# Starting from GRangesList
if (
  require(BSgenome.Hsapiens.UCSC.hg19) &&
    require(TxDb.Hsapiens.UCSC.hg38.knownGene) &&
    require(org.Hs.eg.db)
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

  SampleTileMatrices <- MOCHA::getSampleTileMatrix(
    tiles,
    cellPopulations = c("C2", "C5"),
    reproducibilityThreshold = 0 # Take union of all samples' open tiles
  )
}
# }
```
