# Export motif annotations from `addMotifSet()` to BigBed format for visualization in genome browsers.

`exportMotifs` exports a motif set GRanges from running
[`addMotifSet()`](https://aifimmunology.github.io/MOCHA/reference/addMotifSet.md)
(with \`returnSTM=FALSE\`) to bigBed file files for visualization in
genome browsers.

## Usage

``` r
exportMotifs(
  SampleTileObject,
  motifsGRanges,
  motifSetName = "motifs",
  filterByOpenTiles = FALSE,
  outDir,
  verbose = FALSE
)
```

## Arguments

- SampleTileObject:

  The SummarizedExperiment object output from
  [`getSampleTileMatrix()`](https://aifimmunology.github.io/MOCHA/reference/getSampleTileMatrix.md)

- motifsGRanges:

  A GRanges containing motif annotations, typically from
  [`addMotifSet()`](https://aifimmunology.github.io/MOCHA/reference/addMotifSet.md)
  (with \`returnSTM=FALSE\`)

- motifSetName:

  Optional, a name indicating the motif set. Used to name files in the
  specified `outdir`. Default is "motifs".

- filterByOpenTiles:

  Boolean. If TRUE, a bigBed file will be exported for each cell
  population with motifs filtered to those occurring only in open tiles.

- outDir:

  Desired output directory where bigBed files will be saved

- verbose:

  Set TRUE to display additional messages. Default is FALSE.

## Value

outList A List of output filepaths

## Examples

``` r
if (FALSE) { # \dontrun{
MOCHA::exportMotifs(
  SampleTileObject = SampleTileMatrices,
  motifsGRanges,
  motifSetName = "CISBP",
  filterByOpenTiles = FALSE,
  outDir = tempdir(),
  verbose = TRUE
)
} # }
```
