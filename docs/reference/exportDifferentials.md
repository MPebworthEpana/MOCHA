# Export differential peaks from `getDifferentialAccessibleTiles()` to BigBed format for visualization in genome browsers.

`exportDifferentials` exports the differential peaks output GRangesList
output from
[`getDifferentialAccessibleTiles()`](https://aifimmunology.github.io/MOCHA/reference/getDifferentialAccessibleTiles.md)
to bigBed format for visualization in genome browsers.

## Usage

``` r
exportDifferentials(
  SampleTileObject,
  DifferentialsGRList,
  outDir,
  verbose = FALSE
)
```

## Arguments

- SampleTileObject:

  The SummarizedExperiment object output from
  [`getSampleTileMatrix()`](https://aifimmunology.github.io/MOCHA/reference/getSampleTileMatrix.md)

- DifferentialsGRList:

  GRangesList output from
  [`getDifferentialAccessibleTiles()`](https://aifimmunology.github.io/MOCHA/reference/getDifferentialAccessibleTiles.md)

- outDir:

  Desired output directory where bigBed files will be saved

- verbose:

  Set TRUE to display additional messages. Default is FALSE.

## Value

outList A List of output filepaths

## Examples

``` r
if (FALSE) { # \dontrun{
MOCHA::exportDifferentials(
  SampleTileObject = SampleTileMatrices,
  DifferentialsGRList,
  outDir = tempdir(),
  verbose = TRUE
)
} # }
```
