# Export open tiles of a given cell population to BigBed format for visualization in genome browsers.

`exportOpenTiles` exports the open tiles of a given cell population to
bigBed file for visualization in genome browsers.

## Usage

``` r
exportOpenTiles(SampleTileObject, cellPopulation, outDir, verbose = FALSE)
```

## Arguments

- SampleTileObject:

  The SummarizedExperiment object output from
  [`getSampleTileMatrix()`](https://aifimmunology.github.io/MOCHA/reference/getSampleTileMatrix.md)

- cellPopulation:

  The name of the cell population to export

- outDir:

  Desired output directory where bigBed files will be saved

- verbose:

  Set TRUE to display additional messages. Default is FALSE.

## Value

outList A List of output filepaths

## Examples

``` r
if (FALSE) { # \dontrun{
MOCHA::exportOpenTiles(
  SampleTileObject = SampleTileObject,
  cellPopulation,
  outDir = tempdir(),
  verbose = TRUE
)
} # }
```
