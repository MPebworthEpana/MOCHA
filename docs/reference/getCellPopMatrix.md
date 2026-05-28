# Get the SampleTileMatrix of the given cell population

`getCellPopMatrix` pulls out the SampleTileMatrix of tiles called in one
given cell population.

## Usage

``` r
getCellPopMatrix(
  SampleTileObj,
  cellPopulation,
  dropSamples = TRUE,
  NAtoZero = TRUE
)
```

## Arguments

- SampleTileObj:

  The output from getSampleTileMatrix, a SummarizedExperiment of
  pseudobulk intensities across all tiles & cell types.

- cellPopulation:

  The cell population you want to pull out.

- dropSamples:

  Boolean flag to determine whether to drop samples that were too small
  for peak calling.

- NAtoZero:

  Boolean flag to determine whether to replace NAs with zero

## Value

sampleTileMatrix a matrix of samples by called tiles for a given cell
population.
