# `singlePopulationSampleTileMatrix`

`singlePopulationSampleTileMatrix` is a function that can transform a
set of tile intensities into peak X sample matrix for a custom set of
tiles

## Usage

``` r
singlePopulationSampleTileMatrix(
  peaksExperiment,
  consensusTiles,
  NAtoZero = FALSE
)
```

## Arguments

- peaksExperiment:

  peakset RaggedExperiment, one celltype from the output of
  callOpenTiles

- consensusTiles:

  a vector containing the tileIDs to subset the sample-tile matrix

## Value

sampleTileIntensityMat a sample X peak matrix containing observed
measurements for each sample at each peak.

## Details

The technical details of the algorithm are found in XX.

## References

XX
