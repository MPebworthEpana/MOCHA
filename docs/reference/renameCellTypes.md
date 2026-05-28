# Modify the cell population names in a Sample-Tile Object from `getSampleTileMatrix()`

`renameCellTypes` Allows you to modify the cell type names for a MOCHA
SampleTileObject, from the assay names, GRanges column names, and
summarizedData (within the metadata), all at once.

## Usage

``` r
renameCellTypes(MOCHAObject, oldNames, newNames)
```

## Arguments

- MOCHAObject:

  A RangedSummarizedExperiment,

- oldNames:

  A list of cell type names that you want to change.

- newNames:

  A list of new cell type names to replace the old names with.

## Value

A MOCHA SampleTile object with new cell types.
