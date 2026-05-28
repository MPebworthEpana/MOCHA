# Extract Sample-celltype specific metadata

`getSampleCellTypeMetadata` Extract Sample-celltype specific metadata
like fragment and cell counts

## Usage

``` r
getSampleCellTypeMetadata(object)
```

## Arguments

- object:

  tileResults object from callOpenTiles or SummarizedExperiment from
  getSampleTileMatrix

## Value

a SummarizedExperiment where each assay is a different type of metadata.
