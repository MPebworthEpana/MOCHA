# Get open tiles from a MochaTileResults object

Get open tiles from a MochaTileResults object

## Usage

``` r
openTiles(
  x,
  cellPopulations = "all",
  returnType = c("GRangesList", "data.frame")
)
```

## Arguments

- x:

  A `MochaTileResults` object.

- cellPopulations:

  Cell populations to return, or `"all"`.

- returnType:

  `"GRangesList"` or `"data.frame"`.

## Value

Open tiles per cell population.
