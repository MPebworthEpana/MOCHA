# Convert a data.frame or matrix to a GRanges

Convert a data.frame or matrix to a GRanges

## Usage

``` r
differentialsToGRanges(differentials, tileColumn = "Tile")
```

## Arguments

- differentials:

  a matrix/data.frame with a column tileColumn containing region strings
  in the format "chr:start-end"

- tileColumn:

  name of column containing region strings. Default is "Tile".

## Value

a GRanges containing all original information
