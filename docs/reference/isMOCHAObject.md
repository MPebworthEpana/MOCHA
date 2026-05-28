# isMOCHAObject

`isMOCHAObject` Checks with the object is a MOCHA object, and/or returns
the type of object

## Usage

``` r
isMOCHAObject(Object, returnType = FALSE)
```

## Arguments

- Object:

  The output of callOpenTiles or getSampleTileMatrix. Both have a path
  to saved files.

- returnType:

  A Boolean. If True, returns object type (OpenTiles, SampleTileMatrix,
  or NULL). If FALSE, returns True/FALSE

## Value

A string or list of strings in the format 'chr1:100-200' representing
ranges in the input GRanges
