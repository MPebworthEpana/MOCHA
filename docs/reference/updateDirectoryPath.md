# updateDirectoryPath

`updateDirectoryPath` Updated the path to save files, in case the
directory moved.

## Usage

``` r
updateDirectoryPath(Object, directoryPath = NULL)
```

## Arguments

- Object:

  The output of callOpenTiles or getSampleTileMatrix. Both have a path
  to saved files.

- directoryPath:

  A string, containing the absolute path to the directory.

## Value

A string or list of strings in the format 'chr1:100-200' representing
ranges in the input GRanges
