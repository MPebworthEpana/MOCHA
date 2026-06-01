# Zip up a MOCHA object and its linked files for portability between file systems

`packMOCHA` combines a MOCHA object (Sample-Tile Matrix or tileResults)
with its saved coverage tracks into a single zip archive. This allows
MOCHA objects and the necessary coverage files for plotting to be shared
to other file systems. See also:
[unpackMOCHA](https://aifimmunology.github.io/MOCHA/reference/unpackMOCHA.md)

## Usage

``` r
packMOCHA(MOCHAObj, zipfile, verbose = FALSE)
```

## Arguments

- MOCHAObj:

  A MultiAssayExperiment or RangedSummarizedExperiment, from MOCHA

- zipfile:

  Filename and path of the zip archive.

- verbose:

  Set TRUE to display additional messages. Default is FALSE.

## Value

zipfile Path to zip archive.

## Examples

``` r
if (FALSE) { # \dontrun{
# Depends on and manipulates files on filesystem
myOutputDir <- "/home/documents/MOCHA_out"
zipPath <- MOCHA::packMOCHA(
  tileResults, zipfile = file.path(myOutputDir, "testzip.zip")
)
} # }
```
