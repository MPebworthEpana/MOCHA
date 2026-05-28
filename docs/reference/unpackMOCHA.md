# Unzip a MOCHA object and its linked files created by `packMOCHA()` for portability between file systems

`unpackMOCHA` will unpack a zip archive created by unpackMOCHA, setting
the stored MOCHA object's stored directory path to the new location. See
also:
[packMOCHA](https://aifimmunology.github.io/MOCHA/reference/packMOCHA.md)

## Usage

``` r
unpackMOCHA(zipfile, exdir, verbose = FALSE)
```

## Arguments

- zipfile:

  Filepath to the packed MOCHA object.

- exdir:

  The path to the external directory where you want to unpack the MOCHA
  object.

- verbose:

  Display additional messages. Default is FALSE.

## Value

MOCHAObj the MOCHA object (tileResults or Sample-Tile Matrix)

## Examples

``` r
if (FALSE) { # \dontrun{
# Depends on files existing on your system
MOCHA::unpackMOCHA(zipfile = "./mochaobj.zip", exdir = "./newMOCHAdir")
} # }
```
