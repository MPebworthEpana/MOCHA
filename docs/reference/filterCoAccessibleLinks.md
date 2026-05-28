# Filter links by correlation strength

`filterCoAccessibleLinks` will filter the output from
getCoAccessibleLinks by a threshold, retaining links with a absolute
correlation greater than the threshold. This function also adds the chr,
start, and end site of each link to the output table.

## Usage

``` r
filterCoAccessibleLinks(TileCorr, threshold = 0.5)
```

## Arguments

- TileCorr:

  The correlation table output from getCoAccessibleLinks

- threshold:

  Keep

## Value

FilteredTileCorr The filtered correlation table with chr, start, and end
site of each link

## Examples

``` r
if (FALSE) { # \dontrun{
# links is the output of MOCHA::getCoAccessibleLinks
MOCHA::filterCoAccessibleLinks(links, threshold = 0.5)
} # }
```
