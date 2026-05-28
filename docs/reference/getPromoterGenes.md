# Extract the list of promoter genes from a GRanges annotated with `annotateTiles()`

`getPromoterGenes` Takes rowRanges from annotateTiles and extracts a
unique list of genes.

## Usage

``` r
getPromoterGenes(GRangesObj)
```

## Arguments

- GRangesObj:

  a GRanges object with a metadata column for tileType and Gene.

## Value

vector of strings with gene names.
