# Convert a list of strings in the format "chr1:100-200" into a GRanges

`StringsToGRanges` Turns a list of strings defining genomic regions in
the format chr1:100-200 into a GRanges object

## Usage

``` r
StringsToGRanges(regionString)
```

## Arguments

- regionString:

  A string or list of strings each in the format chr1:100-200

## Value

a GRanges object with ranges representing the input string(s)
