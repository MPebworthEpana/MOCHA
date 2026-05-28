# Test for enrichment of motifs against a background

Test for enrichment of motifs within Group1 against a background Group2
using a hypergeometric t-test.

## Usage

``` r
MotifEnrichment(Group1, Group2, motifPosList, type = NULL)
```

## Arguments

- Group1:

  A GRanges object, such as a set of significant differential tiles.

- Group2:

  A GRanges object containing background regions, non-overlapping with
  Group1

- motifPosList:

  A GRangesList of motifs and positions for each motif. Must be named
  for each motif.

- type:

  Optional, name of a metadata column in Group1 and Group2 to test for
  enrichment the number of unique entries in column given by 'type'.
  Default is NULL, which tests the number of Ranges.

## Value

A data.frame containing enrichment for each group
