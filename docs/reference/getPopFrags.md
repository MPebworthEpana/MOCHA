# Extract fragments by populations from an ArchR Project

`getPopFrags` returns a list of sample-specific fragments per cell
population as a GRangesList.

## Usage

``` r
getPopFrags(
  ArchRProj,
  cellPopLabel,
  cellSubsets = "ALL",
  poolSamples = FALSE,
  numCores = 1,
  returnGRangesList = TRUE,
  verbose = FALSE
)
```

## Arguments

- ArchRProj:

  The ArchR Project.

- cellPopLabel:

  The name of the metadata column of the ArchR Project that contains the
  populations of cells you want to extract fragments from.

- cellSubsets:

  Default is 'ALL'. If you want to export only some populations, then
  give it a list of group names. This needs to be unique - no duplicated
  names. This list of group names must be identical to names that appear
  in the given cellPopLabel metadata column of the ArchR Project.

- poolSamples:

  Set TRUE to pool sample-specific fragments by cell population. By
  default this is FALSE and sample-specific fragments are returned.

- numCores:

  Number of cores to use.

- returnGRangesList:

  Set FALSE to return as a list of GRanges. Default, TRUE, returns a
  GRangesList

- verbose:

  Set TRUE to display additional messages. Default is FALSE.

## Value

A list of GRanges containing fragments. Each GRanges corresponds to a
population defined by cellSubsets and sample.
