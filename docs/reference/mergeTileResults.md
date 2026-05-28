# Merge tileResults from `callOpenTiles()` that share cell populations into a single object containing all samples

`mergeTileResults` merges a list of tileResults that each contain unique
samples into a single object encompassing all samples. Only cell
populations shared among all input tileResults will be retained. This
function can merge MultiAssayExperiment objects from callOpenTiles that
are created with the same TxDb, OrgDb, and Genome assembly.

## Usage

``` r
mergeTileResults(
  tileResultsList,
  numCores = 1,
  verbose = TRUE,
  returnClass = c("legacy", "mocha")
)
```

## Arguments

- tileResultsList:

  List of MultiAssayExperiments objects returned by callOpenTiles
  containing containing peak calling results.

- numCores:

  Optional, the number of cores to use with multiprocessing. Default is
  1.

- verbose:

  Set TRUE to display additional messages. Default is FALSE.

- returnClass:

  Return type: `"legacy"` (default, `MultiAssayExperiment`) or `"mocha"`
  (`MochaTileResults`).

## Value

tileResults a single MultiAssayExperiment containing a sample-tile
intensity matrix for each sample and common cell population in the input
tileResultsList.

## Examples

``` r
if (FALSE) { # \dontrun{
# Depends on local MOCHA tileResults
MOCHA::mergeTileResults(
  list(tileResultsCelltypesABC, tileResultsCelltypesBCD)
)
} # }
```
