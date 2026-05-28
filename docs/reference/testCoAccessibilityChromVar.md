# Test input tile pairs against a ChromVAR background

\`r lifecycle::badge("deprecated")\` This function is deprecated due as
we no longer recommend using a chromVAR background set.
`testCoAccessibilityChromVar` takes an input set of tile pairs and tests
whether they are significantly different compared to a background set
found via ChromVAR

## Usage

``` r
testCoAccessibilityChromVar(
  SampleTileObj,
  tile1,
  tile2,
  numCores = 1,
  ZI = TRUE,
  backNumber = 1000,
  returnBackGround = FALSE,
  highMem = FALSE,
  verbose = TRUE
)
```

## Arguments

- SampleTileObj:

  The SummarizedExperiment object output from getSampleTileMatrix
  containing your sample-tile matrices

- tile1:

  vector of indices or tile names (chrX:100-2000) for tile pairs to test
  (first tile in each pair)

- tile2:

  vector of indices or tile names (chrX:100-2000) for tile pairs to test
  (second tile in each pair)

- numCores:

  Optional, the number of cores to use with multiprocessing. Default is
  1.

- ZI:

  boolean flag that enables zero-inflated (ZI) Spearman correlations to
  be used. Default is TRUE. If FALSE, skip zero-inflation and calculate
  the normal Spearman.

- backNumber:

  number of ChromVAR-matched background pairs. Default is 1000.

- returnBackGround:

  Boolean, if TRUE return the background correlations as well as
  foreground. Default is FALSE.

- highMem:

  Boolean to control memory usage. Default is FALSE. Only set highMem to
  TRUE if you have plenty of memory and want to run this function
  faster.s

- verbose:

  Set TRUE to display additional messages. Default is FALSE.

## Value

foreGround A data.frame with Tile1, Tile2, Correlation, and p-value for
that correlation compared to the background
