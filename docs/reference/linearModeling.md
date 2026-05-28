# `linearModeling`

`linearModeling` Runs linearModeling on MOCHA TSAM

## Usage

``` r
linearModeling(
  Obj,
  formula,
  CellType,
  threshold = 0,
  NAtoZero = FALSE,
  numCores = 1
)
```

## Arguments

- Obj:

  A RangedSummarizedExperment generated from getSampleTileMatrix

- formula:

  Formula used for linear modeling.

- CellType:

  The name of the celltype you wish to model. Should align with
  assayNames of the Obj.

- threshold:

  A number greater than 0 and less then or equal to 1. The threshold
  used to determine whether to model a region or not, based on the
  fraction of non-zero measurements across samples at that location.

- NAtoZero:

  Boolean, whether to convert NA region (no accessibility measurement)
  to zero.

- numCores:

  Number of threads to parallelize modeling over. Default is 1.

## Value

A list of lmer model objects
