# Retrain MOCHA peak-calling models at a chosen tile size

Fits cell-count-specific logistic regression models and Youden-optimal
probability thresholds following the MOCHA publication Methods (NK-cell
training on MACS2 pseudo-bulk labels). The returned `MOCHAPeakModel` can
be passed to
[`callOpenTiles`](https://aifimmunology.github.io/MOCHA/reference/callOpenTiles-methods.md)
via `peakModel`.

## Usage

``` r
trainPeakModel(
  ATACFragments,
  cellColData,
  blackList,
  groundTruthPeaks = NULL,
  tileSize = .MOCHA_DEFAULT_TILE_SIZE,
  cellCol = "RG",
  sampleColumn = "Sample",
  cellSubsetSizes = c(5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000,
    50000, 1e+05, 170000),
  replicatesFn = function(n) {
     if (n < 50000) 
         10L
     else 5L
 },
  threshMethod = c("youden", "f1", "topleft"),
  smoothMethod = c("auto", "loess", "linear"),
  callPeaks = FALSE,
  macs2Args = list(g = "hs", f = "BED", nolambda = TRUE, shift = -75, extsize = 150,
    broad = TRUE),
  trimPeakBp = 75L,
  numCores = 1L,
  seed = 1L,
  verbose = FALSE
)
```

## Arguments

- ATACFragments:

  A `GRangesList` of single-cell ATAC fragments for one training cell
  population (typically pooled across samples).

- cellColData:

  Cell metadata; used to compute the training median fragments per cell
  when column `nFrags` is present.

- blackList:

  `GRanges` of regions to exclude when tiling.

- groundTruthPeaks:

  `GRanges` of accessible regions for training labels. Required unless
  `callPeaks = TRUE`.

- tileSize:

  Width of genomic tiles in base pairs. Default `500L`.

- cellCol:

  Column in fragment `mcols` with cell barcodes. Default `"RG"`.

- sampleColumn:

  Unused for training labels; reserved for compatibility.

- cellSubsetSizes:

  Integer vector of cell counts at which to fit models.

- replicatesFn:

  Function of cell count `n` returning the number of training
  replicates. Default: 10 below 50k cells, 5 otherwise.

- threshMethod:

  Method for optimal probability cutoffs per cell count: `"youden"`
  (requires cutpointr), `"f1"`, or `"topleft"`.

- smoothMethod:

  Coefficient smoothing: `"auto"` (LOESS + linear, as in
  `make_prediction`), `"loess"`, or `"linear"`.

- callPeaks:

  If `TRUE`, run MACS2 on pooled pseudobulk fragments to derive
  `groundTruthPeaks`. Requires `macs2` on `PATH` and rtracklayer.

- macs2Args:

  Named list of MACS2 arguments passed to `macs2 callpeak` (defaults
  mirror the MOCHA paper / ArchR settings).

- trimPeakBp:

  Base pairs trimmed from each end of broad peaks before overlay on
  tiles. Default `75L`.

- numCores:

  Cores for parallel replicate fitting. Default `1L`.

- seed:

  Random seed for cell subsampling.

- verbose:

  Print progress messages.

## Value

A `MOCHAPeakModel` object for use with `peakModel` in
[`callOpenTiles`](https://aifimmunology.github.io/MOCHA/reference/callOpenTiles-methods.md).

## Details

Optional dependencies:

- cutpointr when `threshMethod = "youden"`

- rtracklayer and a `macs2` binary when `callPeaks = TRUE`
