# Conduct a differential test between open regions of two sample groups

`getDifferentialAccessibleTiles` allows you to determine whether regions
of chromatin are differentially accessible between groups by conducting
a test

## Usage

``` r
getDifferentialAccessibleTiles(
  SampleTileObj,
  cellPopulations = NULL,
  groupColumn,
  foreground,
  background,
  minZeroDiff = 0.5,
  qValueMethod = "experimental",
  signalThreshold = NULL,
  qValueThreshold = 0.2,
  outputGRanges = TRUE,
  numCores = 1,
  verbose = FALSE,
  dropoutAdjustment = c("none", "biological_only", "weighted"),
  bioThreshold = 0.8,
  techThreshold = NULL,
  fdrToDisplay = NULL,
  method = c("wilcoxon", "paired_wilcoxon", "polr"),
  pairColumn = NULL,
  cellPopulation = NULL
)
```

## Arguments

- SampleTileObj:

  The SummarizedExperiment object output from getSampleTileMatrix

- cellPopulations:

  Character vector of cell-population names (matching assay names in
  `SampleTileObj`) to test. Length \\\ge 1\\.

- groupColumn:

  The column containing sample group labels

- foreground:

  The foreground group of samples for differential comparison

- background:

  The background group of samples for differential comparison

- minZeroDiff:

  Minimum difference in average dropout rates across groups require to
  keep tiles for differential testing. Default is 0.5 (50%).

- qValueMethod:

  String describing qvalue method. Can be 'standard', or 'experimental'.
  See methods in the MOCHA manuscript for description of the
  experimental. Otherwise, 'standard' applies standard q value.

- signalThreshold:

  Optional value. This is Minimum median intensity required to keep
  tiles for differential testing to increase statistical power in small
  sample cohorts. Default is NULL, at which point the optimal threshold
  will be found.

- qValueThreshold:

  A number greater than 0 and less than 1, used to optimize the noise
  thresholding parameter by cell type. Default is 0.2. It is recommended
  to keep this paramter between 0.05 and 0.3.

- outputGRanges:

  Outputs a GRanges if TRUE and a data.frame if FALSE. Default is TRUE.

- numCores:

  The number of cores to use with multiprocessing. Default is 1.

- verbose:

  Set TRUE to display additional messages. Default is FALSE.

- dropoutAdjustment:

  How to incorporate dropout estimates when computing zero rates for
  tile filtering. `"none"` uses raw zero fractions (default).
  `"biological_only"` counts only zeros with `P(zero) >= bioThreshold`
  (biologically expected closure). `"weighted"` weights zeros by
  `P(zero)`. Requires a fittable `DropoutProb_*` assay from
  [`assessDropout`](https://aifimmunology.github.io/MOCHA/reference/assessDropout.md);
  if absent,
  [`assessDropout()`](https://aifimmunology.github.io/MOCHA/reference/assessDropout.md)
  is run automatically for the requested cell population.

- bioThreshold:

  Minimum `P(zero)` for a zero to count toward biological zero rates
  when `dropoutAdjustment = "biological_only"`. Default is 0.8.

- techThreshold:

  Deprecated alias for `bioThreshold`.

- fdrToDisplay:

  Deprecated alias for `qValueThreshold`.

- method:

  Test method. One of `"wilcoxon"` (default unpaired two-part Wilcoxon),
  `"paired_wilcoxon"` (paired two-part test using the column named in
  `pairColumn`), or `"polr"` (proportional- odds cumulative-logit on
  ordinal-binned values; requires the MASS package).

- pairColumn:

  Column in `colData(SampleTileObj)` that identifies matched pairs
  across the foreground/background groups. Required when
  `method = "paired_wilcoxon"`; ignored otherwise.

- cellPopulation:

  Deprecated alias for `cellPopulations`. Use `cellPopulations` instead.

## Value

full_results The differential accessibility results as a GRanges or
matrix data.frame depending on the flag \`outputGRanges\`.

## Examples

``` r
# \donttest{
if (
  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE) &&
    requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)
) {
  tiles <- MOCHA::callOpenTiles(
    ATACFragments = MOCHA::exampleFragments,
    cellColData = MOCHA::exampleCellColData,
    blackList = MOCHA::exampleBlackList,
    genome = "BSgenome.Hsapiens.UCSC.hg19",
    TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    OrgDb = "org.Hs.eg.db",
    outDir = tempdir(),
    cellPopLabel = "Clusters",
    cellPopulations = "C2",
    numCores = 1
  )
  stm <- MOCHA::getSampleTileMatrix(
    tiles,
    cellPopulations = "C2",
    threshold = 0
  )
  diffs <- MOCHA::getDifferentialAccessibleTiles(
    SampleTileObj = stm,
    cellPopulations = "C2",
    groupColumn = "Sample",
    foreground = unique(SummarizedExperiment::colData(stm)$Sample)[1],
    background = unique(SummarizedExperiment::colData(stm)$Sample)[2],
    numCores = 1
  )
}
#> Error in MOCHA::getDifferentialAccessibleTiles(SampleTileObj = stm, cellPopulations = "C2",     groupColumn = "Sample", foreground = unique(SummarizedExperiment::colData(stm)$Sample)[1],     background = unique(SummarizedExperiment::colData(stm)$Sample)[2],     numCores = 1): Provided groupCol '{groupColumn}' not found in the provided SampleTileObj
# }
```
