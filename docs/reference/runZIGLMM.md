# Run Zero-inflated Generalized Linear Mixed Modeling on pseudobulked scATAC data

\`r lifecycle::badge("deprecated")\` This function is deprecated -
improved modeling functions can be found in the package "ChAI" at
https://github.com/aifimmunology/ChAI `runZIGLMM` Runs linear
mixed-effects modeling for zero-inflated data using
[`glmmTMB`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html).

## Usage

``` r
runZIGLMM(
  TSAM_Object,
  cellPopulation = "all",
  continuousFormula = NULL,
  ziformula = NULL,
  zi_threshold = 0,
  initialSampling = 5,
  verbose = FALSE,
  numCores = 1
)
```

## Arguments

- TSAM_Object:

  A SummarizedExperiment object generated from getSampleTileMatrix.

- cellPopulation:

  Name of a cell type(s), or 'all'. The function will combine the cell
  types mentioned into one matrix before running the model.

- continuousFormula:

  The formula for the continuous data that should be used within
  glmmTMB. It should be in the format (exp ~ factors). All factors must
  be found in column names of the TSAM_Object metadata, except for
  CellType, FragNumber and CellCount, which will be extracted from the
  TSAM_Object. modelFormula must start with 'exp' as the response. See
  [glmmTMB](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html).

- ziformula:

  The formula for the zero-inflated data that should be used within
  glmmTMB. It should be in the format ( ~ factors). All factors must be
  found in column names of the TSAM_Object colData metadata, except for
  CellType, FragNumber and CellCount, which will be extracted from the
  TSAM_Object.

- zi_threshold:

  Zero-inflated threshold ( range = 0-1), representing the fraction of
  samples with zeros. When the percentage of zeros in the tile is
  between 0 and zi_threshold, samples with zeroes are dropped and only
  the continous formula is used. Use this parameter at your own risk.
  Default is 0.

- initialSampling:

  Size of data to use for pilot

- verbose:

  Set TRUE to display additional messages. Default is FALSE.

- numCores:

  integer. Number of cores to parallelize across.

## Value

results a SummarizedExperiment containing LMEM results

## Examples

``` r
if (FALSE) { # \dontrun{
modelList <- runZIGLMM(STM[c(1:1000), ],
  cellPopulation = "CD16 Mono",
  continuousFormula = exp ~ Age + Sex + days_since_symptoms + (1 | PTID),
  ziformula = ~ FragNumber + Age,
  verbose = TRUE,
  numCores = 35
)
} # }
```
