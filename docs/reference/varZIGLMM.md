# Zero-inflated Variance Decomposition for pseudobulked scATAC data

\`r lifecycle::badge("deprecated")\` This function is deprecated -
improved modeling functions can be found in the package "ChAI" at
https://github.com/aifimmunology/ChAI `varZIGLMM` Identified variance
decomposition on a given cell type across both zero-inflated and
continuous space using a zero-inflated general linear mixed model
[`glmmTMB`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html)

## Usage

``` r
varZIGLMM(
  TSAM_Object,
  cellPopulation = NULL,
  continuousRandom = NULL,
  ziRandom = NULL,
  zi_threshold = 0.1,
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

- continuousRandom:

  Random effects to test in the continuous portion. All factors must be
  found in column names of the TSAM_Object metadata, except for
  FragNumber and CellCount, which will be extracted from the
  TSAM_Object's metadata.

- ziRandom:

  Random effects to test in the zero-inflated portion. All factors must
  be found in column names of the TSAM_Object colData metadata, except
  for FragNumber and CellCount, which will be extracted from the
  TSAM_Object's metadata.

- zi_threshold:

  Zero-inflated threshold ( range = 0-1), representing the fraction of
  samples with zeros. When the percentage of zeros in the tile is
  between 0 and zi_threshold, samples with zeroes are dropped and only
  the continous formula is used. Use this parameter at your own risk.
  Default is 0.

- verbose:

  Set TRUE to display additional messages. Default is FALSE.

- numCores:

  integer. Number of cores to parallelize across.

## Value

results a SummarizedExperiment containing results from ZIGLMM (Fixed
effect estiamtes, P-values, and Std Error)

## Examples

``` r
if (FALSE) { # \dontrun{
modelList <- runZIGLMM(STM[c(1:1000), ],
  cellPopulation = "CD16 Mono",
  continuousRandom = c("Age", "Sex", "Days"),
  ziRandom = c("FragNumber", "Days"),
  verbose = TRUE,
  numCores = 35
)
} # }
```
