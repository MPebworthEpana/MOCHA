# Run Linear Mixed-Effects Modeling for continuous, non-zero inflated data

`runLMEM` Runs linear mixed-effects modeling for continuous, non-zero
inflated data using [`lmer`](https://rdrr.io/pkg/lmerTest/man/lmer.html)

## Usage

``` r
runLMEM(
  ExperimentObj,
  assayName,
  modelFormula,
  initialSampling = 5,
  verbose = FALSE,
  numCores = 2
)
```

## Arguments

- ExperimentObj:

  A SummarizedExperiment object generated from getSampleTileMatrix,
  chromVAR, or other. It is expected to contain only one assay, or only
  the first assay will be used for the model. Data should not be
  zero-inflated.

- assayName:

  The name of the assay to model within the SummarizedExperiment.

- modelFormula:

  The formula to use with lmerTest::lmer, in the format (exp ~ factors).
  All factors must be found in column names of the ExperimentObj
  metadata. modelFormula must start with 'exp' as the response. See
  [lmer](https://rdrr.io/pkg/lmerTest/man/lmer.html).

- initialSampling:

  Size of data to use for pilot

- verbose:

  Set TRUE to display additional messages. Default is FALSE.

- numCores:

  integer. Number of cores to parallelize across.

## Value

results a SummarizedExperiment containing LMEM results. Assays are
metrics related to the model coefficients, including the Estimate,
Std_Error, df, t_value, p_value. Within each assay, each row corresponds
to each row of the SummarizedExperiment and columns correspond to each
fixed effect variable within the model. Any row metadata from the
ExperimentObject (see rowData(ExperimentObj)) is preserved in the
output. The Residual matrix and the variance of the random effects are
saved in the metadata slot of the output.

## Examples
