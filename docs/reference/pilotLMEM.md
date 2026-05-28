# Execute a pilot run of single linear model on a subset of data

\`r lifecycle::badge("deprecated")\` This function is deprecated -
improved modeling functions can be found in the package "ChAI" at
https://github.com/aifimmunology/ChAI `pilotLMEM` Runs linear
mixed-effects modeling for continuous, non-zero inflated data using
[`lmer`](https://rdrr.io/pkg/lmerTest/man/lmer.html)

## Usage

``` r
pilotLMEM(
  ExperimentObj,
  assayName,
  modelFormula,
  pilotIndices = 1:10,
  verbose = FALSE
)
```

## Arguments

- ExperimentObj:

  A SummarizedExperiment-type object generated from chromVAR,
  makePseudobulkRNA, or other. Objects from getSampleTileMatrix can
  work, but we recommend runZIGLMM for those objects, not runLMEM\>

- assayName:

  a character string, matching the name of an assay within the
  SummarizedExperiment. The assay named will be used for modeling.

- modelFormula:

  The formula to use with lmerTest::lmer, in the format (exp ~ factors).
  All factors must be found in column names of the ExperimentObj
  metadata.

- pilotIndices:

  A vector of integers defining the subset of the ExperimentObj matrix.
  Default is 1:10.

- verbose:

  Set TRUE to display additional messages. Default is FALSE.

## Value

modelList a list of outputs from lmerTest::lmer
