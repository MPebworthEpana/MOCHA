# Execute a pilot run of model on a subset of data

\`r lifecycle::badge("deprecated")\` This function is deprecated -
improved modeling functions can be found in the package "ChAI" at
https://github.com/aifimmunology/ChAI `pilotLMEM` Runs linear
mixed-effects modeling for zero inflated data using
[`glmmTMB`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html). TryCatch will
catch errors, and return the error and dataframe for troubleshooting.

## Usage

``` r
pilotZIGLMM(
  TSAM_Object,
  cellPopulation = NULL,
  continuousFormula = NULL,
  ziformula = NULL,
  zi_threshold = 0,
  verbose = FALSE,
  pilotIndices = 1:10
)
```

## Arguments

- TSAM_Object:

  A SummarizedExperiment object generated from getSampleTileMatrix,
  chromVAR, or other.

- cellPopulation:

  A single cell population on which to run this pilot model

- continuousFormula:

  The formula, see
  [`glmmTMB`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html). Combined
  fixed and random effects formula, following lme4 syntax.

- ziformula:

  The zero-inflated formula, see
  [`glmmTMB`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html). a one-sided
  (i.e., no response variable) formula for zero-inflation combining
  fixed and random effects: the default ~0 specifies no zero-inflation.
  Specifying ~. sets the zero-inflation formula identical to the
  right-hand side of formula (i.e., the conditional effects formula);
  terms can also be added or subtracted. When using ~. as the
  zero-inflation formula in models where the conditional effects formula
  contains an offset term, the offset term will automatically be
  dropped. The zero-inflation model uses a logit link.

- zi_threshold:

  Zero-inflated threshold (range = 0-1), representing the fraction of
  samples with zeros. When the percentage of zeros in the tile is
  between 0 and zi_threshold, samples with zeroes are dropped and only
  the continous formula is used. Use this parameter at your own risk.
  Default is 0.

- verbose:

  Set TRUE to display additional messages. Default is FALSE.

- pilotIndices:

  A vector of integers defining the subset of the ExperimentObj matrix.
  Default is 1:10.

## Value

modelList a list of outputs from glmmTMB::glmmTMB
