# Get predicted P(zero) dropout scores for a cell population

Get predicted P(zero) dropout scores for a cell population

## Usage

``` r
getDropoutProb(TSAM_Object, cellPopulation)
```

## Arguments

- TSAM_Object:

  A SummarizedExperiment with `DropoutProb_*` assays from
  [`assessDropout`](https://aifimmunology.github.io/MOCHA/reference/assessDropout.md).

- cellPopulation:

  Cell population (intensity assay) name.

## Value

A matrix of predicted `P(zero)` for each tile and sample (lower values
on observed zeros indicate more likely technical dropout).
