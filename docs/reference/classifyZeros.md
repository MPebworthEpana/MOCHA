# Classify observed zeros as technical, biological, or ambiguous

Classify observed zeros as technical, biological, or ambiguous

## Usage

``` r
classifyZeros(
  TSAM_Object,
  cellPopulation,
  techThreshold = 0.2,
  bioThreshold = 0.8
)
```

## Arguments

- TSAM_Object:

  A SummarizedExperiment with dropout probabilities from
  [`assessDropout`](https://aifimmunology.github.io/MOCHA/reference/assessDropout.md).

- cellPopulation:

  Cell population name.

- techThreshold:

  See
  [`assessDropout`](https://aifimmunology.github.io/MOCHA/reference/assessDropout.md).

- bioThreshold:

  See
  [`assessDropout`](https://aifimmunology.github.io/MOCHA/reference/assessDropout.md).

## Value

A matrix of the same dimensions as the intensity matrix with values
`"technical"`, `"biological"`, `"ambiguous"`, or `NA` for non-zero
observations.
