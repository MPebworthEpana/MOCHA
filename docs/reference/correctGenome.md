# `correctGenome`

`correctGenome` Runs ChromVAR on a MOCHA Tile-Sample Accessibility
Matrix object

## Usage

``` r
correctGenome(TSAM, genome)
```

## Arguments

- TSAM:

  RangedSummarizedExperiment object, containing the TSAM from
  getSampleTileMatrix, with a motif set added via addMotifSet()

- genome:

  genome to correct with.

## Value

chromVAR object

## Details

This is a wrapper for basic SummarizedExperiment-based subsetting.
