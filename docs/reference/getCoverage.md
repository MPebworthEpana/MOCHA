# Get sample-specific coverage files for each sample-cell population

getCoverage takes the output of MOCHA::getPopFrags and returns a GRanges
of singe-basepair resolution coverage.

## Usage

``` r
getCoverage(
  popFrags,
  normFactor,
  TxDb,
  cl,
  filterEmpty = FALSE,
  verbose = FALSE
)
```

## Arguments

- popFrags:

  GRangesList of fragments for all sample/cell populations

- normFactor:

  Normalization factor. Can be either be one, in which case all coverage
  files will be normalized by the same value, or the same length as the
  GRangesList

- TxDb:

  The TxDb-class transcript annotation package for your organism (e.g.
  "TxDb.Hsapiens.UCSC.hg38.knownGene"). This must be installed. See
  [Bioconductor AnnotationData
  Packages](https://bioconductor.org/packages/release/data/annotation/).

- cl:

  cl argument to
  [`pblapply`](https://peter.solymos.org/pbapply/reference/pbapply.html)

- filterEmpty:

  True/False flag on whether or not to carry forward regions without
  coverage.

- verbose:

  Boolean variable to determine verbosity of output.

## Value

popCounts A GRangesList of coverage for each sample and cell population
