# Get motifset for a given cell type

Extracts an existing motif set and filters it down to only motifs found
within regions open in that celltype A necessary step for accurate motif
enrichment analysis.

## Usage

``` r
getCellTypeMotifs(
  STM,
  cellPopulation,
  MotifSetName = "Motifs",
  specMotif = NULL,
  asGRangesList = TRUE
)
```

## Arguments

- STM:

  A Sample-Tile Object from MOCHA

- cellPopulation:

  The name of a cell type within MOCHA

- MotifSetName:

  The name of the total motif set, saved as metadata within the MOCHA
  object

- specMotif:

  Name of a specific motif or sets of motif that you wish to pull out.
  Default is NULL, with returns all motifs.

- asGRangesList:

  Default is FALSE. If TRUE, then converts list of GRanges to
  GRangesList, which takes quite a while.

## Value

A GRangesList containing the motifs that are present in open regions for
a given cell type
