# Test the enrichment of a given TF motif set against a motif set downstream of multiple ligands

This analogous to Gene Set Enrichment Analysis. Instead of testing for
enrichment of a geneset with a given gene set in a pathway, we are
testing the enrichment of a given TF motif set against a motif set
downstream of multiple ligands. If there is enrichment, it's a sign that
that ligand could drive that set of motifs.

## Usage

``` r
MotifSetEnrichmentAnalysis(
  ligandTFMatrix,
  motifEnrichmentDF,
  motifColumn,
  ligands,
  statColumn,
  statThreshold,
  annotationName = "CellType",
  annotation = "none",
  numCores = 1,
  verbose = FALSE
)
```

## Arguments

- ligandTFMatrix:

  NicheNet Ligand-TF matrix

- motifEnrichmentDF:

  Dataframe (unfiltered) from ArchR's peakAnnoEnrich step. Expected to
  have a column with motif names, and a column with the -log10 adjusted
  p-values.

- motifColumn:

  Column name within the motifEnrichmentDF that has motif names.

- ligands:

  Vector of ligands to test

- statColumn:

  Column name in motifEnrichmentDF containing the statistic to test

- statThreshold:

  Significance threshold used to select significant motif set

- annotationName:

  Optional column name for the annotation. Default is "CellType".

- annotation:

  Optional annotation value added to all rows of the output motif
  dataframe. Can be character vector or numeric. Default is "none".

- numCores:

  The number of cores to use with multiprocessing. Default is 1.

- verbose:

  Set TRUE to display additional messages. Default is FALSE.

## Value

specDF A dataframe containing enrichment analysis results
