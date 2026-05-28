# Identify motifs within a peakset

`addMotifSet` Identify motifs within your peakset.

## Usage

``` r
addMotifSet(
  SampleTileObj,
  motifPWMs,
  w = 7,
  returnSTM = TRUE,
  motifSetName = "Motifs"
)
```

## Arguments

- SampleTileObj:

  A SummarizedExperiment, specifically the output of getSampleTileMatrix

- motifPWMs:

  A pwms object for the motif database. Either PFMatrix, PFMatrixList,
  PWMatrix, or PWMatrixList

- w:

  Parameter for motifmatchr controlling size in basepairs of window for
  filtration. Default is 7.

- returnSTM:

  If TRUE, return the modified SampleTileObj with motif set added to
  metadata (default). If FALSE, return the motifs from motifmatchr as a
  GRanges.

- motifSetName:

  Name to give motifList in the SampleTileObj's metadata if
  \`returnSTM=TRUE\`. Default is 'Motifs'.

## Value

the modified SampleTileObj with motifs added to the metadata

## Examples

``` r
if (FALSE) { # \dontrun{
# load a curated motif set from library(chromVARmotifs)
# included with ArchR installation
data(human_pwms_v2)
SE_with_motifs <- addMotifSet(
  SampleTileObj,
  motifPWMs = human_pwms_v2,
  returnSTM = TRUE, motifSetName = "Motifs", w = 7
)
} # }
```
