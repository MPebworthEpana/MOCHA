# Convert Seurat/Signac object to sample-level MOCHA fragment inputs.

Converts a Seurat object with a Signac `ChromatinAssay` into a
sample-level `GRangesList` and matching `cellColData` for use with
[`callOpenTiles`](https://aifimmunology.github.io/MOCHA/reference/callOpenTiles-methods.md).

## Usage

``` r
seuratToMOCHAInputs(
  seuratObj,
  assay = NULL,
  sampleColumn = "Sample",
  cellPopLabel,
  cellCol = "RG",
  fragmentPathColumn = NULL,
  dropEmptySamples = TRUE,
  verbose = FALSE
)
```

## Arguments

- seuratObj:

  A Seurat object with a ChromatinAssay.

- assay:

  Name of the ChromatinAssay. Default uses
  `Seurat::DefaultAssay(seuratObj)`.

- sampleColumn:

  Metadata column with biological sample IDs.

- cellPopLabel:

  Metadata column with cell population labels.

- cellCol:

  Column name for cell barcodes in fragment `mcols`. Default is `"RG"`
  (MOCHA/ArchR convention).

- fragmentPathColumn:

  Optional metadata column with per-cell fragment file paths. If
  provided, takes precedence over Signac `Fragments()`. File barcodes
  must match `rownames(meta)` unless `Fragment@cells` is available on
  the assay (then file barcodes are translated automatically).

- dropEmptySamples:

  If `TRUE`, drop samples with no fragments and warn.

- verbose:

  Print progress messages.

## Value

A list with `ATACFragments` (`GRangesList`) and `cellColData`
(`data.frame`).
