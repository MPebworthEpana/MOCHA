# Add a column to the sample-level colData of a MOCHA object

`addCellColData` adds a new column to the sample-level colData of a
MOCHA tileResults (`MultiAssayExperiment` from `callOpenTiles`) or
SampleTileMatrix (`RangedSummarizedExperiment` from
`getSampleTileMatrix`). MOCHA pseudobulks by sample x cell-population,
so colData rows are biological samples; the name mirrors
`ArchR::addCellColData` for API familiarity.

## Usage

``` r
addCellColData(object, name, value, samples = NULL, force = FALSE)
```

## Arguments

- object:

  A MOCHA tileResults or SampleTileMatrix object.

- name:

  Character scalar. Name of the column to add.

- value:

  Vector of values to add. Either length `nrow(colData)` (in which case
  it is assumed aligned with `samples`) or a named vector with names
  matching sample identifiers.

- samples:

  Optional character vector of sample identifiers that `value`
  corresponds to. Defaults to `rownames(colData(object))`. Missing
  samples will receive `NA`.

- force:

  Logical. If `TRUE`, an existing column with the same `name` is
  overwritten. Default `FALSE`.

## Value

The input object with the new colData column attached.
