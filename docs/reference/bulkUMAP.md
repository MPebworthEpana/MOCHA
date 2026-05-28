# Generate UMAP from pseudobulk LSI results

`bulkUMAP` generates UMAP from pseudobulk LSIObj object from
[bulkDimReduction](https://aifimmunology.github.io/MOCHA/reference/bulkDimReduction.md)

## Usage

``` r
bulkUMAP(
  SEObj,
  assay = "LSI",
  components = c(1:30),
  nNeighbors = 15,
  returnModel = FALSE,
  verbose = TRUE,
  seed = 1,
  ...
)
```

## Arguments

- SEObj:

  The SummarizedExperiment object output from bulkDimReduction, or an
  STM, subset to just one cell type.

- assay:

  A string, describing the name of the assay within SEObj to run UMAP
  ('PCA', 'LSI', or 'counts').

- components:

  A vector of integers. Number of components to include in LSI (1:30
  typically).

- nNeighbors:

  See [umap](https://jlmelville.github.io/uwot/reference/umap.html). The
  size of local neighborhood (in terms of number of neighboring sample
  points) used for manifold approximation. Default is 15.

- returnModel:

  A boolean. Default is FALSE. If set to true, it will return a list,
  where the first is the UMAP coordinates with metadata for plotting,
  and the second is the full UMAP model so further projection can occur.

- verbose:

  Set TRUE to display additional messages. Default is FALSE.

- seed:

  an integer. Represents the random seed to pass to the UMAP. Default
  seed is 1.

- ...:

  Additional arguments to be passed to
  [umap](https://jlmelville.github.io/uwot/reference/umap.html).

## Value

fullUMAP data.frame of UMAP values with metadata attached.

## Examples

``` r
if (FALSE) { # \dontrun{
UMAPvalues <- MOCHA::bulkUMAP(LSIObj)
} # }
```
