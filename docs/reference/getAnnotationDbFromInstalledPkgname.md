# Loads and attaches an installed TxDb or OrgDb-class Annotation database package.

See
[getBSgenome](https://rdrr.io/pkg/BSgenome/man/available.genomes.html)

## Usage

``` r
getAnnotationDbFromInstalledPkgname(dbName, type)
```

## Arguments

- dbName:

  Exact name of installed annotation data package.

- type:

  Expected class of the annotation data package, must be either "OrgDb"
  or "TxDb".

## Value

the loaded Annotation database object.#' @noRd
