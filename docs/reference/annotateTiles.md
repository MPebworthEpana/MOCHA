# Annotate tiles with gene annotations

`annotateTiles` annotates a set of sample-tile matrices given with gene
annotations. Details on TxDb and Org annotation packages and available
annotations can be found at Bioconductor:
https://bioconductor.org/packages/3.15/data/annotation/

## Usage

``` r
annotateTiles(Obj, TxDb = NULL, Org = NULL, promoterRegion = c(2000, 100))
```

## Arguments

- Obj:

  A RangedSummarizedExperiment generated from getSampleTileMatrix,
  containing TxDb and Org in the metadata. This may also be a GRanges
  object.

- TxDb:

  The annotation package for TxDb object for your genome. Optional, only
  required if Obj is a GRanges.

- Org:

  The genome-wide annotation for your organism. Optional, only required
  if Obj is a GRanges.

- promoterRegion:

  Optional list containing the window size in basepairs defining the
  promoter region. The format is (upstream, downstream). Default is
  (2000, 100).

## Value

Obj, the input data structure with added gene annotations (whether
GRanges or SampleTileObj)

## Examples

``` r
# \donttest{
if (
  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE) &&
    requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)
) {
  tiles <- MOCHA::callOpenTiles(
    ATACFragments = MOCHA::exampleFragments,
    cellColData = MOCHA::exampleCellColData,
    blackList = MOCHA::exampleBlackList,
    genome = "BSgenome.Hsapiens.UCSC.hg19",
    TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    OrgDb = "org.Hs.eg.db",
    outDir = tempdir(),
    cellPopLabel = "Clusters",
    cellPopulations = "C2",
    numCores = 1
  )
  stm <- MOCHA::getSampleTileMatrix(tiles, cellPopulations = "C2", threshold = 0)
  annotated <- MOCHA::annotateTiles(stm)
}
#> 
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> 
#> Attaching package: ‘Biostrings’
#> The following object is masked from ‘package:base’:
#> 
#>     strsplit
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> Warning: GRanges object contains 1342 out-of-bound ranges located on sequences chr1
#>   and chr2. Note that ranges located on a sequence whose length is unknown (NA)
#>   or on a circular sequence are not considered out-of-bound (use seqlengths()
#>   and isCircular() to get the lengths and circularity flags of the underlying
#>   sequences). You can use trim() to trim these ranges. See
#>   ?`trim,GenomicRanges-method` for more information.
#> Warning: GRanges object contains 1490 out-of-bound ranges located on sequences chr1
#>   and chr2. Note that ranges located on a sequence whose length is unknown (NA)
#>   or on a circular sequence are not considered out-of-bound (use seqlengths()
#>   and isCircular() to get the lengths and circularity flags of the underlying
#>   sequences). You can use trim() to trim these ranges. See
#>   ?`trim,GenomicRanges-method` for more information.
#> Warning: The `threshold` argument of `getSampleTileMatrix()` is deprecated as of MOCHA
#> 2.0.0.
#> ℹ Please use the `reproducibilityThreshold` argument instead.
#> 'select()' returned 1:1 mapping between keys and columns
#> Warning: GRanges object contains 222 out-of-bound ranges located on sequences
#>   chr19_KI270868v1_alt, chr1_KI270762v1_alt, chr4_GL000257v2_alt,
#>   chr16_KI270728v1_random, chr22_KI270731v1_random, chr13_KI270838v1_alt,
#>   chr14_KI270847v1_alt, chr17_JH159146v1_alt, chr3_KI270781v1_alt,
#>   chr4_ML143349v1_fix, chr3_GL000221v1_random, chr1_KI270706v1_random,
#>   chr16_GL383556v1_alt, chr7_KI270809v1_alt, chr19_KI270922v1_alt,
#>   chr19_KI270929v1_alt, chr2_GL383522v1_alt, chr4_KI270788v1_alt,
#>   chr12_GL383553v2_alt, chr12_KI270834v1_alt, chrUn_KI270748v1,
#>   chr22_KI270879v1_alt, chr3_KI270777v1_alt, chr5_KV575244v1_fix,
#>   chr16_KI270854v1_alt, chr21_GL383581v2_alt, chr15_KI270851v1_alt,
#>   chr17_KI270910v1_alt, chr9_GL383540v1_alt, chr7_KI270806v1_alt,
#>   chr15_ML143370v1_fix, chr11_KI270831v1_alt, chr17_KI270857v1_alt,
#>   chr6_KI270798v1_alt, chr21_KI270872v1_alt, chr11_KI270902v1_alt,
#>   chr10_MU273367v1_fix, chr2_KI270774v1_alt, chr4_KV766193v1_alt,
#>   chr16_ML143373v1_fix, chr17_KV766196v1_fix, chr7_GL383534v2_alt,
#>   chr17_MU273383v1_fix, chr5_KI270795v1_alt, chr5_KI270898v1_alt,
#>   chr6_KI270801v1_alt, chr15_KI270850v1_alt, chr20_KI270869v1_alt,
#>   chr17_KI270860v1_alt, chr19_GL383575v2_alt, chr7_KZ208912v1_fix,
#>   chr19_KI270923v1_alt, chr17_JH159147v1_alt, chr16_GL383557v1_alt,
#>   chr17_KV575245v1_fix, chr7_KI270803v1_alt, chr8_KI270815v1_alt,
#>   chr5_GL339449v2_alt, and chr19_KI270866v1_alt. Note that ranges located on a
#>   sequence whose length is unknown (NA) or on a circular sequence are not
#>   considered out-of-bound (use seqlengths() and isCircular() to get the lengths
#>   and circularity flags of the underlying sequences). You can use trim() to
#>   trim these ranges. See ?`trim,GenomicRanges-method` for more information.
# }
```
