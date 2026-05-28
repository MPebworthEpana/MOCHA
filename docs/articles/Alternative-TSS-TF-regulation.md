# Alternative TSS usage and TF regulation

## Introduction

After differential testing, complete the [MOCHA workflow
tutorial](https://aifimmunology.github.io/MOCHA/articles/MOCHA-workflow-tutorial.md)
through
[`getDifferentialAccessibleTiles()`](https://aifimmunology.github.io/MOCHA/reference/getDifferentialAccessibleTiles.md)
so you have `sampleTileMatrices` and `differentials` (`GRanges` with
`FDR` and `Log2FC_C`). MOCHA then offers three downstream-analysis
modules that correspond to Figures 4–5 of the manuscript:

1.  **Motif enrichment** within differentially accessible regions.
2.  **Alternative TSS regulation** — identifying genes whose multiple
    TSSs respond differently across conditions.
3.  **Motif footprinting** — average insertion profiles around TF motif
    instances within open chromatin.

All chunks in this vignette are shown with `eval = FALSE` because the
bundled example data is intentionally minimal. They are meant as a
template you can adapt to your own MOCHA objects.

## Setup

``` r

library(MOCHA)
```

We assume you already have:

- `sampleTileMatrices` — a `RangedSummarizedExperiment` from
  [`getSampleTileMatrix()`](https://aifimmunology.github.io/MOCHA/reference/getSampleTileMatrix.md).
- `differentials` — a `GRanges` from
  [`getDifferentialAccessibleTiles()`](https://aifimmunology.github.io/MOCHA/reference/getDifferentialAccessibleTiles.md)
  with `FDR` and `Log2FC_C` columns.
- TxDb / OrgDb references registered in your MOCHA object metadata.

## Motif annotation with `addMotifSet()`

Before enrichment or footprinting, attach PWM positions to your
sample-tile matrix (requires `chromVARmotifs`, `motifmatchr`, and a
BSgenome):

``` r

sampleTileMatrices <- addMotifSet(
  SampleTileObj = sampleTileMatrices,
  motifPWMs = chromVARmotifs::human_pwms_v2,
  motifSetName = "CISBP"
)
# Motif positions are stored in metadata(sampleTileMatrices)$CISBP
```

## 1. Motif enrichment in differential regions

[`MotifEnrichment()`](https://aifimmunology.github.io/MOCHA/reference/MotifEnrichment.md)
runs a hypergeometric enrichment test for each motif in a `GRangesList`
of motif positions against a background of non-target regions. A typical
workflow is:

1.  Split `differentials` into a target set (e.g. FDR ≤ 0.1, \|Log2FC\|
    ≥ 1) and a background set (the remaining tested tiles).
2.  Provide a `motifPosList` — a named `GRangesList` keyed by motif ID
    (e.g. JASPAR identifiers). MOCHA’s
    [`addMotifSet()`](https://aifimmunology.github.io/MOCHA/reference/addMotifSet.md)
    annotates a SampleTileMatrix with motif positions; the positions can
    then be pulled out for enrichment testing.

``` r

sig <- differentials[!is.na(differentials$FDR) &
                     differentials$FDR <= 0.1 &
                     abs(differentials$Log2FC_C) >= 1]
bg <- differentials[!(differentials %in% sig)]

# motifPosList is a named GRangesList of motif binding sites.
# In practice you would build this once and reuse it.
enr <- MotifEnrichment(
  Group1       = sig,
  Group2       = bg,
  motifPosList = motifPosList
)
head(enr[order(enr$adjp_val), ])
```

[`MotifSetEnrichmentAnalysis()`](https://aifimmunology.github.io/MOCHA/reference/MotifSetEnrichmentAnalysis.md)
is a higher-level wrapper that takes a ligand–TF matrix and an existing
motif-enrichment table, and identifies upstream regulators whose motif
set is over-represented in the target peaks. Use it when you want to go
from “these peaks moved” to “these TFs are likely driving the change”:

``` r

upstream <- MotifSetEnrichmentAnalysis(
  ligandTFMatrix     = ligandTFMatrix,
  motifEnrichmentDF  = enr,
  # additional parameters control significance thresholds and TF filtering
)
head(upstream)
```

## 2. Alternative TSS regulation (`getAltTSS`)

Genes can have multiple TSSs that respond independently to a treatment.
[`getAltTSS()`](https://aifimmunology.github.io/MOCHA/reference/getAltTSS.md)
annotates each differential peak with its overlapping TSS (if any),
groups peaks by gene, and flags genes where the TSS-level behaviour is
mixed:

- **Type i** — only some of the open TSSs are significantly more or less
  accessible.
- **Type ii** — multiple TSSs are significant but with *opposite* signs.

``` r

altTSS <- getAltTSS(
  completeDAPs = differentials,   # GRanges of all tested tiles with FDR + Log2FC_C
  threshold    = 0.2,             # FDR cutoff for "significant"
  TxDb         = "TxDb.Hsapiens.UCSC.hg38.knownGene",
  OrgDb        = "org.Hs.eg.db"
)
table(altTSS$type)               # type i / type ii counts
head(altTSS[altTSS$type == "ii", ])  # opposite-direction TSS pairs
```

Set `returnAllTSS = TRUE` to get the annotated long table of every TSS
overlap (one row per TSS-tile pair) without filtering. That’s the form
you want when feeding the result into
[`plotRegion()`](https://aifimmunology.github.io/MOCHA/reference/plotRegion.md)
to visualise a specific locus across conditions:

``` r

# Pick a candidate alternative TSS gene from the table above:
candidate <- "MYD88"
regionGR <- GenomicRanges::GRanges("chr3:38179000-38186000")

countSE <- getCoverage(
  sampleTileMatrices,
  cellPopulations = "C2",
  regions         = regionGR,
  groupColumn     = "COVID_status"
)

plotRegion(countSE = countSE, whichGene = candidate)
```

## 3. Motif footprinting

[`motifFootprint()`](https://aifimmunology.github.io/MOCHA/reference/motifFootprint.md)
returns a per-position average of normalised insertions in a window
around motif centres, optionally split by group. You typically run it
after
[`addMotifSet()`](https://aifimmunology.github.io/MOCHA/reference/addMotifSet.md)
has annotated the `SampleTileMatrix` with motif positions:

``` r

fp <- motifFootprint(
  SampleTileObj   = sampleTileMatrices,
  motifName       = "Motifs",          # the slot populated by addMotifSet()
  specMotif       = "MA0080.4_SPI1",   # JASPAR ID of the motif of interest
  cellPopulations = "CD16 Mono",
  windowSize      = 500,
  normTn5         = TRUE,
  smoothTn5       = 10,
  groupColumn     = "COVID_status"
)
fp
```

The returned `SummarizedExperiment` carries a per-position-by-sample
matrix of normalised insertion counts. Aggregate it with
`SummarizedExperiment::assays(fp)` and plot the mean profile across
groups; or pass `sampleSpecific = TRUE` to get one curve per sample.

[`plotMotifs()`](https://aifimmunology.github.io/MOCHA/reference/plotMotifs.md)
visualizes footprint profiles and can return summary statistics when
`returnDF = TRUE`:

``` r

fp_stats <- plotMotifs(
  fp,
  footprint = "MA0080.4_SPI1",
  groupColumn = "COVID_status",
  returnDF = TRUE,
  plotIndividualRegions = FALSE
)
head(fp_stats)
```

## Putting it together

A common Figure 4–5 narrative threads these three modules:

1.  Find differential tiles →
    [`getDifferentialAccessibleTiles()`](https://aifimmunology.github.io/MOCHA/reference/getDifferentialAccessibleTiles.md).
2.  Score motif enrichment in the differentials →
    [`MotifEnrichment()`](https://aifimmunology.github.io/MOCHA/reference/MotifEnrichment.md) +
    [`MotifSetEnrichmentAnalysis()`](https://aifimmunology.github.io/MOCHA/reference/MotifSetEnrichmentAnalysis.md)
    to identify candidate upstream TFs.
3.  Drill into specific candidate loci →
    [`getAltTSS()`](https://aifimmunology.github.io/MOCHA/reference/getAltTSS.md)
    for alternative TSS regulation,
    [`plotRegion()`](https://aifimmunology.github.io/MOCHA/reference/plotRegion.md)
    for locus-level visualisation.
4.  Validate TF activity at the binding-site level →
    [`motifFootprint()`](https://aifimmunology.github.io/MOCHA/reference/motifFootprint.md)
    to show insertion-profile shifts at the motif’s binding sites.

## Session information

``` r

sessionInfo()
#> R version 4.5.3 (2026-03-11)
#> Platform: x86_64-conda-linux-gnu
#> Running under: Ubuntu 24.04.2 LTS
#> 
#> Matrix products: default
#> BLAS/LAPACK: /home/enki/miniforge3/envs/mocha-test/lib/libopenblasp-r0.3.33.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: America/Los_Angeles
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] MOCHA_2.0.0      BiocStyle_2.38.0
#> 
#> loaded via a namespace (and not attached):
#>   [1] RColorBrewer_1.1-3          jsonlite_2.0.0             
#>   [3] MultiAssayExperiment_1.36.1 magrittr_2.0.5             
#>   [5] spatstat.utils_3.2-3        farver_2.1.2               
#>   [7] rmarkdown_2.31              fs_2.1.0                   
#>   [9] ragg_1.5.1                  vctrs_0.7.3                
#>  [11] ROCR_1.0-12                 spatstat.explore_3.8-0     
#>  [13] htmltools_0.5.9             S4Arrays_1.10.1            
#>  [15] BiocBaseUtils_1.12.0        SparseArray_1.10.8         
#>  [17] sass_0.4.10                 sctransform_0.4.3          
#>  [19] parallelly_1.47.0           KernSmooth_2.23-26         
#>  [21] bslib_0.11.0                htmlwidgets_1.6.4          
#>  [23] desc_1.4.3                  ica_1.0-3                  
#>  [25] plyr_1.8.9                  plotly_4.12.0              
#>  [27] zoo_1.8-15                  cachem_1.1.0               
#>  [29] igraph_2.3.1                mime_0.13                  
#>  [31] lifecycle_1.0.5             pkgconfig_2.0.3            
#>  [33] Matrix_1.7-5                R6_2.6.1                   
#>  [35] fastmap_1.2.0               MatrixGenerics_1.22.0      
#>  [37] fitdistrplus_1.2-6          future_1.70.0              
#>  [39] shiny_1.13.0                digest_0.6.39              
#>  [41] patchwork_1.3.2             S4Vectors_0.48.0           
#>  [43] tensor_1.5.1                Seurat_5.5.0               
#>  [45] RSpectra_0.16-2             irlba_2.3.7                
#>  [47] textshaping_1.0.5           GenomicRanges_1.62.1       
#>  [49] progressr_0.19.0            spatstat.sparse_3.1-0      
#>  [51] polyclip_1.10-7             httr_1.4.8                 
#>  [53] abind_1.4-8                 compiler_4.5.3             
#>  [55] S7_0.2.2                    fastDummies_1.7.6          
#>  [57] MASS_7.3-65                 DelayedArray_0.36.0        
#>  [59] tools_4.5.3                 lmtest_0.9-40              
#>  [61] otel_0.2.0                  httpuv_1.6.17              
#>  [63] future.apply_1.20.2         goftest_1.2-3              
#>  [65] glue_1.8.1                  nlme_3.1-169               
#>  [67] promises_1.5.0              grid_4.5.3                 
#>  [69] Rtsne_0.17                  cluster_2.1.8.2            
#>  [71] reshape2_1.4.5              generics_0.1.4             
#>  [73] gtable_0.3.6                spatstat.data_3.1-9        
#>  [75] tidyr_1.3.2                 data.table_1.17.8          
#>  [77] sp_2.2-1                    XVector_0.50.0             
#>  [79] spatstat.geom_3.7-3         BiocGenerics_0.56.0        
#>  [81] RcppAnnoy_0.0.23            ggrepel_0.9.8              
#>  [83] RANN_2.6.2                  pillar_1.11.1              
#>  [85] stringr_1.6.0               spam_2.11-3                
#>  [87] RcppHNSW_0.6.0              later_1.4.8                
#>  [89] splines_4.5.3               dplyr_1.2.1                
#>  [91] lattice_0.22-9              deldir_2.0-4               
#>  [93] survival_3.8-6              tidyselect_1.2.1           
#>  [95] miniUI_0.1.2                pbapply_1.7-4              
#>  [97] knitr_1.51                  gridExtra_2.3              
#>  [99] bookdown_0.46               IRanges_2.44.0             
#> [101] Seqinfo_1.0.0               SummarizedExperiment_1.40.0
#> [103] scattermore_1.2             stats4_4.5.3               
#> [105] xfun_0.57                   Biobase_2.70.0             
#> [107] matrixStats_1.5.0           stringi_1.8.7              
#> [109] lazyeval_0.2.3              yaml_2.3.12                
#> [111] evaluate_1.0.5              codetools_0.2-20           
#> [113] tibble_3.3.1                BiocManager_1.30.27        
#> [115] cli_3.6.6                   uwot_0.2.4                 
#> [117] xtable_1.8-8                reticulate_1.46.0          
#> [119] systemfonts_1.3.2           jquerylib_0.1.4            
#> [121] dichromat_2.0-0.1           Rcpp_1.1.1-1.1             
#> [123] spatstat.random_3.4-5       globals_0.19.1             
#> [125] png_0.1-9                   spatstat.univar_3.1-7      
#> [127] parallel_4.5.3              pkgdown_2.2.0              
#> [129] ggplot2_4.0.3               dotCall64_1.2              
#> [131] listenv_0.10.1              viridisLite_0.4.3          
#> [133] scales_1.4.0                ggridges_0.5.7             
#> [135] SeuratObject_5.4.0          purrr_1.2.2                
#> [137] rlang_1.2.0                 cowplot_1.2.0
```
