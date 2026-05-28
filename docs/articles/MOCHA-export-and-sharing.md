# Exporting and sharing MOCHA results

## Introduction

After the [MOCHA workflow
tutorial](https://aifimmunology.github.io/MOCHA/articles/MOCHA-workflow-tutorial.md)
and optional [downstream
workflows](https://aifimmunology.github.io/MOCHA/articles/MOCHA-downstream-workflows.md),
you may want to share results with collaborators or visualize tracks in
IGV or the MOCHA Dash app.

This vignette ports the narrative from
`inst/notebooks/MOCHA-Export-Vignette.ipynb` into BiocStyle R Markdown.
All chunks are reference-only because exports involve file I/O, large
fragment data, and optional heavy dependencies.

``` r

library(MOCHA)
```

**Prerequisites:** `tileResults` from
[`callOpenTiles()`](https://aifimmunology.github.io/MOCHA/reference/callOpenTiles-methods.md)
and `sampleTileMatrices` from
[`getSampleTileMatrix()`](https://aifimmunology.github.io/MOCHA/reference/getSampleTileMatrix.md).

## Pack and unpack MOCHA objects

[`packMOCHA()`](https://aifimmunology.github.io/MOCHA/reference/packMOCHA.md)
archives a MOCHA object and its on-disk outputs into a single zip for
transfer. New MOCHA runs store per-sample accessibility and insertion
tracks as bigWig files under `tracks/{cellPop}/{Accessibility|Insertions}/`
inside the object directory; legacy `{cellPop}_CoverageFiles.RDS` bundles
are still read automatically when present.
[`unpackMOCHA()`](https://aifimmunology.github.io/MOCHA/reference/unpackMOCHA.md)
restores the object and rewrites stored paths.

``` r

outDir <- tempdir()
zipPath <- packMOCHA(
  MOCHAObj = tileResults,
  zipfile = file.path(outDir, "mocha_results.zip")
)
unpacked <- unpackMOCHA(
  zipfile = zipPath,
  exdir = file.path(outDir, "unpacked_mocha")
)
```

If you move the project directory after unpacking, update stored paths:

``` r

tileResults <- updateDirectoryPath(
  tileResults,
  directoryPath = file.path(outDir, "unpacked_mocha")
)
```

## Migrate legacy coverage bundles

If you have an older MOCHA output directory with
`{cellPop}_CoverageFiles.RDS` bundles, convert to the structured bigWig
layout with [`migrateCoverageToTracks()`](https://aifimmunology.github.io/MOCHA/reference/migrateCoverageToTracks.md):

``` r

tileResults <- migrateCoverageToTracks(
  MOCHAObj = tileResults,
  cellPopulations = "ALL",
  removeLegacy = FALSE,
  verbose = TRUE
)
```

If structured tracks are incomplete and no legacy RDS remains, re-run
`callOpenTiles(force = TRUE)` or `migrateCoverageToTracks(force = TRUE)`.

## Structured track filenames

On-disk bigWig files under `tracks/` use percent-encoded sample names
(for example `Sample%20A.bw` for sample `"Sample A"`). When
[`exportCoverage()`](https://aifimmunology.github.io/MOCHA/reference/exportCoverage.md)
writes sample-specific bigWigs for IGV or other browsers, export filenames
use the legacy double-underscore convention (for example
`C3__Sample__A_Coverage.bw`). MOCHA copies existing structured tracks when
possible instead of recomputing coverage.

If both structured tracks and legacy `{cellPop}_CoverageFiles.RDS` bundles
exist, MOCHA reads the structured layout and warns that the legacy file
can be removed after migration.

## Export sample-specific coverage

[`exportCoverage()`](https://aifimmunology.github.io/MOCHA/reference/exportCoverage.md)
writes bigWig files from fragment data in a tile-results object. This
step can be slow and disk-intensive; use `numCores` and a dedicated
output directory.

``` r

exportCoverage(
  SampleTileObject = tileResults,
  dir = "./data/sample_specific_coverage/",
  cellPopulations = c("C2", "C5"),
  sampleSpecific = TRUE,
  saveFile = TRUE,
  numCores = 4
)
```

## Export open tiles

``` r

exportOpenTiles(
  SampleTileObject = sampleTileMatrices,
  outDir = "./data/tiles_samplespecific/"
)
```

## Export differential tiles

Load or compute differential peaks
([`getDifferentialAccessibleTiles()`](https://aifimmunology.github.io/MOCHA/reference/getDifferentialAccessibleTiles.md)),
then export accessibility tracks for foreground/background contrasts.

``` r

# diffsGR <- readRDS("./data/my_differential_peaks.rds")
exportDifferentials(
  SampleTileObject = sampleTileMatrices,
  DifferentialsGRList = list(C2 = diffs),
  outDir = "./data/tiles_differential/"
)
```

## Export motifs

Annotate peaks with PWMs
([`addMotifSet()`](https://aifimmunology.github.io/MOCHA/reference/addMotifSet.md)),
then export motif positions for visualization. See also [Alternative TSS
usage and TF
regulation](https://aifimmunology.github.io/MOCHA/articles/Alternative-TSS-TF-regulation.md).

``` r

library(chromVARmotifs)
motifsGRanges <- addMotifSet(
  SampleTileObj = sampleTileMatrices,
  motifPWMs = chromVARmotifs::human_pwms_v2,
  returnSTM = FALSE,
  motifSetName = "CISBP"
)

exportMotifs(
  SampleTileObject = tileResults,
  motifsGRanges = unlist(motifsGRanges),
  motifSetName = "CISBP",
  outDir = "./data/motifs/",
  filterCellTypePeaks = TRUE,
  verbose = TRUE
)
```

## Export local footprints

``` r

exportLocalFootprints(
  SampleTileObj = sampleTileMatrices,
  cellPopulations = "C2",
  outDir = "./data/footprints/",
  numCores = 4
)
```

## Interactive tracks (Dash app)

For an interactive genome browser built on exported tracks, see
`inst/notebooks/Tracks-Dash-App-Generic.ipynb` in the MOCHA source
repository.

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
