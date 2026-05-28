# MOCHA workflow tutorial

## Introduction

MOCHA (Model-based single-cell Open Chromatin Analysis) provides peak
calling, sample-tile matrices, and downstream analyses for scATAC-seq
after cell-type identification. This vignette runs an end-to-end
workflow on small bundled example data (`exampleFragments`,
`exampleCellColData`, `exampleBlackList`) so it can be built on the
Bioconductor build system without an ArchR project.

For importing from ArchR, Signac, SnapATAC2, and other sources and
calling open tiles, see [Data Import Tutorial: Signac, ArchR, SnapATAC
and
more](https://aifimmunology.github.io/MOCHA/articles/Data-Import-Tutorial.md).

## Setup

``` r

library(MOCHA)
library(MultiAssayExperiment)
library(SummarizedExperiment)
```

Annotation packages used below are listed in `DESCRIPTION` under
`Suggests`.

## Peak calling with `callOpenTiles`

MOCHA accepts a per-sample `GRangesList` of fragments plus cell
metadata.

``` r

tileResults <- callOpenTiles(
  ATACFragments = exampleFragments,
  cellColData = exampleCellColData,
  blackList = exampleBlackList,
  genome = "BSgenome.Hsapiens.UCSC.hg19",
  TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
  OrgDb = "org.Hs.eg.db",
  outDir = tempdir(),
  cellPopLabel = "Clusters",
  cellPopulations = c("C2", "C5"),
  numCores = 1,
  verbose = TRUE
)
tileResults
#> A MultiAssayExperiment object of 2 listed
#>  experiments with user-defined names and respective classes.
#>  Containing an ExperimentList class object of length 2:
#>  [1] C2: RaggedExperiment with 71764 rows and 1 columns
#>  [2] C5: RaggedExperiment with 66915 rows and 1 columns
#> Functionality:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample coordination DataFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
#>  exportClass() - save data to flat files
```

Inspect called open tiles for one population:

``` r

openC2 <- getOpenTiles(tileResults, cellPopulations = "C2")
length(openC2[[1]])
#> [1] 15665
```

## Tuning parameters

Before building sample-tile matrices, inspect how peak reproducibility
varies across samples.
[`plotConsensus()`](https://aifimmunology.github.io/MOCHA/reference/plotConsensus.md)
and
[`suggestConsensusThreshold()`](https://aifimmunology.github.io/MOCHA/reference/suggestConsensusThreshold.md)
help you choose `reproducibilityThreshold` for
[`getSampleTileMatrix()`](https://aifimmunology.github.io/MOCHA/reference/getSampleTileMatrix.md).

``` r

consensusDF <- plotConsensus(
  tileResults,
  cellPopulations = c("C2", "C5"),
  returnDFs = TRUE,
  numCores = 1
)
head(consensusDF$C2)
#>   Reproducibility PeakNumber
#> 1               0      71764
#> 2               1      15665
```

``` r

suggested <- suggestConsensusThreshold(
  tileResults,
  cellPopulations = c("C2", "C5"),
  method = "kneedle",
  numCores = 1
)
suggested
#> [1] CellPopulation  Reproducibility PeakNumber      Method         
#> <0 rows> (or 0-length row.names)
```

The table below maps common tuning knobs to the functions that help you
inspect them. For dropout-aware differentials, see [MOCHA downstream
workflows](https://aifimmunology.github.io/MOCHA/articles/MOCHA-downstream-workflows.md).

| Knob | Where set | Helper to inspect |
|----|----|----|
| `reproducibilityThreshold` | [`getSampleTileMatrix()`](https://aifimmunology.github.io/MOCHA/reference/getSampleTileMatrix.md) | [`plotConsensus()`](https://aifimmunology.github.io/MOCHA/reference/plotConsensus.md), [`suggestConsensusThreshold()`](https://aifimmunology.github.io/MOCHA/reference/suggestConsensusThreshold.md) |
| `signalThreshold` | tile filtering in STM pipeline | [`plotIntensityDistribution()`](https://aifimmunology.github.io/MOCHA/reference/plotIntensityDistribution.md) |
| `method` / `pairColumn` | [`getDifferentialAccessibleTiles()`](https://aifimmunology.github.io/MOCHA/reference/getDifferentialAccessibleTiles.md) | — |
| `dropoutAdjustment` | [`getDifferentialAccessibleTiles()`](https://aifimmunology.github.io/MOCHA/reference/getDifferentialAccessibleTiles.md) | [`assessDropout()`](https://aifimmunology.github.io/MOCHA/reference/assessDropout.md) |

## Sample-tile matrices

Consensus tiles across samples are summarized with
[`getSampleTileMatrix()`](https://aifimmunology.github.io/MOCHA/reference/getSampleTileMatrix.md).
We use `reproducibilityThreshold = 0` here for a permissive union; in
practice you may use the value from
[`suggestConsensusThreshold()`](https://aifimmunology.github.io/MOCHA/reference/suggestConsensusThreshold.md)
above.

``` r

sampleTileMatrices <- getSampleTileMatrix(
  tileResults,
  cellPopulations = c("C2", "C5"),
  reproducibilityThreshold = 0
)
sampleTileMatrices
#> class: RangedSummarizedExperiment 
#> dim: 25112 1 
#> metadata(6): summarizedData Genome ... Directory History
#> assays(2): C2 C5
#> rownames(25112): chr1:1000000-1000499 chr1:10002000-10002499 ...
#>   chr2:99954000-99954499 chr2:99954500-99954999
#> rowData names(2): C2 C5
#> colnames(1): PBMCSmall
#> colData names(1): X
```

``` r

# Example using a non-zero reproducibility threshold from suggestConsensusThreshold():
# thr <- suggested$Reproducibility[suggested$CellPopulation == "C2"]
# sampleTileMatrices <- getSampleTileMatrix(
#   tileResults,
#   cellPopulations = c("C2", "C5"),
#   reproducibilityThreshold = thr,
#   numCores = 1
# )
```

Subset to one cell population:

``` r

c2Matrix <- subsetMOCHAObject(
  sampleTileMatrices,
  subsetBy = "celltype",
  groupList = "C2"
)
dim(c2Matrix)
#> [1] 15665     1
```

``` r

plotIntensityDistribution(c2Matrix, cellPopulation = "C2", returnDF = TRUE)
```

## Annotate tiles

Gene annotations are added from TxDb/OrgDb stored in object metadata.

``` r

annotated <- annotateTiles(c2Matrix)
head(SummarizedExperiment::rowData(annotated))
#> DataFrame with 6 rows and 4 columns
#>                               C2        C5    tileType        Gene
#>                        <logical> <logical> <character> <character>
#> chr1:10002500-10002999      TRUE      TRUE  Intragenic        RBP7
#> chr1:10003000-10003499      TRUE      TRUE  Intragenic        RBP7
#> chr1:10003500-10003999      TRUE      TRUE  Intragenic        RBP7
#> chr1:10010000-10010499      TRUE     FALSE  Intragenic        RBP7
#> chr1:10010500-10010999      TRUE      TRUE  Intragenic        RBP7
#> chr1:10011000-10011499      TRUE      TRUE  Intragenic        RBP7
```

### Motif annotations (optional)

To attach PWM motif positions to a sample-tile matrix for enrichment or
footprinting, use
[`addMotifSet()`](https://aifimmunology.github.io/MOCHA/reference/addMotifSet.md)
(requires `chromVARmotifs`, `motifmatchr`, and a BSgenome). See
[Alternative TSS usage and TF
regulation](https://aifimmunology.github.io/MOCHA/articles/Alternative-TSS-TF-regulation.md).

``` r

# annotatedWithMotifs <- addMotifSet(
#   SampleTileObj = annotated,
#   motifPWMs = chromVARmotifs::human_pwms_v2,
#   motifSetName = "CISBP"
# )
```

## Differential accessibility (illustration)

When sample metadata includes a grouping column, differential testing
compares foreground and background sample sets. The example data uses
the `Sample` column; we compare the first two samples as a minimal
demonstration.

``` r

samples <- unique(SummarizedExperiment::colData(c2Matrix)$Sample)
if (length(samples) >= 2L) {
  diffs <- getDifferentialAccessibleTiles(
    SampleTileObj = c2Matrix,
    cellPopulations = "C2",
    groupColumn = "Sample",
    foreground = samples[1],
    background = samples[2],
    numCores = 1,
    verbose = TRUE
  )
  head(diffs)
}
```

### Advanced differential options

[`getDifferentialAccessibleTiles()`](https://aifimmunology.github.io/MOCHA/reference/getDifferentialAccessibleTiles.md)
supports additional arguments for paired designs and dropout handling:

``` r

# diffs_paired <- getDifferentialAccessibleTiles(
#   SampleTileObj = c2Matrix,
#   cellPopulations = "C2",
#   groupColumn = "Sample",
#   foreground = samples[1],
#   background = samples[2],
#   method = "wilcoxon_paired",
#   pairColumn = "Subject",       # column linking paired samples
#   dropoutAdjustment = "biological_only",
#   bioThreshold = 0.8,
#   numCores = 4
# )
```

When `dropoutAdjustment` is set and dropout probabilities are not
already stored on the object, MOCHA runs
[`assessDropout()`](https://aifimmunology.github.io/MOCHA/reference/assessDropout.md)
automatically.

## Next steps

- [MOCHA downstream
  workflows](https://aifimmunology.github.io/MOCHA/articles/MOCHA-downstream-workflows.md)
  — co-accessibility, dropout modeling, pseudobulk exploration, merging
  studies.
- [Alternative TSS usage and TF
  regulation](https://aifimmunology.github.io/MOCHA/articles/Alternative-TSS-TF-regulation.md)
  — motif enrichment, alternative TSSs, footprinting.
- [Exporting and sharing MOCHA
  results](https://aifimmunology.github.io/MOCHA/articles/MOCHA-export-and-sharing.md)
  —
  [`packMOCHA()`](https://aifimmunology.github.io/MOCHA/reference/packMOCHA.md),
  coverage export, IGV-oriented outputs.
- [Advanced
  modeling](https://aifimmunology.github.io/MOCHA/articles/MOCHA-advanced-modeling.md)
  — mixed models and custom peak training with
  [`trainPeakModel()`](https://aifimmunology.github.io/MOCHA/reference/trainPeakModel.md).

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
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] org.Hs.eg.db_3.22.0                     
#>  [2] TxDb.Hsapiens.UCSC.hg38.knownGene_3.22.0
#>  [3] GenomicFeatures_1.62.0                  
#>  [4] AnnotationDbi_1.72.0                    
#>  [5] BSgenome.Hsapiens.UCSC.hg19_1.4.3       
#>  [6] BSgenome_1.78.0                         
#>  [7] rtracklayer_1.70.1                      
#>  [8] BiocIO_1.20.0                           
#>  [9] Biostrings_2.78.0                       
#> [10] XVector_0.50.0                          
#> [11] MultiAssayExperiment_1.36.1             
#> [12] SummarizedExperiment_1.40.0             
#> [13] Biobase_2.70.0                          
#> [14] GenomicRanges_1.62.1                    
#> [15] Seqinfo_1.0.0                           
#> [16] IRanges_2.44.0                          
#> [17] S4Vectors_0.48.0                        
#> [18] BiocGenerics_0.56.0                     
#> [19] generics_0.1.4                          
#> [20] MatrixGenerics_1.22.0                   
#> [21] matrixStats_1.5.0                       
#> [22] MOCHA_2.0.0                             
#> [23] BiocStyle_2.38.0                        
#> 
#> loaded via a namespace (and not attached):
#>   [1] RcppAnnoy_0.0.23         splines_4.5.3            later_1.4.8             
#>   [4] bitops_1.0-9             RaggedExperiment_1.34.0  tibble_3.3.1            
#>   [7] polyclip_1.10-7          XML_3.99-0.22            fastDummies_1.7.6       
#>  [10] lifecycle_1.0.5          globals_0.19.1           lattice_0.22-9          
#>  [13] MASS_7.3-65              magrittr_2.0.5           plotly_4.12.0           
#>  [16] sass_0.4.10              rmarkdown_2.31           jquerylib_0.1.4         
#>  [19] yaml_2.3.12              httpuv_1.6.17            otel_0.2.0              
#>  [22] Seurat_5.5.0             sctransform_0.4.3        spam_2.11-3             
#>  [25] sp_2.2-1                 spatstat.sparse_3.1-0    reticulate_1.46.0       
#>  [28] cowplot_1.2.0            pbapply_1.7-4            DBI_1.3.0               
#>  [31] RColorBrewer_1.1-3       abind_1.4-8              Rtsne_0.17              
#>  [34] purrr_1.2.2              RCurl_1.98-1.17          ggrepel_0.9.8           
#>  [37] irlba_2.3.7              listenv_0.10.1           spatstat.utils_3.2-3    
#>  [40] goftest_1.2-3            RSpectra_0.16-2          spatstat.random_3.4-5   
#>  [43] fitdistrplus_1.2-6       parallelly_1.47.0        pkgdown_2.2.0           
#>  [46] codetools_0.2-20         DelayedArray_0.36.0      tidyselect_1.2.1        
#>  [49] UCSC.utils_1.6.1         farver_2.1.2             spatstat.explore_3.8-0  
#>  [52] GenomicAlignments_1.46.0 jsonlite_2.0.0           progressr_0.19.0        
#>  [55] ggridges_0.5.7           survival_3.8-6           systemfonts_1.3.2       
#>  [58] tools_4.5.3              ragg_1.5.1               ica_1.0-3               
#>  [61] Rcpp_1.1.1-1.1           glue_1.8.1               gridExtra_2.3           
#>  [64] SparseArray_1.10.8       BiocBaseUtils_1.12.0     xfun_0.57               
#>  [67] GenomeInfoDb_1.46.2      dplyr_1.2.1              withr_3.0.2             
#>  [70] BiocManager_1.30.27      fastmap_1.2.0            digest_0.6.39           
#>  [73] R6_2.6.1                 mime_0.13                textshaping_1.0.5       
#>  [76] scattermore_1.2          tensor_1.5.1             dichromat_2.0-0.1       
#>  [79] spatstat.data_3.1-9      RSQLite_3.52.0           cigarillo_1.0.0         
#>  [82] tidyr_1.3.2              data.table_1.17.8        httr_1.4.8              
#>  [85] htmlwidgets_1.6.4        S4Arrays_1.10.1          uwot_0.2.4              
#>  [88] pkgconfig_2.0.3          gtable_0.3.6             blob_1.3.0              
#>  [91] lmtest_0.9-40            S7_0.2.2                 htmltools_0.5.9         
#>  [94] dotCall64_1.2            bookdown_0.46            plyranges_1.30.1        
#>  [97] SeuratObject_5.4.0       scales_1.4.0             png_0.1-9               
#> [100] spatstat.univar_3.1-7    knitr_1.51               reshape2_1.4.5          
#> [103] rjson_0.2.23             nlme_3.1-169             curl_7.1.0              
#> [106] cachem_1.1.0             zoo_1.8-15               stringr_1.6.0           
#> [109] KernSmooth_2.23-26       parallel_4.5.3           miniUI_0.1.2            
#> [112] restfulr_0.0.16          desc_1.4.3               pillar_1.11.1           
#> [115] grid_4.5.3               vctrs_0.7.3              RANN_2.6.2              
#> [118] promises_1.5.0           xtable_1.8-8             cluster_2.1.8.2         
#> [121] evaluate_1.0.5           cli_3.6.6                compiler_4.5.3          
#> [124] Rsamtools_2.26.0         rlang_1.2.0              crayon_1.5.3            
#> [127] future.apply_1.20.2      plyr_1.8.9               fs_2.1.0                
#> [130] stringi_1.8.7            viridisLite_0.4.3        deldir_2.0-4            
#> [133] BiocParallel_1.44.0      lazyeval_0.2.3           spatstat.geom_3.7-3     
#> [136] Matrix_1.7-5             RcppHNSW_0.6.0           patchwork_1.3.2         
#> [139] bit64_4.8.0              future_1.70.0            ggplot2_4.0.3           
#> [142] KEGGREST_1.50.0          shiny_1.13.0             ROCR_1.0-12             
#> [145] memoise_2.0.1            igraph_2.3.1             bslib_0.11.0            
#> [148] bit_4.6.0
```
