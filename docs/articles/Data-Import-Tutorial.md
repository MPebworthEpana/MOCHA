# Data Import Tutorial: Signac, ArchR, SnapATAC and more

## Introduction

MOCHA peak calling with
[`callOpenTiles()`](https://aifimmunology.github.io/MOCHA/reference/callOpenTiles-methods.md)
runs **after** you have clustered cells and assigned cell-type or
cluster labels. This vignette shows how to reach that step from seven
common inputs:

1.  Bundled `GRangesList` example data (runnable on Bioconductor
    builders).
2.  An [ArchR](https://www.archrproject.com/) `ArchRProject`
    (Hematopoiesis tutorial data).
3.  Direct [Signac](https://stuartlab.org/signac/) / Seurat input
    (indexed fragments + labels).
4.  [Signac](https://stuartlab.org/signac/) 10k PBMC ATAC tutorial
    (Seurat object end-to-end).
5.  Signac multiome — manual fragment extraction (advanced).
6.  [SnapATAC2](https://scverse.org/SnapATAC2/) AnnData object (PBMC 5k
    tutorial).
7.  Legacy SnapATAC (v1) manual extraction from `.snap` files.

For downstream steps (`getSampleTileMatrix`, differentials, motifs), see
the [MOCHA workflow
tutorial](https://aifimmunology.github.io/MOCHA/articles/MOCHA-workflow-tutorial.md).

Reference sections (ArchR, Signac paths, SnapATAC2 shortcut, multiome
manual, and legacy SnapATAC v1) use `eval = FALSE` code chunks so the
vignette builds on Bioconductor without large downloads.

### Prerequisites

| Input | R packages | Other requirements |
|----|----|----|
| Bundled `GRangesList` | MOCHA + annotation `Suggests` | None |
| ArchR | MOCHA, ArchR (install separately) | Pre-clustered `ArchRProject` |
| Direct Seurat/Signac | MOCHA, Seurat, Signac, BSgenome | `ChromatinAssay` with fragment paths; `Sample` + cell-type column |
| Signac 10k ATAC | MOCHA, Seurat, Signac, BSgenome | 10x fragment TSV + index; clustering |
| Signac multiome (manual) | MOCHA, Seurat, Signac, SeuratDisk, EnsDb, BSgenome | Multiome H5 + fragments; label transfer or existing labels |
| SnapATAC2 | MOCHA, reticulate | Python `anndata` or `snapatac2`; clustered `.h5ad` |
| Legacy SnapATAC v1 | MOCHA, SnapATAC (patched fork) | `.snap` files + metadata; large RAM for `extractReads` |

Every path needs **`cellPopLabel`** (cluster/cell-type column),
**`Sample`** (or `sampleColumn`), **`blackList`**, **`genome`**,
**`TxDb`**, and **`OrgDb`**.

## Bundled `GRangesList` (Bioconductor-safe)

The package ships small example objects derived from ArchR `PBMCSmall`
(`exampleFragments`, `exampleCellColData`, `exampleBlackList`).
Fragments use the `RG` cell barcode column; metadata includes `Sample`
and `Clusters`.

Bundled code is maintained in `inst/tutorials/02-import-bundled.R`.

``` r

library(MOCHA)
```

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

``` r

openC2 <- getOpenTiles(tileResults, cellPopulations = "C2")
length(openC2[[1]])
#> [1] 15665
```

## ArchR `ArchRProject`

ArchR is not on Bioconductor; the chunks below are reference code
(`eval = FALSE`). Install ArchR from [the ArchR
website](https://www.archrproject.com/) before running locally.

**Tutorial data:** `ArchR::getTutorialData("Hematopoiesis")` downloads
three 10x fragment files (BMMC, PBMC, CD34-BMMC; ~35k cells total).
Follow the [ArchR
tutorial](https://www.archrproject.com/articles/Articles/tutorial.html)
through clustering before MOCHA.

**Quick local test:** `ArchR::getTestProject()` builds a small
`PBMCSmall` project (see `tests/testthat/TESTING.md` in the MOCHA
repository).

``` r

library(ArchR)
library(MOCHA)

addArchRGenome("hg19")
inputFiles <- ArchR::getTutorialData("Hematopoiesis")

ArrowFiles <- ArchR::createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4,
  filterFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

proj <- ArchR::ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "HemeTutorial",
  copyArrows = TRUE
)
proj <- ArchR::filterDoublets(ArchRProj = proj)
proj <- ArchR::addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI"
)
proj <- ArchR::addClusters(input = proj, reducedDims = "IterativeLSI")
proj <- ArchR::saveArchRProject(ArchRProj = proj)
```

MOCHA reads metadata, blacklist, and genome from the project. You supply
TxDb/OrgDb and optionally restrict populations:

``` r

# proj <- ArchR::loadArchRProject("HemeTutorial")
# Or for a quick test: proj <- ArchR::loadArchRProject("PBMCSmall")

tileResults <- callOpenTiles(
  ATACFragments = proj,
  cellPopLabel = "Clusters",
  cellPopulations = c("C2", "C5"),
  TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
  OrgDb = "org.Hs.eg.db",
  outDir = NULL,
  numCores = 4
)
```

## Direct Seurat / Signac

If your Seurat object already has a labeled `ChromatinAssay` with
fragment paths (standard Signac workflow), you can call MOCHA directly
without manually building a `GRangesList`:

``` r

library(MOCHA)
library(Signac)

# seurat_obj: clustered, labeled Seurat object with ATAC assay + Fragment paths
# Ensure metadata includes Sample and your cell type column (e.g. predicted.id)

tileResults <- MOCHA::callOpenTiles(
  ATACFragments = seurat_obj,
  blackList = blacklist_hg38_unified,
  genome = "BSgenome.Hsapiens.UCSC.hg38",
  cellPopLabel = "predicted.id",
  cellPopulations = "ALL",
  TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
  OrgDb = "org.Hs.eg.db",
  outDir = file.path(getwd(), "MOCHA_Out"),
  numCores = 4
)
```

To inspect or debug the converted inputs before peak calling:

``` r

inputs <- MOCHA::seuratToMOCHAInputs(
  seuratObj = seurat_obj,
  cellPopLabel = "predicted.id"
)
# inputs$ATACFragments  # sample-level GRangesList
# inputs$cellColData    # cell metadata
```

The manual fragment extraction workflow below remains available for
advanced use cases. \## Signac 10k PBMC ATAC {#signac-atac}

This section follows the Signac [PBMC 10k ATAC
vignette](https://stuartlab.org/signac/articles/pbmc_vignette.html)
through clustering, then calls MOCHA directly on the Seurat object. Code
is reference-only (`eval = FALSE`) because of large downloads and
optional packages.

``` r

# 10x PBMC 10k ATACv2 (Signac tutorial URLs)
download.file(
  "https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.h5",
  destfile = "filtered_peak_bc_matrix.h5"
)
download.file(
  "https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_singlecell.csv",
  destfile = "singlecell.csv"
)
download.file(
  "https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz",
  destfile = "fragments.tsv.gz"
)
download.file(
  "https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz.tbi",
  destfile = "fragments.tsv.gz.tbi"
)
```

``` r

library(Signac)
library(Seurat)
library(MOCHA)

set.seed(1234)
counts <- Read10X_h5(filename = "filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "singlecell.csv",
  header = TRUE,
  row.names = 1
)
metadata <- metadata[colnames(counts), , drop = FALSE]

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = "fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

pbmc <- NucleosomeSignal(object = pbmc)
pbmc <- TSSEnrichment(object = pbmc)
pbmc <- subset(
  x = pbmc,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    TSS.enrichment > 2 &
    nucleosome_signal < 4
)

pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q0")
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(object = pbmc, reduction = "lsi", dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = "lsi", dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3, resolution = 0.5)

# Single 10x library: MOCHA still needs a Sample column
pbmc$Sample <- "pbmc10k"
```

**Recommended:** call MOCHA on the Seurat object (no manual
`GRangesList`):

``` r

tileResults <- callOpenTiles(
  ATACFragments = pbmc,
  blackList = blacklist_hg38_unified,
  genome = "BSgenome.Hsapiens.UCSC.hg38",
  cellPopLabel = "seurat_clusters",
  cellPopulations = c("0", "1"),
  sampleColumn = "Sample",
  TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
  OrgDb = "org.Hs.eg.db",
  outDir = file.path(tempdir(), "MOCHA_Signac"),
  numCores = 4
)
```

Inspect converted inputs before peak calling:

``` r

inputs <- seuratToMOCHAInputs(
  seuratObj = pbmc,
  cellPopLabel = "seurat_clusters",
  sampleColumn = "Sample"
)
names(inputs$ATACFragments)
head(inputs$cellColData)
```

## Signac multiome — manual fragment extraction

For the lighter **10k PBMC ATAC-only** tutorial and direct
`callOpenTiles(Seurat)`, see [Data Import Tutorial: Signac 10k PBMC
ATAC](#signac-atac).

This section follows the [Signac multiome
tutorial](https://stuartlab.org/signac/articles/pbmc_multiomic.html) to
build a labeled object and manually extract fragments. If you already
have a Signac object with cell types labelled, skip to [Extract
Fragments from Signac](#extract-frags-signac).

``` r

library(Signac)
library(Seurat)
library(SeuratDisk)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1234)
```

### Signac Tutorial Through Clustering

``` r

# Download files
system('wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5')
system('wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz')
system('wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi')

# Load in ATAC and RNA data
counts <- Read10X_h5("pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
fragpath <- "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
```

``` r

# create a Seurat object containing the RNA adata
pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"

pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

DefaultAssay(pbmc) <- "ATAC"
pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)

pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)

# Transform RNA 
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc)
pbmc <- RunPCA(pbmc)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
```

Label cell types:

``` r

# load PBMC reference
system('wget https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat')
reference <- LoadH5Seurat("pbmc_multimodal.h5seurat")
```

``` r

DefaultAssay(pbmc) <- "SCT"

# transfer cell type labels from reference to query
transfer_anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc,
  normalization.method = "SCT",
  reference.reduction = "spca",
  recompute.residuals = FALSE,
  dims = 1:50
)

predictions <- TransferData(
  anchorset = transfer_anchors, 
  refdata = reference$celltype.l2,
  weight.reduction = pbmc[['pca']],
  dims = 1:50
)

pbmc <- AddMetaData(
  object = pbmc,
  metadata = predictions
)
```

``` r

# set the cell identities to the cell type predictions
Idents(pbmc) <- "predicted.id"

# set a reasonable order for cell types to be displayed when plotting
levels(pbmc) <- c("CD4 Naive", "CD4 TCM", "CD4 CTL", "CD4 TEM", "CD4 Proliferating",
                  "CD8 Naive", "dnT",
                 "CD8 TEM", "CD8 TCM", "CD8 Proliferating", "MAIT", "NK", "NK_CD56bright",
                 "NK Proliferating", "gdT",
                 "Treg", "B naive", "B intermediate", "B memory", "Plasmablast",
                 "CD14 Mono", "CD16 Mono",
                 "cDC1", "cDC2", "pDC", "HSPC", "Eryth", "ASDC", "ILC", "Platelet")

saveRDS(pbmc, 'pbmc_Signac_tutorial.rds')
```

### Extract Fragments from Signac Object

``` r

pbmc <- readRDS('pbmc_Signac_tutorial.rds')

DefaultAssay(pbmc) <- 'ATAC'
fragObj <- Fragments(pbmc)
```

Signac’s training dataset only has one sample, so to simulate multiple
samples, we will duplicate the data here. ***This should not be
necessary with a large dataset.***

``` r

full_pbmc <- merge(
  pbmc, 
  y = c(pbmc, pbmc), 
  add.cell.ids = c('Sample1', 'Sample2', 'Sample3'),
  project = 'DuplicateData'
) 
```

For a larger dataset, you would interact over the fragObj list and
extract each fragments.tsv.gz file.

Instead, we’re just duplicating data for the sake of this tutorial.

More importantly, you do need to modify the barcode in the fragment file
so it matches the barcodes in Seurat’s metadata.

``` r

fragList <- parallel::mclapply(1:3, function(x){
    frags <- read.table(GetFragmentData(fragObj[[1]]))
    names(frags) <- c('chr', 'start', 'end', 'barcode', 'val')
    frags$barcode <- paste("Sample",x,"_", frags$barcode, sep ='')
    frags
}, mc.cores = 3)
names(fragList) <- c('Sample1', 'Sample2', 'Sample3')
```

### Format metadata and fragments for MOCHA

Generate your sample/cell type list by finding all combinations of
Samples and Cell Populations

Here we must also rename the GRangesList to have names in the format
`CellPopulation#Sample`.

``` r

celltype_sample_list <- apply(expand.grid(unique(pbmc@meta.data$predicted.id), names(fragList)), 1, paste, collapse = "#")
# Extract metadata, add the Sample column, as well as CellBarcode. 
fullMeta <- full_pbmc@meta.data
fullMeta$Sample = gsub("_.*","", rownames(full_pbmc@meta.data))
fullMeta$CellBarcode = gsub(".*_","", rownames(full_pbmc@meta.data))

# Change the CellPopulation_Sample format to CellPopulation#Sample for MOCHA
rownames(fullMeta) <- gsub("_","#", rownames(fullMeta))
CellType_GRanges <- pbapply::pblapply(cl = 30, celltype_sample_list, function(x){
  celltype <- gsub("#.*","", x)
  sample <- gsub(".*#","", x)
  sortedMeta <- dplyr::filter(fullMeta, predicted.id == celltype, Sample == sample)
  sortedFrags <- dplyr::filter(fragList[[sample]], barcode %in% rownames(sortedMeta))
  makeGRangesFromDataFrame(sortedFrags, keep.extra.columns = TRUE)
})

names(CellType_GRanges) <- celltype_sample_list
```

Calculate the study signal (the median number of fragments per cell)

``` r

avg_reads <- lapply(fragList, function(x){
                filtFrag <- dplyr::filter(as.data.frame(x), barcode  %in% rownames(fullMeta))
                as.vector(table(filtFrag$barcode))
    })
studySignal <- median(unlist(avg_reads))
```

### Call Open Tiles with MOCHA

``` r

# Our blacklist comes included with Signac
blackList <- blacklist_hg38_unified
# Call Open Tiles
tileResults <- MOCHA::callOpenTiles(ATACFragments = CellType_GRanges,
                             cellColData = fullMeta, 
                             blackList = blackList, 
                             genome = 'BSgenome.Hsapiens.UCSC.hg38',
                             cellPopLabel = 'predicted.id',
                             cellPopulations = fullMeta$predicted.id,
                             studySignal = studySignal,
                             cellCol = 'barcode',
                             TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
                             OrgDb = "org.Hs.eg.db", 
                            outDir = paste(getwd(),'/MOCHA_Out', sep = ''), numCores= 35)

TSAM <- MOCHA::getSampleTileMatrix(tileResults, threshold = 0.2, numCores = 3, verbose = TRUE)
```

## SnapATAC2

SnapATAC2 examples use the [PBMC 5k
tutorial](https://scverse.org/SnapATAC2/tutorials/pbmc.html). MOCHA
imports fragment matrices from `.h5ad` via
[`import_snap_atac()`](https://aifimmunology.github.io/MOCHA/reference/import_snap_atac.md)
(requires **reticulate** and Python **anndata**). Reference code only
(`eval = FALSE`).

``` r

# In R:
if (!requireNamespace("reticulate", quietly = TRUE)) {
  stop("Install reticulate, then configure a Python env with anndata or snapatac2.")
}
reticulate::py_install(c("anndata", "snapatac2"), pip = TRUE)
```

**Shortcut:** use the packaged annotated dataset (includes Leiden
clusters):

``` r

library(MOCHA)
library(reticulate)

snapatac2 <- reticulate::import("snapatac2")
h5ad_path <- snapatac2$datasets$pbmc5k(type = "annotated_h5ad")
adata <- reticulate::import("anndata")$read_h5ad(h5ad_path)

# PBMC 5k is one library; add a Sample column for MOCHA
adata$obs["Sample"] <- "pbmc5k"
```

Or run the full tutorial through Leiden clustering and save `pbmc.h5ad`,
then load it in R with `import_snap_atac("pbmc.h5ad", ...)`.

``` r

inputs <- import_snap_atac(
  adata,
  cellPopLabel = "leiden",
  sampleColumn = "Sample",
  chromPrefix = "add"
)

black_list <- Signac::blacklist_hg38_unified

tileResults <- callOpenTiles(
  ATACFragments = inputs$ATACFragments,
  cellColData = inputs$cellColData,
  blackList = black_list,
  genome = "BSgenome.Hsapiens.UCSC.hg38",
  cellPopLabel = "leiden",
  cellCol = inputs$cellCol,
  cellPopulations = c("0", "1"),
  TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
  OrgDb = "org.Hs.eg.db",
  outDir = file.path(tempdir(), "MOCHA_SnapATAC2"),
  numCores = 4
)
```

[`import_snap_atac()`](https://aifimmunology.github.io/MOCHA/reference/import_snap_atac.md)
does not infer genome or blacklist; harmonize chromosome names with
`chromPrefix` if needed (`"add"` prepends `chr` when missing).

## Legacy SnapATAC (v1)

In this example, we follow the [SnapATAC 10X PBMC
tutorial](https://github.com/r3fang/SnapATAC/tree/master/examples/10X_PBMC_15K#data_download)
through the clustering step before extracting the fragments and cell
metadata necessary for MOCHA.

If you have a fully-formed Snap file for your analysis with clustering
results and sample information added to the metadata, skip to section
[Formatting Snap File Metadata](#format-metadata). If following the
vignette, we STRONGLY recommending installing this patched version of
SnapATAC from [this repository](https://github.com/imran-aifi/SnapATAC)
with `devtools::install_github("imran-aifi/SnapATAC")`.

### SnapATAC Tutorial through Clustering

Download the all the SnapATAC Tutorial Data (linked above) to your
working directory, and load it:

``` bash
# CLI tool for installing from Google Drive share links
pip install gdown
# https://drive.google.com/file/d/1YiYd_Ydes3tqsJGpNuqQquUOoVj2EEjE/view?usp=share_link
gdown https://drive.google.com/uc?id=1YiYd_Ydes3tqsJGpNuqQquUOoVj2EEjE
# https://drive.google.com/file/d/1NvGn4M2_HD06PL5Nj2if5xO-uVA0y8Q5/view?usp=share_link
gdown https://drive.google.com/uc?id=1NvGn4M2_HD06PL5Nj2if5xO-uVA0y8Q5
# https://drive.google.com/file/d/1LUOqsXoQN6lVx-4RNlgH90e5oXQ0y9Bd/view?usp=share_link
gdown https://drive.google.com/uc?id=1LUOqsXoQN6lVx-4RNlgH90e5oXQ0y9Bd
# https://drive.google.com/file/d/1oMJ6wFsfS-q-sY_yaLEtnYebM7RrBix6/view?usp=share_link
gdown https://drive.google.com/uc?id=1oMJ6wFsfS-q-sY_yaLEtnYebM7RrBix6
# https://drive.google.com/file/d/1SEFZ5CJgcmoAkmo4_1kCMYS60YOFP379/view?usp=share_link
gdown https://drive.google.com/uc?id=1SEFZ5CJgcmoAkmo4_1kCMYS60YOFP379
# https://drive.google.com/file/d/1RlBvTCqz6mhaTAkfYiCdeojD-U2mN2wp/view?usp=share_link
gdown https://drive.google.com/uc?id=1RlBvTCqz6mhaTAkfYiCdeojD-U2mN2wp
```

``` r

library(SnapATAC);
snap.files = c(
  "atac_pbmc_5k_nextgem.snap", 
  "atac_pbmc_10k_nextgem.snap"
);
sample.names = c(
  "PBMC 5K",
  "PBMC 10K"
);
barcode.files = c(
  "atac_pbmc_5k_nextgem_singlecell.csv",
  "atac_pbmc_10k_nextgem_singlecell.csv"
);
x.sp.ls = lapply(seq(snap.files), function(i){
  createSnap(
      file=snap.files[i],
      sample=sample.names[i]
  );
})
names(x.sp.ls) = sample.names;
barcode.ls = lapply(seq(snap.files), function(i){
  barcodes = read.csv(
      barcode.files[i], 
      head=TRUE
  );
  barcodes = barcodes[2:nrow(barcodes),];
  barcodes$logUMI = log10(barcodes$passed_filters + 1);
  barcodes$promoter_ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1);
  barcodes
})
x.sp.ls
```

``` r

# for both datasets, we identify usable barcodes using [3.5-5] for log10(UMI) and [0.4-0.8] for promoter ratio as cutoff.
cutoff.logUMI.low = c(3.5, 3.5);
cutoff.logUMI.high = c(5, 5);
cutoff.FRIP.low = c(0.4, 0.4);
cutoff.FRIP.high = c(0.8, 0.8);
barcode.ls = lapply(seq(snap.files), function(i){
  barcodes = barcode.ls[[i]];
  idx = which(
      barcodes$logUMI >= cutoff.logUMI.low[i] & 
      barcodes$logUMI <= cutoff.logUMI.high[i] & 
      barcodes$promoter_ratio >= cutoff.FRIP.low[i] &
      barcodes$promoter_ratio <= cutoff.FRIP.high[i]
  );
  barcodes[idx,]
});
x.sp.ls = lapply(seq(snap.files), function(i){
  barcodes = barcode.ls[[i]];
  x.sp = x.sp.ls[[i]];
  barcode.shared = intersect(x.sp@barcode, barcodes$barcode);
  x.sp = x.sp[match(barcode.shared, x.sp@barcode),];
  barcodes = barcodes[match(barcode.shared, barcodes$barcode),];
  x.sp@metaData = barcodes;
  x.sp
})
names(x.sp.ls) = sample.names;
x.sp.ls
```

``` r

# combine two snap object
x.sp = Reduce(snapRbind, x.sp.ls);
x.sp@metaData["Sample"] = x.sp@sample;
print(table(x.sp@sample))
x.sp
```

``` r

# Step 2. Add cell-by-bin matrix
x.sp = addBmatToSnap(x.sp, bin.size=5000);
# Step 3. Matrix binarization
x.sp = makeBinary(x.sp, mat="bmat");
# Step 4. Bin filtering
library(GenomicRanges);
black_list = read.table("hg19.blacklist.bed.gz");
black_list.gr = GRanges(
  black_list[,1], 
  IRanges(black_list[,2], black_list[,3])
);
idy = queryHits(
  findOverlaps(x.sp@feature, black_list.gr)
);
if(length(idy) > 0){
  x.sp = x.sp[,-idy, mat="bmat"];
};
x.sp

# Remove unwanted chromosomes
chr.exclude = seqlevels(x.sp@feature)[grep("random|chrM", seqlevels(x.sp@feature))];
idy = grep(paste(chr.exclude, collapse="|"), x.sp@feature);
if(length(idy) > 0){
  x.sp = x.sp[,-idy, mat="bmat"]
};
x.sp
```

``` r

# The coverage of bins roughly obeys a log normal distribution. We remove the top 5% bins that overlap with invariant features such as the house keeping gene promoters.
bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);
hist(
  bin.cov[bin.cov > 0], 
  xlab="log10(bin cov)", 
  main="log10(Bin Cov)", 
  col="lightblue", 
  xlim=c(0, 5)
);
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95);
idy = which(bin.cov <= bin.cutoff & bin.cov > 0);
x.sp = x.sp[, idy, mat="bmat"];
x.sp

# We will further remove any cells of bin coverage less than 1,000. The rational behind this is that some cells may have high number of unique fragments but end up with low bin coverage after filtering. This step is optional but highly recommended.
idx = which(Matrix::rowSums(x.sp@bmat) > 1000);
x.sp = x.sp[idx,];
x.sp
```

``` r

# Step 5. Dimensionality reduction
row.covs.dens <- density(
  x = x.sp@metaData[,"logUMI"], 
  bw = 'nrd', adjust = 1
);
sampling_prob <- 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = x.sp@metaData[,"logUMI"])$y + .Machine$double.eps);
set.seed(1);
idx.landmark.ds <- base::sort(sample(x = seq(nrow(x.sp)), size = 10000, prob = sampling_prob));
x.landmark.sp = x.sp[idx.landmark.ds,];
x.query.sp = x.sp[-idx.landmark.ds,];
x.landmark.sp = runDiffusionMaps(
  obj= x.landmark.sp,
  input.mat="bmat", 
  num.eigs=50
);
x.landmark.sp@metaData$landmark = 1;
x.query.sp = runDiffusionMapsExtension(
  obj1=x.landmark.sp, 
  obj2=x.query.sp,
  input.mat="bmat"
);
x.query.sp@metaData$landmark = 0;
x.sp = snapRbind(x.landmark.sp, x.query.sp);
x.sp = x.sp[order(x.sp@metaData["sample"])];
x.sp = runKNN(
  obj=x.sp,
  eigs.dims=1:20,
  k=15
);
x.sp=runCluster(
  obj=x.sp,
  tmp.folder=tempdir(),
  louvain.lib="R-igraph", #"leiden" preferred, but may cause issues. Requires 'library(leiden)'.
  seed.use=10,
  resolution=0.7
)
```

### Format Snap File Metadata

Add the computed clusters to the Snap object metadata.

The Snap object contains two samples, “PBMC 5K” and “PBMC 10K”. Let’s
add a “Sample” column to the metadata. Let’s also add a column “files”
pointing to the original .snap files from which each cell came.

``` r

# Add clusters (from SnapATAC::runCluster) to metadata
x.sp@metaData$cluster = x.sp@cluster

# Add Sample name to metadata (if not done previously)
x.sp@metaData$Sample = x.sp@sample

# Add files to metadata, indicating the original snap file each cell belongs to.
snap.files <- c(
  "atac_pbmc_5k_nextgem.snap", 
  "atac_pbmc_10k_nextgem.snap"
)
fileList <- unlist(lapply(x.sp@metaData$Sample, function(x){
  ifelse(x == "PBMC 5K", snap.files[[1]],snap.files[[2]])
}))
x.sp@metaData$files <- fileList

# SAVE this metadata to disk
write.csv(x.sp@metaData, "./snapMetadataforMOCHA.csv")
```

### Extract Fragments from Snap File

Now we have a Snap object with metadata containing barcodes, unique cell
ids (column cell_id), sample names, and cell populations (cluster). We
also have our HG19 blackList, `black_list.gr`.

> Note: Following the tutorial from SnapATAC can often result in a
> segfault when extracting fragments with `SnapATAC::extractReads`. We
> recommend running on a machine with large RAM and avoiding
> parallelization.

Next we extract reads by sample and cell population, ensuring our final
GRanges list is named following the pattern `CellPopulation#Sample`.

``` r

snapMetadata <- read.csv("./snapMetadataforMOCHA.csv")
cellCol <- "barcode"
cellPopLabel <- "cluster"
cellPopulations <- unique(snapMetadata$cluster)
allSamples <- unique(snapMetadata$Sample)

fragmentsGRangesList <- unlist(lapply(allSamples, function(sample){
  barcodesList <- lapply(cellPopulations, function(cellPop) {
    snapMetadata[snapMetadata$Sample == sample,]
    # Extract barcodes for a single cell population
    cellPopIdx <- snapMetadata$cluster == cellPop
    cellPopBarcodes <- snapMetadata[cellPopIdx,]$barcode
    
    # Build the file list for the selected cell barcode
    files <- snapMetadata[cellPopIdx,]$files
    
    # Extract fragments
    cellPopFrags <- SnapATAC::extractReads(cellPopBarcodes, files, do.par = FALSE)
  })
  names(barcodesList) <- paste(cellPopulations, sample, sep="#")
  barcodesList
}))
```

Calculate the study signal (the median number of fragments per cell)

``` r

avg_reads <- lapply(fragmentsGRangesList, function(x){
  filtFrag <- dplyr::filter(as.data.frame(x), barcode  %in% snapMetadata$barcode)
  as.vector(table(filtFrag$barcode))
})
studySignal <- median(unlist(avg_reads))
```

### Call Open Tiles with MOCHA

``` r

# Call Open Tiles
tileResults <- MOCHA::callOpenTiles(
  ATACFragments = fragmentsGRangesList,
  cellColData = snapMetadata, 
  blackList = black_list.gr, 
  genome = "BSgenome.Hsapiens.UCSC.hg38",
  cellPopLabel = cellPopLabel,
  cellPopulations = cellPopulations,
  studySignal = studySignal,
  cellCol = cellCol,
  TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
  OrgDb = "org.Hs.eg.db", 
  outDir = paste(getwd(),'/MOCHA_Out', sep = ''), 
  numCores = 5
)

TSAM <- MOCHA::getSampleTileMatrix(tileResults, threshold = 0.2, numCores = 3, verbose = TRUE)
```

## Utilities after import

After
[`callOpenTiles()`](https://aifimmunology.github.io/MOCHA/reference/callOpenTiles-methods.md),
these helpers inspect or adjust objects without running a full
downstream workflow:

``` r

# Per-sample fragment GRanges from an ArchR project (requires ArchR):
# frags <- getPopFrags(ArchRProj, cellPopLabel = "Clusters")

# Quick access to called peaks from tileResults
openC2_util <- getOpenTiles(tileResults, cellPopulations = "C2")
# equivalent: getOpenTiles(tileResults, cellPopulations = "C2")
```

See
[`?getPopFrags`](https://aifimmunology.github.io/MOCHA/reference/getPopFrags.md),
[`?correctGenome`](https://aifimmunology.github.io/MOCHA/reference/correctGenome.md),
[`?openTiles`](https://aifimmunology.github.io/MOCHA/reference/openTiles.md),
and
[`?getOpenTiles`](https://aifimmunology.github.io/MOCHA/reference/getOpenTiles.md).

## Next steps

After
[`callOpenTiles()`](https://aifimmunology.github.io/MOCHA/reference/callOpenTiles-methods.md)
returns a `MultiAssayExperiment`:

- Build sample-tile matrices with
  [`getSampleTileMatrix()`](https://aifimmunology.github.io/MOCHA/reference/getSampleTileMatrix.md).
- Continue with annotation, differentials, and motifs in the [MOCHA
  workflow
  tutorial](https://aifimmunology.github.io/MOCHA/articles/MOCHA-workflow-tutorial.md).

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
#>  [5] Biobase_2.70.0                          
#>  [6] BSgenome.Hsapiens.UCSC.hg19_1.4.3       
#>  [7] BSgenome_1.78.0                         
#>  [8] rtracklayer_1.70.1                      
#>  [9] BiocIO_1.20.0                           
#> [10] Biostrings_2.78.0                       
#> [11] XVector_0.50.0                          
#> [12] GenomicRanges_1.62.1                    
#> [13] Seqinfo_1.0.0                           
#> [14] IRanges_2.44.0                          
#> [15] S4Vectors_0.48.0                        
#> [16] BiocGenerics_0.56.0                     
#> [17] generics_0.1.4                          
#> [18] MOCHA_2.0.0                             
#> [19] BiocStyle_2.38.0                        
#> 
#> loaded via a namespace (and not attached):
#>   [1] RcppAnnoy_0.0.23            splines_4.5.3              
#>   [3] later_1.4.8                 bitops_1.0-9               
#>   [5] RaggedExperiment_1.34.0     tibble_3.3.1               
#>   [7] polyclip_1.10-7             XML_3.99-0.22              
#>   [9] fastDummies_1.7.6           lifecycle_1.0.5            
#>  [11] globals_0.19.1              lattice_0.22-9             
#>  [13] MASS_7.3-65                 MultiAssayExperiment_1.36.1
#>  [15] magrittr_2.0.5              plotly_4.12.0              
#>  [17] sass_0.4.10                 rmarkdown_2.31             
#>  [19] jquerylib_0.1.4             yaml_2.3.12                
#>  [21] httpuv_1.6.17               otel_0.2.0                 
#>  [23] Seurat_5.5.0                sctransform_0.4.3          
#>  [25] spam_2.11-3                 sp_2.2-1                   
#>  [27] spatstat.sparse_3.1-0       reticulate_1.46.0          
#>  [29] cowplot_1.2.0               pbapply_1.7-4              
#>  [31] DBI_1.3.0                   RColorBrewer_1.1-3         
#>  [33] abind_1.4-8                 Rtsne_0.17                 
#>  [35] purrr_1.2.2                 RCurl_1.98-1.17            
#>  [37] ggrepel_0.9.8               irlba_2.3.7                
#>  [39] listenv_0.10.1              spatstat.utils_3.2-3       
#>  [41] goftest_1.2-3               RSpectra_0.16-2            
#>  [43] spatstat.random_3.4-5       fitdistrplus_1.2-6         
#>  [45] parallelly_1.47.0           pkgdown_2.2.0              
#>  [47] codetools_0.2-20            DelayedArray_0.36.0        
#>  [49] tidyselect_1.2.1            UCSC.utils_1.6.1           
#>  [51] farver_2.1.2                matrixStats_1.5.0          
#>  [53] spatstat.explore_3.8-0      GenomicAlignments_1.46.0   
#>  [55] jsonlite_2.0.0              progressr_0.19.0           
#>  [57] ggridges_0.5.7              survival_3.8-6             
#>  [59] systemfonts_1.3.2           tools_4.5.3                
#>  [61] ragg_1.5.1                  ica_1.0-3                  
#>  [63] Rcpp_1.1.1-1.1              glue_1.8.1                 
#>  [65] gridExtra_2.3               SparseArray_1.10.8         
#>  [67] BiocBaseUtils_1.12.0        xfun_0.57                  
#>  [69] MatrixGenerics_1.22.0       GenomeInfoDb_1.46.2        
#>  [71] dplyr_1.2.1                 withr_3.0.2                
#>  [73] BiocManager_1.30.27         fastmap_1.2.0              
#>  [75] digest_0.6.39               R6_2.6.1                   
#>  [77] mime_0.13                   textshaping_1.0.5          
#>  [79] scattermore_1.2             tensor_1.5.1               
#>  [81] dichromat_2.0-0.1           spatstat.data_3.1-9        
#>  [83] RSQLite_3.52.0              cigarillo_1.0.0            
#>  [85] tidyr_1.3.2                 data.table_1.17.8          
#>  [87] httr_1.4.8                  htmlwidgets_1.6.4          
#>  [89] S4Arrays_1.10.1             uwot_0.2.4                 
#>  [91] pkgconfig_2.0.3             gtable_0.3.6               
#>  [93] blob_1.3.0                  lmtest_0.9-40              
#>  [95] S7_0.2.2                    htmltools_0.5.9            
#>  [97] dotCall64_1.2               bookdown_0.46              
#>  [99] plyranges_1.30.1            SeuratObject_5.4.0         
#> [101] scales_1.4.0                png_0.1-9                  
#> [103] spatstat.univar_3.1-7       knitr_1.51                 
#> [105] reshape2_1.4.5              rjson_0.2.23               
#> [107] nlme_3.1-169                curl_7.1.0                 
#> [109] cachem_1.1.0                zoo_1.8-15                 
#> [111] stringr_1.6.0               KernSmooth_2.23-26         
#> [113] parallel_4.5.3              miniUI_0.1.2               
#> [115] restfulr_0.0.16             desc_1.4.3                 
#> [117] pillar_1.11.1               grid_4.5.3                 
#> [119] vctrs_0.7.3                 RANN_2.6.2                 
#> [121] promises_1.5.0              xtable_1.8-8               
#> [123] cluster_2.1.8.2             evaluate_1.0.5             
#> [125] cli_3.6.6                   compiler_4.5.3             
#> [127] Rsamtools_2.26.0            rlang_1.2.0                
#> [129] crayon_1.5.3                future.apply_1.20.2        
#> [131] plyr_1.8.9                  fs_2.1.0                   
#> [133] stringi_1.8.7               deldir_2.0-4               
#> [135] viridisLite_0.4.3           BiocParallel_1.44.0        
#> [137] lazyeval_0.2.3              spatstat.geom_3.7-3        
#> [139] Matrix_1.7-5                RcppHNSW_0.6.0             
#> [141] patchwork_1.3.2             bit64_4.8.0                
#> [143] future_1.70.0               ggplot2_4.0.3              
#> [145] KEGGREST_1.50.0             shiny_1.13.0               
#> [147] SummarizedExperiment_1.40.0 ROCR_1.0-12                
#> [149] igraph_2.3.1                memoise_2.0.1              
#> [151] bslib_0.11.0                bit_4.6.0
```
