# Perform peak-calling on a set of fragments or an ArchR Project

`callOpenTiles` is the main peak-calling function in MOCHA that serves
as a wrapper function to call peaks provided a set of fragment files and
an ArchR Project for meta-data purposes

## Usage

``` r
callOpenTiles(
  ATACFragments,
  cellColData,
  blackList,
  genome,
  cellPopLabel,
  cellPopulations = "ALL",
  sampleColumn = "Sample",
  studySignal = NULL,
  generalizeStudySignal = FALSE,
  cellCol = "RG",
  TxDb,
  OrgDb,
  outDir,
  numCores = 30,
  verbose = FALSE,
  force = FALSE,
  peakModel = NULL,
  returnClass = c("legacy", "mocha")
)

.callOpenTiles_default(
  ATACFragments,
  cellColData,
  blackList,
  genome,
  cellPopLabel,
  cellPopulations = "ALL",
  sampleColumn = "Sample",
  studySignal = NULL,
  generalizeStudySignal = FALSE,
  cellCol = "RG",
  TxDb,
  OrgDb,
  outDir,
  numCores = 30,
  verbose = FALSE,
  force = FALSE,
  peakModel = NULL,
  returnClass = c("legacy", "mocha")
)

# S4 method for class 'GRangesList'
callOpenTiles(
  ATACFragments,
  cellColData,
  blackList,
  genome,
  cellPopLabel,
  cellPopulations = "ALL",
  sampleColumn = "Sample",
  studySignal = NULL,
  generalizeStudySignal = FALSE,
  cellCol = "RG",
  TxDb,
  OrgDb,
  outDir,
  numCores = 30,
  verbose = FALSE,
  force = FALSE,
  peakModel = NULL,
  returnClass = c("legacy", "mocha")
)

# S4 method for class 'list'
callOpenTiles(
  ATACFragments,
  cellColData,
  blackList,
  genome,
  cellPopLabel,
  cellPopulations = "ALL",
  sampleColumn = "Sample",
  studySignal = NULL,
  generalizeStudySignal = FALSE,
  cellCol = "RG",
  TxDb,
  OrgDb,
  outDir,
  numCores = 30,
  verbose = FALSE,
  force = FALSE,
  peakModel = NULL,
  returnClass = c("legacy", "mocha")
)

.callOpenTiles_ArchR(
  ATACFragments,
  cellPopLabel,
  cellPopulations = "ALL",
  sampleColumn = "Sample",
  studySignal = NULL,
  generalizeStudySignal = FALSE,
  TxDb,
  OrgDb,
  outDir = NULL,
  numCores = 30,
  verbose = FALSE,
  force = FALSE,
  peakModel = NULL,
  returnClass = c("legacy", "mocha")
)
```

## Arguments

- ATACFragments:

  an ArchR Project, a Seurat object with a Signac ChromatinAssay
  (fragment paths in the assay), or a GRangesList of fragments. Two
  input formats are supported for non-ArchR GRangesList input:

  - \*\*Legacy:\*\* each list element is named \`CellPopulation#Sample\`
    (optional \`\_\_normalizationFactor\` suffix), with mixed cell types
    already split per sample.

  - \*\*Sample-level:\*\* each list element is named by sample ID
    (matching \`cellColData\[\[sampleColumn\]\]\`), containing fragments
    from all cell types in that sample. Cell populations are derived
    from \`cellColData\` using \`cellPopLabel\` and \`cellCol\`.

  Each GRanges must contain unique cell IDs in the column given by
  \`cellCol\`.

- cellColData:

  A DataFrame containing cell-level metadata. This must contain both a
  column 'Sample' with unique sample IDs and the column specified by
  'cellPopLabel'.

- blackList:

  A GRanges of blacklisted regions

- genome:

  A BSgenome object, or the full name of an installed BSgenome data
  package, or a short string specifying the name of an NCBI assembly
  (e.g. "GRCh38", "TAIR10.1", etc...) or UCSC genome (e.g. "hg38",
  "bosTau9", "galGal6", "ce11", etc...). The supplied short string must
  refer unambiguously to an installed BSgenome data package. See
  [getBSgenome](https://rdrr.io/pkg/BSgenome/man/available.genomes.html).

- cellPopLabel:

  string indicating which column in the ArchRProject metadata contains
  the cell population label.

- cellPopulations:

  vector of strings. Cell subsets for which to call peaks. This list of
  group names must be identical to names that appear in the ArchRProject
  metadata. Optional, if cellPopulations='ALL', then peak calling is
  done on all cell populations in the ArchR project metadata. Default is
  'ALL'.

- sampleColumn:

  The name of the metadata column with biological sample names. Standard
  is 'Sample', as fixed by ArchR

- studySignal:

  The median signal (number of fragments) in your study. If not set,
  this will be calculated using the input ArchR project but relies on
  the assumption that the ArchR project encompasses your whole study
  (i.e. is not a subset).

- generalizeStudySignal:

  If \`studySignal\` is not provided, calculate the signal as the mean
  of the mean & median number of fragments for of individual samples
  within each cell population. This may improve MOCHA's ability to
  generalize to datasets with XXXXXX \#TODO. Default is FALSE, use the
  median number of fragments.

- cellCol:

  The column in cellColData specifying unique cell ids or barcodes.
  Default is "RG", the unique cell identifier used by ArchR.

- TxDb:

  The exact package name of a TxDb-class transcript annotation package
  for your organism (e.g. "TxDb.Hsapiens.UCSC.hg38.knownGene"). This
  must be installed. See [Bioconductor AnnotationData
  Packages](https://bioconductor.org/packages/release/data/annotation/).

- OrgDb:

  The exact package name of a OrgDb-class genome wide annotation package
  for your organism (e.g. "org.Hs.eg.db"). This must be installed. See
  [Bioconductor AnnotationData
  Packages](https://bioconductor.org/packages/release/data/annotation/)

- outDir:

  is a string describing the output directory for coverage files. Must
  be a complete directory string. With ArchR input, set outDir to NULL
  to create a directory within the input ArchR project directory named
  MOCHA for saving files. Coverage tracks are written under
  \`tracks/cellPop/Accessibility\|Insertions/sample.bw\` with a
  \`tracks/manifest.rds\` index. Legacy \`cellPop_CoverageFiles.RDS\`
  bundles are still read when present.

- numCores:

  integer. Number of cores to parallelize peak-calling across multiple
  cell populations.

- verbose:

  Set TRUE to display additional messages. Default is FALSE.

- force:

  Optional, whether to force creation of coverage files if they already
  exist. Default is FALSE.

- peakModel:

  Optional `MOCHAPeakModel` from
  [`trainPeakModel`](https://aifimmunology.github.io/MOCHA/reference/trainPeakModel.md).
  When `NULL`, the bundled 500 bp model is used. Custom models set the
  tile width and study-signal calibration (`trainingMedian`).

- returnClass:

  Return type: `"legacy"` (default, `MultiAssayExperiment`) or `"mocha"`
  (`MochaTileResults`).

## Value

tileResults A MultiAssayExperiment object containing ranged data for
each tile

## Examples

``` r
# \donttest{
# ArchR Project input requires the ArchR package (not run here).
# Starting from GRangesList with bundled example data:
# Starting from GRangesList:
if (
  requireNamespace("BSgenome.Hsapiens.UCSC.hg19") &&
  requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene") &&
  requireNamespace("org.Hs.eg.db")
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
    cellPopulations = c("C2", "C5"),
    numCores = 1
  )
}
#> Warning: GRanges object contains 2447 out-of-bound ranges located on sequences chr1
#>   and chr2. Note that ranges located on a sequence whose length is unknown (NA)
#>   or on a circular sequence are not considered out-of-bound (use seqlengths()
#>   and isCircular() to get the lengths and circularity flags of the underlying
#>   sequences). You can use trim() to trim these ranges. See
#>   ?`trim,GenomicRanges-method` for more information.
#> Warning: GRanges object contains 2560 out-of-bound ranges located on sequences chr1
#>   and chr2. Note that ranges located on a sequence whose length is unknown (NA)
#>   or on a circular sequence are not considered out-of-bound (use seqlengths()
#>   and isCircular() to get the lengths and circularity flags of the underlying
#>   sequences). You can use trim() to trim these ranges. See
#>   ?`trim,GenomicRanges-method` for more information.
# Sample-level GRangesList (one element per sample):
# sampleFragments <- GRangesList(Sample1 = frags1, Sample2 = frags2)
# callOpenTiles(sampleFragments, cellColData = meta, cellPopLabel = "Clusters", ...)
# }
```
