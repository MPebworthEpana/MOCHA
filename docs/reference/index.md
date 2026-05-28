# Package index

## Open Chromatin Identification

- [`callOpenTiles()`](https://aifimmunology.github.io/MOCHA/reference/callOpenTiles-methods.md)
  [`.callOpenTiles_default()`](https://aifimmunology.github.io/MOCHA/reference/callOpenTiles-methods.md)
  [`.callOpenTiles_ArchR()`](https://aifimmunology.github.io/MOCHA/reference/callOpenTiles-methods.md)
  : Perform peak-calling on a set of fragments or an ArchR Project
- [`getSampleTileMatrix()`](https://aifimmunology.github.io/MOCHA/reference/getSampleTileMatrix.md)
  : Get consensus sample-tile matrices containing the signal intensity
  at each tile
- [`import_snap_atac()`](https://aifimmunology.github.io/MOCHA/reference/import_snap_atac.md)
  : Import SnapATAC2 AnnData / AnnDataSet for MOCHA peak calling
- [`trainPeakModel()`](https://aifimmunology.github.io/MOCHA/reference/trainPeakModel.md)
  : Retrain MOCHA peak-calling models at a chosen tile size

## Plotting

- [`plotConsensus()`](https://aifimmunology.github.io/MOCHA/reference/plotConsensus.md)
  : Plot to determine the reproducibility threshold
- [`plotDropoutDiagnostics()`](https://aifimmunology.github.io/MOCHA/reference/plotDropoutDiagnostics.md)
  : Plot dropout model diagnostics
- [`plotIntensityDistribution()`](https://aifimmunology.github.io/MOCHA/reference/plotIntensityDistribution.md)
  : Plots the distribution of sample-tile intensities for a given cell
  population
- [`plotRegion()`](https://aifimmunology.github.io/MOCHA/reference/plotRegion.md)
  : Plot a given region summarized across all cell groupings
- [`suggestConsensusThreshold()`](https://aifimmunology.github.io/MOCHA/reference/suggestConsensusThreshold.md)
  : Suggest a reproducibility threshold for consensus tile selection

## Downstream Analyses

- [`MotifEnrichment()`](https://aifimmunology.github.io/MOCHA/reference/MotifEnrichment.md)
  : Test for enrichment of motifs against a background

- [`MotifSetEnrichmentAnalysis()`](https://aifimmunology.github.io/MOCHA/reference/MotifSetEnrichmentAnalysis.md)
  : Test the enrichment of a given TF motif set against a motif set
  downstream of multiple ligands

- [`addAccessibilityShift()`](https://aifimmunology.github.io/MOCHA/reference/addAccessibilityShift.md)
  : Add difference in accessibility between two conditions

- [`addInsertionBias()`](https://aifimmunology.github.io/MOCHA/reference/addInsertionBias.md)
  : Add InsertionBias to MOCHA object

- [`assessDropout()`](https://aifimmunology.github.io/MOCHA/reference/assessDropout.md)
  : Assess technical vs biological zeros on a sample-tile matrix

- [`bulkDimReduction()`](https://aifimmunology.github.io/MOCHA/reference/bulkDimReduction.md)
  : Run PCA or LSI dimensionality reduction on tiles

- [`bulkUMAP()`](https://aifimmunology.github.io/MOCHA/reference/bulkUMAP.md)
  : Generate UMAP from pseudobulk LSI results

- [`estimateDropoutModel()`](https://aifimmunology.github.io/MOCHA/reference/estimateDropoutModel.md)
  : Estimate a logistic dropout model for a cell population

- [`filterCoAccessibleLinks()`](https://aifimmunology.github.io/MOCHA/reference/filterCoAccessibleLinks.md)
  : Filter links by correlation strength

- [`getAltTSS()`](https://aifimmunology.github.io/MOCHA/reference/getAltTSS.md)
  : Annotate peaks falling in Transcription Start Sites (TSSs) and
  identify alternatively regulated TSSs for each gene

- [`getCellTypeMotifs()`](https://aifimmunology.github.io/MOCHA/reference/getCellTypeMotifs.md)
  : Get motifset for a given cell type

- [`getCoAccessibleLinks()`](https://aifimmunology.github.io/MOCHA/reference/getCoAccessibleLinks.md)
  : Find co-accessible neighboring regions

- [`getDifferentialAccessibleTiles()`](https://aifimmunology.github.io/MOCHA/reference/getDifferentialAccessibleTiles.md)
  : Conduct a differential test between open regions of two sample
  groups

- [`linearModeling()`](https://aifimmunology.github.io/MOCHA/reference/linearModeling.md)
  :

  `linearModeling`

- [`motifFootprint()`](https://aifimmunology.github.io/MOCHA/reference/motifFootprint.md)
  : Generate motif footprints

- [`pilotLMEM()`](https://aifimmunology.github.io/MOCHA/reference/pilotLMEM.md)
  : Execute a pilot run of single linear model on a subset of data

- [`pilotZIGLMM()`](https://aifimmunology.github.io/MOCHA/reference/pilotZIGLMM.md)
  : Execute a pilot run of model on a subset of data

- [`runLMEM()`](https://aifimmunology.github.io/MOCHA/reference/runLMEM.md)
  : Run Linear Mixed-Effects Modeling for continuous, non-zero inflated
  data

- [`runZIGLMM()`](https://aifimmunology.github.io/MOCHA/reference/runZIGLMM.md)
  : Run Zero-inflated Generalized Linear Mixed Modeling on pseudobulked
  scATAC data

- [`testCoAccessibility()`](https://aifimmunology.github.io/MOCHA/reference/testCoAccessibility.md)
  : Test if tile pairs are significantly different against a random,
  non-overlapping background set

- [`testCoAccessibilityChromVar()`](https://aifimmunology.github.io/MOCHA/reference/testCoAccessibilityChromVar.md)
  : Test input tile pairs against a ChromVAR background

- [`testCoAccessibilityRandom()`](https://aifimmunology.github.io/MOCHA/reference/testCoAccessibilityRandom.md)
  : Test input tile pairs against a random background

- [`varZIGLMM()`](https://aifimmunology.github.io/MOCHA/reference/varZIGLMM.md)
  : Zero-inflated Variance Decomposition for pseudobulked scATAC data

## Exporting/Sharing

- [`exportCoverage()`](https://aifimmunology.github.io/MOCHA/reference/exportCoverage.md)
  : Export normalized coverage files to per-sample or sample-averaged
  (per cell population) BigWig files.

- [`exportDifferentials()`](https://aifimmunology.github.io/MOCHA/reference/exportDifferentials.md)
  :

  Export differential peaks from
  [`getDifferentialAccessibleTiles()`](https://aifimmunology.github.io/MOCHA/reference/getDifferentialAccessibleTiles.md)
  to BigBed format for visualization in genome browsers.

- [`exportLocalFootprints()`](https://aifimmunology.github.io/MOCHA/reference/exportLocalFootprints.md)
  : Export insertion counts to per-sample BigWig files after applying a
  rolling sum and rolling median smoothing filter.

- [`exportMotifs()`](https://aifimmunology.github.io/MOCHA/reference/exportMotifs.md)
  :

  Export motif annotations from
  [`addMotifSet()`](https://aifimmunology.github.io/MOCHA/reference/addMotifSet.md)
  to BigBed format for visualization in genome browsers.

- [`exportOpenTiles()`](https://aifimmunology.github.io/MOCHA/reference/exportOpenTiles.md)
  : Export open tiles of a given cell population to BigBed format for
  visualization in genome browsers.

- [`packMOCHA()`](https://aifimmunology.github.io/MOCHA/reference/packMOCHA.md)
  : Zip up a MOCHA object and its linked files for portability between
  file systems

- [`plotMotifs()`](https://aifimmunology.github.io/MOCHA/reference/plotMotifs.md)
  : plotMotifs - plots motif footprints, exported from motifFootprint

- [`unpackMOCHA()`](https://aifimmunology.github.io/MOCHA/reference/unpackMOCHA.md)
  :

  Unzip a MOCHA object and its linked files created by
  [`packMOCHA()`](https://aifimmunology.github.io/MOCHA/reference/packMOCHA.md)
  for portability between file systems

## Utilities

- [`GRangesToString()`](https://aifimmunology.github.io/MOCHA/reference/GRangesToString.md)
  : Convert a GRanges object to a string in the format 'chr1:100-200'

- [`StringsToGRanges()`](https://aifimmunology.github.io/MOCHA/reference/StringsToGRanges.md)
  : Convert a list of strings in the format "chr1:100-200" into a
  GRanges

- [`addCellColData()`](https://aifimmunology.github.io/MOCHA/reference/addCellColData.md)
  : Add a column to the sample-level colData of a MOCHA object

- [`addMotifSet()`](https://aifimmunology.github.io/MOCHA/reference/addMotifSet.md)
  : Identify motifs within a peakset

- [`annotateTiles()`](https://aifimmunology.github.io/MOCHA/reference/annotateTiles.md)
  : Annotate tiles with gene annotations

- [`asMochaSTM()`](https://aifimmunology.github.io/MOCHA/reference/asMochaSTM.md)
  : Coerce to MochaSampleTileMatrix

- [`asMochaTileResults()`](https://aifimmunology.github.io/MOCHA/reference/asMochaTileResults.md)
  : Coerce to MochaTileResults

- [`cellTypes()`](https://aifimmunology.github.io/MOCHA/reference/cellTypes.md)
  : Get cell population names from a MOCHA object

- [`classifyZeros()`](https://aifimmunology.github.io/MOCHA/reference/classifyZeros.md)
  : Classify observed zeros as technical, biological, or ambiguous

- [`combineSampleTileMatrix()`](https://aifimmunology.github.io/MOCHA/reference/combineSampleTileMatrix.md)
  : Merge the TSAM from multiple cell populations into a single matrix

- [`correctGenome()`](https://aifimmunology.github.io/MOCHA/reference/correctGenome.md)
  :

  `correctGenome`

- [`differentialsToGRanges()`](https://aifimmunology.github.io/MOCHA/reference/differentialsToGRanges.md)
  : Convert a data.frame or matrix to a GRanges

- [`extractRegion()`](https://aifimmunology.github.io/MOCHA/reference/extractRegion.md)
  : Extract accessibility coverage for a given region

- [`getCellPopMatrix()`](https://aifimmunology.github.io/MOCHA/reference/getCellPopMatrix.md)
  : Get the SampleTileMatrix of the given cell population

- [`getCellTypeTiles()`](https://aifimmunology.github.io/MOCHA/reference/getCellTypeTiles.md)
  : Extract the GRanges for a particular cell population

- [`getCellTypes()`](https://aifimmunology.github.io/MOCHA/reference/getCellTypes.md)
  : Extract cell population names from a Tile Results or Sample Tile
  object.

- [`getCoverage()`](https://aifimmunology.github.io/MOCHA/reference/getCoverage.md)
  : Get sample-specific coverage files for each sample-cell population

- [`getDropoutProb()`](https://aifimmunology.github.io/MOCHA/reference/getDropoutProb.md)
  : Get predicted P(zero) dropout scores for a cell population

- [`getModelValues()`](https://aifimmunology.github.io/MOCHA/reference/getModelValues.md)
  : Get a data.frame of model values from the output of linear modeling

- [`getOpenTiles()`](https://aifimmunology.github.io/MOCHA/reference/getOpenTiles.md)
  : Get per-cell-population open tiles from a MOCHA tileResults object

- [`getPopFrags()`](https://aifimmunology.github.io/MOCHA/reference/getPopFrags.md)
  : Extract fragments by populations from an ArchR Project

- [`getPromoterGenes()`](https://aifimmunology.github.io/MOCHA/reference/getPromoterGenes.md)
  :

  Extract the list of promoter genes from a GRanges annotated with
  [`annotateTiles()`](https://aifimmunology.github.io/MOCHA/reference/annotateTiles.md)

- [`getSampleCellTypeMetadata()`](https://aifimmunology.github.io/MOCHA/reference/getSampleCellTypeMetadata.md)
  : Extract Sample-celltype specific metadata

- [`isMOCHAObject()`](https://aifimmunology.github.io/MOCHA/reference/isMOCHAObject.md)
  : isMOCHAObject

- [`mergeTileResults()`](https://aifimmunology.github.io/MOCHA/reference/mergeTileResults.md)
  :

  Merge tileResults from
  [`callOpenTiles()`](https://aifimmunology.github.io/MOCHA/reference/callOpenTiles-methods.md)
  that share cell populations into a single object containing all
  samples

- [`openTiles()`](https://aifimmunology.github.io/MOCHA/reference/openTiles.md)
  : Get open tiles from a MochaTileResults object

- [`renameCellTypes()`](https://aifimmunology.github.io/MOCHA/reference/renameCellTypes.md)
  :

  Modify the cell population names in a Sample-Tile Object from
  [`getSampleTileMatrix()`](https://aifimmunology.github.io/MOCHA/reference/getSampleTileMatrix.md)

- [`seuratToMOCHAInputs()`](https://aifimmunology.github.io/MOCHA/reference/seuratToMOCHAInputs.md)
  : Convert Seurat/Signac object to sample-level MOCHA fragment inputs.

- [`subsetMOCHAObject()`](https://aifimmunology.github.io/MOCHA/reference/subsetMOCHAObject.md)
  : Subset a tileResults object by metadata

- [`updateDirectoryPath()`](https://aifimmunology.github.io/MOCHA/reference/updateDirectoryPath.md)
  : updateDirectoryPath
