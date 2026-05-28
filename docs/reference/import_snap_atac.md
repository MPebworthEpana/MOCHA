# Import SnapATAC2 AnnData / AnnDataSet for MOCHA peak calling

Converts a SnapATAC2 \`AnnData\` or \`AnnDataSet\` object (or path to
on-disk \`.h5ad\`) into MOCHA-native inputs for \[callOpenTiles()\].
Requires optional Python packages via \*\*reticulate\*\* (\`anndata\`
preferred; \`snapatac2\` for \`AnnDataSet\` objects).

## Usage

``` r
import_snap_atac(
  x,
  cellPopLabel,
  sampleColumn = NULL,
  cellCol = "RG",
  fragmentsLayer = c("auto", "fragment_paired", "fragment_single"),
  chromPrefix = c("auto", "add", "strip", "none"),
  pythonModule = c("auto", "anndata", "snapatac2"),
  verbose = FALSE
)
```

## Arguments

- x:

  A reticulate Python object (\`AnnData\` / \`AnnDataSet\`), a character
  path to an \`.h5ad\` file, or a named list of paths (names used as
  sample IDs when \`sampleColumn\` is absent from \`obs\`).

- cellPopLabel:

  Column in \`obs\` containing cell population / cluster labels.

- sampleColumn:

  Column in \`obs\` with biological sample IDs. Default \`NULL\`: uses
  \`"sample"\` for \`AnnDataSet\`, otherwise must be supplied for
  \`AnnData\`.

- cellCol:

  Name of the metadata column stamped onto each fragment GRanges
  (default \`"RG"\`, matching MOCHA / ArchR conventions).

- fragmentsLayer:

  Which \`.obsm\` layer to use: \`"auto"\`, \`"fragment_paired"\`, or
  \`"fragment_single"\`.

- chromPrefix:

  Chromosome name harmonization: \`"auto"\` (no change), \`"add"\`
  (prepend \`chr\` if missing), \`"strip"\`, or \`"none"\`.

- pythonModule:

  Python module for I/O: \`"auto"\`, \`"anndata"\`, or \`"snapatac2"\`.

- verbose:

  Logical; print progress messages.

## Value

A named list with elements suitable for splatting into
\[callOpenTiles()\]:

- ATACFragments:

  Sample-level `GRangesList`.

- cellColData:

  `data.frame` with rownames = cell barcodes.

- cellCol:

  Character; column name on fragment metadata.

- seqinfo:

  `Seqinfo` from SnapATAC2 `reference_sequences`.

You must still provide `genome` and `blackList` to \[callOpenTiles()\].

## Details

SnapATAC2 stores fragments in \`.obsm\` as a CSR-like sparse matrix
(rows = cells, columns = global genomic start positions, values =
fragment length). Decode offsets use \`uns\['reference_sequences'\]\`.
See SnapATAC2 documentation: <https://scverse.org/SnapATAC2/>.
