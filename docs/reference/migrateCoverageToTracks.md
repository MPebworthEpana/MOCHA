# Migrate legacy coverage RDS bundles to structured bigWig tracks

Converts \`cellPop_CoverageFiles.RDS\` bundles written by older MOCHA
versions into the structured \`tracks/\` layout. Legacy files are kept
unless \`removeLegacy = TRUE\`.

## Usage

``` r
migrateCoverageToTracks(
  MOCHAObj,
  cellPopulations = "ALL",
  removeLegacy = FALSE,
  force = FALSE,
  verbose = FALSE
)
```

## Arguments

- MOCHAObj:

  A MOCHA tileResults or Sample-Tile Matrix object with
  \`metadata\$Directory\` set.

- cellPopulations:

  Character vector of cell populations to migrate, or \`"ALL"\` for all
  assays.

- removeLegacy:

  If TRUE, delete legacy RDS bundles after successful migration. Default
  is FALSE.

- force:

  If TRUE, overwrite existing structured track files. Default is FALSE.

- verbose:

  Display additional messages. Default is FALSE.

## Value

The input MOCHA object with \`metadata\$CoverageLayout\` set to
\`"tracks"\`.
