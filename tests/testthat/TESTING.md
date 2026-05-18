# MOCHA test profiles

## Default (CI / `R CMD check`)

Runs deterministic unit tests using bundled internal data (`testTileResults`,
`testTileResultsMultisample`, `exampleFragments`, fixtures under `fixtures/`).

```bash
R -e 'testthat::test_local()'
```

## Heavy integration

Requires local ArchR projects and optional coverage directories. Enable with:

```bash
export MOCHA_HEAVY_TESTS=true
```

Optional paths (searched relative to package root and parent directories):

| Resource | Typical location |
|----------|------------------|
| `PBMCSmall` | `../PBMCSmall` |
| `FullCovid` | `../FullCovid` |
| `HemeTutorial/MOCHA` | `../HemeTutorial/MOCHA` |

Heavy-only tests: `test_COVID_data_pipeline.R`, `test_plotRegion.R`, ArchR
branches in `test_callOpenTiles.R`, `test_packMOCHA.R`, `test_getPopFrags.R`.

Optional MSEA regression (when reference CSVs exist):

```bash
export MOCHA_MSEA_INPUT=/path/to/input_motifenrichment.csv
export MOCHA_MSEA_EXPECTED=/path/to/results_MSEA.csv
```

## Coverage guard

`test_exported_api_guard.R` fails if a newly exported function is not referenced
in any test file and is not listed in the allowlist inside that test.

See [../COVERAGE_MATRIX.md](../COVERAGE_MATRIX.md) for the full API matrix.
