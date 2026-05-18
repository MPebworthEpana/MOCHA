# MOCHA test profiles

## Default (CI / `R CMD check`)

Runs deterministic unit tests using bundled internal data (`testTileResults`,
`testTileResultsMultisample`, `exampleFragments`, fixtures under `fixtures/`).

`tests/testthat/setup.R` sets `NOT_CRAN=true` when it is unset, so
`skip_on_cran()` tests run under local `devtools::test()` / `test_local()` without
exporting `NOT_CRAN`. To reproduce CRAN-skip behavior locally:
`NOT_CRAN=false R -e 'testthat::test_local()'`.

The `mocha-test` conda env ([`environment-mocha-test.yml`](../../environment-mocha-test.yml))
includes `r-cutpointr` for `trainPeakModel(..., threshMethod = "youden")` tests.
Optional CRAN suggests can also be installed via
[`tests/scripts/install-test-deps.R`](../scripts/install-test-deps.R).

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

## WSL / Seurat–Signac ingest

Use the **`mocha-test`** conda environment (`/home/enki/miniforge3/envs/mocha-test`):

```bash
source /home/enki/miniforge3/etc/profile.d/conda.sh
conda activate mocha-test
cd /path/to/MOCHA
mamba env update -n mocha-test -f environment-mocha-test.yml
Rscript tests/scripts/install-seurat-test-deps.R   # Signac (not on conda)
```

`tests/testthat/setup.R` prepends `$CONDA_PREFIX/bin` to `PATH` so `tabix`/`bgzip` from `htslib` are found.

```bash
# Documentation (preload ggbio)
Rscript -e 'library(ggbio); library(ensembldb); devtools::document()'

# Seurat ingest tests
NOT_CRAN=true Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test_seurat_fragments.R")'

# End-to-end smoke (synthetic Seurat or MOCHA_SEURAT_RDS=/path/to/obj.rds)
NOT_CRAN=true Rscript tests/scripts/smoke_seurat_ingest.R
```

Indexed fragment fixture: `tests/testthat/fixtures/seurat/mini_fragments.tsv.gz` (+ `.tbi`).

## Coverage guard

`test_exported_api_guard.R` fails if a newly exported function is not referenced
in any test file and is not listed in the allowlist inside that test.

See [../COVERAGE_MATRIX.md](../COVERAGE_MATRIX.md) for the full API matrix.
