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

### Full suite (serial, recommended locally)

Parallel test workers can clash on `makeCluster` ports. Run serially:

```bash
TESTTHAT_PARALLEL=false Rscript tests/scripts/run_test_local.R
```

Equivalent one-liner:

```bash
TESTTHAT_PARALLEL=false Rscript -e 'testthat::set_max_fails(Inf); pkgload::load_all(); testthat::test_local()'
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
branches in `test_callOpenTiles.R`, `test_packMOCHA.R`, `test_getPopFrags.R`,
`test_extractRegion.R` (HemeTutorial coverage branch).

### MSEA heavy regression

Bundled defaults (no setup required when heavy tests are enabled):

- `tests/testthat/fixtures/msea_motif_enrichment.csv`
- `tests/testthat/fixtures/msea_expected_results.csv`
- `tests/testthat/fixtures/ligand_tf_matrix_mini.rds` (offline ligand–TF matrix; no network)

```bash
export MOCHA_HEAVY_TESTS=true
# optional overrides:
export MOCHA_MSEA_INPUT=/path/to/input_motifenrichment.csv
export MOCHA_MSEA_EXPECTED=/path/to/results_MSEA.csv
```

If `ligand_tf_matrix_mini.rds` is absent, the heavy test falls back to downloading
`ligand_tf_matrix.rds` from Zenodo (requires network and `curl`).

### Synthetic multisample fixtures (no ArchR)

These tests use `make_synthetic_sample_tile_matrix()` and do not need local projects:

- `test_filterCoAccessibleLinks.R` (co-accessibility correlations across samples; optional `ZI = TRUE` path)
- `test_modeling_ziglmm.R` (`linearModeling` with `exp ~ GroupA + (1 | GroupA)`)

### Coverage for `extractRegion`

**Default** `test_extractRegion.R` always uses `make_synthetic_coverage_dir()` (no
`HemeTutorial` dependency). Real ArchR/MOCHA output at `../HemeTutorial/MOCHA` is
exercised only when `MOCHA_HEAVY_TESTS=true` (same layout as `callOpenTiles` exports:
`C3_CoverageFiles.RDS`, etc.).

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

# End-to-end smoke: seuratToMOCHAInputs + callOpenTiles(Seurat)
# Uses MOCHA::exampleBlackList when available, else chr1:1-750100 mini blacklist
NOT_CRAN=true Rscript tests/scripts/smoke_seurat_ingest.R
```

Indexed fragment fixture: `tests/testthat/fixtures/seurat/mini_fragments.tsv.gz` (+ `.tbi`).

Smoke helpers: `tests/scripts/smoke_seurat_helpers.R` (`smoke_load_blacklist`, dependency checks).
Requires `BSgenome.Hsapiens.UCSC.hg19`, `TxDb.Hsapiens.UCSC.hg38.knownGene`, and `org.Hs.eg.db`.
The synthetic object has &lt;5 cells per sample; MOCHA may warn that samples are ignored while still returning a `MultiAssayExperiment`.

## Coverage guard

`test_exported_api_guard.R` fails if a newly exported function is not referenced
in any test file and is not listed in the allowlist inside that test.

See [../COVERAGE_MATRIX.md](../COVERAGE_MATRIX.md) for the full API matrix.

## Tutorial scripts (`inst/tutorials/`)

Vignette **code** is maintained in `inst/tutorials/*.R` and included in
`vignettes/*.Rmd` via `knitr::read_chunk()`. Runnable scripts are verified in CI
(`.github/workflows/tutorial-check.yml`).

| Tier | Scripts | How to run |
|------|---------|------------|
| A | `01-workflow.R`, `02-import-bundled.R` | `Rscript tests/scripts/run_tutorials.R --tier a` |
| B | `00-fixtures.R`, `03`–`06` | `NOT_CRAN=true Rscript tests/scripts/run_tutorials.R --tier b` |
| C | `reference/archr-call-open-tiles.R`, `reference/signac-smoke.R` | `MOCHA_HEAVY_TESTS=true` + local ArchR/Signac fixtures |

Validate chunk label sync:

```bash
Rscript scripts/validate_tutorial_chunks.R
```

Testthat entry points: `test-tutorial-tier-a.R`, `test-tutorial-tier-b.R`,
`test-tutorial-tier-c.R`.
