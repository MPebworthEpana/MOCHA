# MOCHA runnable tutorial scripts

Tutorial **code** lives here as the single source of truth. Vignettes under
[`vignettes/`](../../vignettes/) include these files with `knitr::read_chunk()`;
pkgdown builds HTML from the vignettes.

## Run locally

From the MOCHA repository root with MOCHA installed or `devtools::load_all()`:

```bash
# All tiers (A + B; C when MOCHA_HEAVY_TESTS=true)
NOT_CRAN=true Rscript tests/scripts/run_tutorials.R

# Tier A only (Bioconductor-safe bundled paths)
Rscript tests/scripts/run_tutorials.R --tier a

# Validate chunk labels match vignettes
Rscript scripts/validate_tutorial_chunks.R
```

## Tiers

| Tier | Scripts | CI |
|------|---------|-----|
| **A** | `01-workflow.R`, `02-import-bundled.R` | Bioconductor Docker / every PR |
| **B** | `03-downstream.R` … `06-alt-tss-motifs.R` | `NOT_CRAN=true` tutorial-check workflow |
| **C** | `reference/*.R` | `MOCHA_HEAVY_TESTS=true` (ArchR, Signac smoke, etc.) |

## Files

| Script | Vignette |
|--------|----------|
| `01-workflow.R` | `MOCHA-workflow-tutorial.Rmd` |
| `02-import-bundled.R` | `Data-Import-Tutorial.Rmd` (bundled section) |
| `03-downstream.R` | `MOCHA-downstream-workflows.Rmd` |
| `04-export.R` | `MOCHA-export-and-sharing.Rmd` |
| `05-advanced-modeling.R` | `MOCHA-advanced-modeling.Rmd` |
| `06-alt-tss-motifs.R` | `Alternative-TSS-TF-regulation.Rmd` |
| `00-fixtures.R` | Shared objects for Tier B scripts |
| `_shared.R` | Helpers (not knitr chunks) |

Reference import paths (ArchR, Signac downloads, SnapATAC) remain in the Data
Import vignette as `eval = FALSE` chunks only.

## Chunk layout

Each script must start with `# ---- libraries ----` (no executable code before the
first label). Tier B scripts use a `# ---- script-init ----` chunk for fixture
loading; that chunk is runner-only (not referenced in vignettes) so Bioconductor
builds only execute `library(MOCHA)` from the `libraries` chunk.

Pkgdown scans `inst/tutorials/` and requires the `rsconnect` package (listed in
`environment-mocha-docs.yml`).
