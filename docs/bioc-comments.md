# Bioconductor submission check notes (MOCHA 2.0.0)

Branch: `bioc/submission`

## Local verification

Run on R-devel with Bioconductor dependencies installed:

``` bash
R CMD build MOCHA
R CMD check --no-manual MOCHA_2.0.0.tar.gz
R -q -e 'BiocCheck::BiocCheck("MOCHA_2.0.0.tar.gz")'
```

Document results below after running locally or from the `Bioc-check`
GitHub Action.

## Expected check status (Phase 1)

- **R CMD check:** target 0 ERROR, 0 WARNING on Linux devel.
- **BiocCheck:** ERRORs must be zero; NOTEs from style (function length,
  `T`/`F`, etc.) are deferred to review iteration.

## Known acceptable items (Phase 1)

- `ArchR` remains in `Suggests`; not available on Bioconductor. Examples
  and vignette use bundled `GRangesList` data only.
- Six vignettes: `MOCHA-workflow-tutorial.Rmd` (runnable end-to-end on
  bundled data), `Data-Import-Tutorial.Rmd` (runnable bundled section
  only; other import paths `eval = FALSE`),
  `MOCHA-downstream-workflows.Rmd`, `MOCHA-export-and-sharing.Rmd`, and
  `MOCHA-advanced-modeling.Rmd` (reference-only, `eval = FALSE`),
  `Alternative-TSS-TF-regulation.Rmd` (mostly `eval = FALSE` templates).
- Only chunks using bundled `exampleFragments` / package data run on
  Bioconductor builders; all other chunks use `eval = FALSE` or
  `NOT_CRAN` guards.
- Tutorial code lives in `inst/tutorials/*.R` (single source of truth);
  vignettes use
  [`knitr::read_chunk()`](https://rdrr.io/pkg/knitr/man/read_chunk.html).
  Tier A/B scripts are exercised by `tutorial-check` CI and
  `tests/scripts/run_tutorials.R`.
- `ImportingFromOtherSources.Rmd` merged into `Data-Import-Tutorial.Rmd`
  (removed).
- Pre-built vignette HTML under `vignettes/` is excluded via
  `.Rbuildignore`.

## Maintainer actions before opening Contributions PR

1.  Subscribe maintainer to
    [bioc-devel](https://stat.ethz.ch/mailman/listinfo/bioc-devel).
2.  Open an issue at <https://github.com/Bioconductor/Contributions>
    using the New Package template.
3.  Point the repository to `https://github.com/aifimmunology/MOCHA`
    branch `bioc/submission`.

## Check results

Local R is not available in the development environment used for this
branch. Run verification via the `Bioc-check` GitHub Action on push to
`bioc/submission`, or locally:

``` bash
docker run --rm -v "$PWD":/workspace -w /workspace \
  bioconductor/bioconductor_docker:devel bash -c '
  apt-get update -qq && apt-get install -y -qq libcurl4-openssl-dev libssl-dev libxml2-dev
  R -q -e "install.packages(\"remotes\"); remotes::install_deps(suggests=TRUE)"
  R CMD build . && R CMD check --no-manual MOCHA_2.0.0.tar.gz
  R -q -e "BiocCheck::BiocCheck(\"MOCHA_2.0.0.tar.gz\")"
'
```

    R CMD check: (pending CI / local docker run)
      ERRORs: pending
      WARNINGs: pending
      NOTEs: pending

    BiocCheck: (pending CI / local docker run)
      ERRORS: pending
      WARNINGS: pending
      NOTES: pending
