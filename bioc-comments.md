# Bioconductor submission check notes (MOCHA 0.99.0)

Branch: `bioc/submission`

## Local verification

Run on R-devel with Bioconductor dependencies installed:

```bash
R CMD build MOCHA
R CMD check --no-manual MOCHA_0.99.0.tar.gz
R -q -e 'BiocCheck::BiocCheck("MOCHA_0.99.0.tar.gz")'
```

Document results below after running locally or from the `Bioc-check` GitHub Action.

## Expected check status (Phase 1)

- **R CMD check:** target 0 ERROR, 0 WARNING on Linux devel.
- **BiocCheck:** ERRORs must be zero; NOTEs from style (function length, `T`/`F`, etc.) are deferred to review iteration.

## Known acceptable items (Phase 1)

- `ArchR` remains in `Suggests`; not available on Bioconductor. Examples and vignette use bundled `GRangesList` data only.
- `ImportingFromOtherSources.Rmd` uses `eval = FALSE` (reference tutorial); primary runnable vignette is `COVID-walkthrough.Rmd`.
- Pre-built vignette HTML under `vignettes/` is excluded via `.Rbuildignore`.

## Maintainer actions before opening Contributions PR

1. Subscribe maintainer to [bioc-devel](https://stat.ethz.ch/mailman/listinfo/bioc-devel).
2. Open an issue at https://github.com/Bioconductor/Contributions using the New Package template.
3. Point the repository to `https://github.com/aifimmunology/MOCHA` branch `bioc/submission`.

## Check results

_Fill in after local or CI run:_

```
R CMD check:
  ERRORs: 
  WARNINGs: 
  NOTEs: 

BiocCheck:
  ERRORS: 
  WARNINGS: 
  NOTES: 
```
