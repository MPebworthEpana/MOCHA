# MOCHA tests

- **Unit / CI:** `tests/testthat/` (see [testthat/TESTING.md](testthat/TESTING.md))
- **Coverage matrix:** [COVERAGE_MATRIX.md](COVERAGE_MATRIX.md)

```bash
# Default
R -e 'devtools::test()'

# Heavy integration (local ArchR projects)
MOCHA_HEAVY_TESTS=true R -e 'devtools::test()'
```
