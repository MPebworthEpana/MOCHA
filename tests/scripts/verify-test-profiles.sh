#!/usr/bin/env bash
# Run MOCHA tests in default and heavy profiles (requires R with package deps).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

echo "== Default CI profile =="
Rscript -e 'testthat::test_local(filter = NULL)'

echo "== Heavy integration profile =="
MOCHA_HEAVY_TESTS=true Rscript -e 'testthat::test_local(filter = NULL)'

echo "Done."
