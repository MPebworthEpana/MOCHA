#!/usr/bin/env bash
# Run MOCHA tests in default and heavy profiles (requires R with package deps).
# Recommended: conda env mocha-test (see environment-mocha-test.yml at repo root).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

if command -v conda >/dev/null 2>&1; then
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate mocha-test 2>/dev/null || true
fi

export RETICULATE_PYTHON="${RETICULATE_PYTHON:-${CONDA_PREFIX}/bin/python}"

echo "== Default CI profile =="
echo "R: $(which Rscript)"
echo "Python (reticulate): ${RETICULATE_PYTHON}"
Rscript -e 'devtools::load_all(); testthat::test_local(filter = NULL)'

echo "== Heavy integration profile =="
MOCHA_HEAVY_TESTS=true Rscript -e 'devtools::load_all(); testthat::test_local(filter = NULL)'

echo "Done."
