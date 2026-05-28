#!/usr/bin/env bash
# Build MOCHA pkgdown articles on WSL/Linux using the mocha-docs conda environment.
#
# Usage (from repo root):
#   ./scripts/build_pkgdown_articles.sh
#   ./scripts/build_pkgdown_articles.sh --full-site
#   ./scripts/build_pkgdown_articles.sh --no-clean --full-site
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

# Default to mocha-test when present (full MOCHA deps); override with MOCHA_DOCS_ENV.
ENV_NAME="${MOCHA_DOCS_ENV:-mocha-docs}"
if [[ -z "${MOCHA_DOCS_ENV:-}" ]] && conda env list | awk '{print $1}' | grep -qx "mocha-test"; then
  ENV_NAME="mocha-test"
  MOCHA_SKIP_ENV_UPDATE="${MOCHA_SKIP_ENV_UPDATE:-1}"
fi
ENV_FILE="${REPO_ROOT}/environment-mocha-docs.yml"

find_conda() {
  if [[ -n "${CONDA_EXE:-}" && -x "${CONDA_EXE}" ]]; then
    echo "${CONDA_EXE}"
    return 0
  fi
  local candidate
  for candidate in \
    "${HOME}/miniforge3/bin/conda" \
    "${HOME}/mambaforge/bin/conda" \
    "${HOME}/miniconda3/bin/conda" \
    "${HOME}/anaconda3/bin/conda" \
    "/opt/conda/bin/conda"
  do
    if [[ -x "${candidate}" ]]; then
      echo "${candidate}"
      return 0
    fi
  done
  return 1
}

CONDA_BIN="$(find_conda)" || {
  echo "ERROR: conda not found. Install Miniforge/Mambaforge or set CONDA_EXE." >&2
  exit 1
}

# shellcheck source=/dev/null
source "$("$CONDA_BIN" info --base)/etc/profile.d/conda.sh"

if [[ "${MOCHA_SKIP_ENV_UPDATE:-}" == "1" ]]; then
  echo "Skipping conda env create/update (MOCHA_SKIP_ENV_UPDATE=1); using existing '${ENV_NAME}'."
elif ! conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  echo "Creating conda environment '${ENV_NAME}' from ${ENV_FILE} ..."
  conda env create -f "${ENV_FILE}"
else
  echo "Updating conda environment '${ENV_NAME}' ..."
  conda env update -f "${ENV_FILE}" --prune
fi

conda activate "${ENV_NAME}"

# Ensure vignette-only Bioconductor deps if conda solve omitted them
Rscript -e '
  pkgs <- c("BiocStyle", "TxDb.Hsapiens.UCSC.hg38.knownGene", "BSgenome.Hsapiens.UCSC.hg19")
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) stop("BiocManager required")
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  }
'

echo "Using R: $(which R)"
R --version | head -1

export R_LIBS_USER="${CONDA_PREFIX}/lib/R/library"
export R_LIBS="${CONDA_PREFIX}/lib/R/library"

# Log outside docs/ — full-site builds delete docs/ before pkgdown runs.
Rscript "${REPO_ROOT}/scripts/build_pkgdown_articles.R" "$@"
