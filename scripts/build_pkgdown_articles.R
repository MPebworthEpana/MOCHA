#!/usr/bin/env Rscript
# Build pkgdown articles (vignettes -> docs/articles/) for MOCHA.
#
# Usage (from repo root, with mocha-docs conda env active):
#   Rscript scripts/build_pkgdown_articles.R
#   Rscript scripts/build_pkgdown_articles.R --full-site
#   Rscript scripts/build_pkgdown_articles.R --no-install
#
args <- commandArgs(trailingOnly = TRUE)
full_site <- "--full-site" %in% args
no_install <- "--no-install" %in% args
clean <- !("--no-clean" %in% args)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg)) {
  script_path <- sub("^--file=", "", file_arg[1])
  repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
} else {
  repo_root <- normalizePath(getwd(), mustWork = TRUE)
}
if (!file.exists(file.path(repo_root, "DESCRIPTION"))) {
  stop("Run from the MOCHA repository root (or via scripts/build_pkgdown_articles.sh).")
}
setwd(repo_root)
message("Repository root: ", repo_root)

if (!no_install) {
  message("Installing MOCHA from source ...")
  status <- system2("R", c("CMD", "INSTALL", ".", "--no-multiarch", "--with-keep.source"))
  if (status != 0) {
    stop("R CMD INSTALL failed (exit ", status, ").")
  }
}

if (!requireNamespace("pkgdown", quietly = TRUE)) {
  stop(
    "Package 'pkgdown' is required. Activate the mocha-docs conda env:\n",
    "  conda activate mocha-docs"
  )
}

if (full_site) {
  message("Building full pkgdown site (reference + articles + ...) ...")
  if (clean && dir.exists("docs")) {
    message("Removing existing docs/ (--no-clean to skip).")
    unlink("docs", recursive = TRUE)
  }
  if (!dir.exists("docs")) {
    dir.create("docs", showWarnings = FALSE)
  }
  pkgdown::build_site_github_pages(
    pkg = ".",
    clean = clean,
    install = FALSE
  )
} else {
  message("Building articles only (docs/articles/) ...")
  if (!dir.exists("docs")) {
    dir.create("docs", showWarnings = FALSE)
  }
  pkgdown::build_articles()
  if (file.exists(file.path("docs", "reference"))) {
    message("Refreshing search index ...")
    pkgdown::build_search(pkg = ".")
  }
}

built <- list.files(
  file.path("docs", "articles"),
  pattern = "\\.html$",
  full.names = FALSE
)
message(
  "Done. Article HTML files: ",
  if (length(built)) paste(built, collapse = ", ") else "(none — check errors above)"
)
