#!/usr/bin/env Rscript
# Run MOCHA tutorial scripts (inst/tutorials/) by verification tier.
#
# Usage (from repo root):
#   Rscript tests/scripts/run_tutorials.R
#   Rscript tests/scripts/run_tutorials.R --tier a
#   NOT_CRAN=true Rscript tests/scripts/run_tutorials.R --tier b
#   MOCHA_HEAVY_TESTS=true Rscript tests/scripts/run_tutorials.R --tier c

args <- commandArgs(trailingOnly = TRUE)
tier <- "all"
if ("--tier" %in% args) {
  idx <- match("--tier", args)
  tier <- tolower(args[idx + 1L])
}

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg)) {
  script_path <- sub("^--file=", "", file_arg[1])
  repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
} else {
  repo_root <- normalizePath(getwd(), mustWork = TRUE)
}
setwd(repo_root)

Sys.setenv(MOCHA_TUTORIAL_DIR = file.path(repo_root, "inst", "tutorials"))

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(repo_root, quiet = TRUE)
} else if (!requireNamespace("MOCHA", quietly = TRUE)) {
  stop("Install devtools or MOCHA before running tutorials")
}

tutorials_path <- file.path(repo_root, "inst", "tutorials")
shared <- file.path(tutorials_path, "_shared.R")
if (file.exists(shared)) {
  source(shared, local = FALSE)
}

source_script <- function(file, env) {
  path <- file.path(tutorials_path, file)
  if (!file.exists(path)) {
    stop("Missing tutorial script: ", path)
  }
  message("==> ", file)
  source(path, local = env)
  invisible(TRUE)
}

run_tier_a <- function() {
  env <- new.env(parent = globalenv())
  source_script("01-workflow.R", env)
  source_script("02-import-bundled.R", new.env(parent = globalenv()))
  invisible(env)
}

run_tier_b <- function() {
  env <- new.env(parent = globalenv())
  source_script("00-fixtures.R", env)
  for (f in c(
    "03-downstream.R",
    "04-export.R",
    "05-advanced-modeling.R",
    "06-alt-tss-motifs.R"
  )) {
    child <- new.env(parent = env)
    source_script(f, child)
  }
  invisible(env)
}

run_tier_c <- function() {
  if (!identical(tolower(Sys.getenv("MOCHA_HEAVY_TESTS", "false")), "true") &&
      !tolower(Sys.getenv("MOCHA_HEAVY_TESTS", "false")) %in% c("1", "yes")) {
    message("Skipping Tier C (set MOCHA_HEAVY_TESTS=true)")
    return(invisible(NULL))
  }
  ref_dir <- file.path(tutorials_path, "reference")
  archr <- file.path(ref_dir, "archr-call-open-tiles.R")
  if (file.exists(archr)) {
    source_script("reference/archr-call-open-tiles.R", new.env(parent = globalenv()))
  }
  signac <- file.path(ref_dir, "signac-smoke.R")
  if (file.exists(signac)) {
    source(signac, local = FALSE)
  }
  invisible(NULL)
}

status <- 0L
tryCatch({
  if (tier %in% c("all", "a")) {
    run_tier_a()
  }
  if (tier %in% c("all", "b")) {
    if (!identical(Sys.getenv("NOT_CRAN", "false"), "true")) {
      message("Skipping Tier B (set NOT_CRAN=true)")
    } else {
      run_tier_b()
    }
  }
  if (tier %in% c("all", "c")) {
    run_tier_c()
  }
  message("Tutorial scripts completed.")
}, error = function(e) {
  message("Tutorial run failed: ", conditionMessage(e))
  status <<- 1L
})

quit(status = status)
