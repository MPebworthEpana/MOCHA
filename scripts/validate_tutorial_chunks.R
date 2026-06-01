#!/usr/bin/env Rscript
# Validate that vignette chunk labels exist in inst/tutorials/*.R scripts.

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg)) {
  repo_root <- normalizePath(
    file.path(dirname(sub("^--file=", "", file_arg[1])), ".."),
    mustWork = TRUE
  )
} else {
  repo_root <- normalizePath(getwd(), mustWork = TRUE)
}

# Map vignette -> script; optional explicit chunk allowlist (reference chunks stay in Rmd).
vignette_map <- list(
  "MOCHA-workflow-tutorial.Rmd" = list(
    script = "01-workflow.R",
    chunks = NULL
  ),
  "Data-Import-Tutorial.Rmd" = list(
    script = "02-import-bundled.R",
    chunks = c(
      "check-deps",
      "bundled-libraries",
      "bundled-call-open-tiles",
      "bundled-get-open-tiles",
      "import-utilities",
      "session-info"
    )
  ),
  "MOCHA-downstream-workflows.Rmd" = list(
    script = "03-downstream.R",
    chunks = NULL
  ),
  "MOCHA-export-and-sharing.Rmd" = list(
    script = "04-export.R",
    chunks = NULL
  ),
  "MOCHA-advanced-modeling.Rmd" = list(
    script = "05-advanced-modeling.R",
    chunks = NULL
  ),
  "Alternative-TSS-TF-regulation.Rmd" = list(
    script = "06-alt-tss-motifs.R",
    chunks = NULL
  )
)

extract_chunk_labels <- function(rmd_path) {
  lines <- readLines(rmd_path, warn = FALSE)
  hits <- grep("^```\\{r ([^,}]+)", lines, value = TRUE)
  unique(sub("^```\\{r ([^,}]+).*$", "\\1", hits))
}

extract_script_labels <- function(r_path) {
  lines <- readLines(r_path, warn = FALSE)
  hits <- grep("^# ---- ([^ ]+) ----", lines, value = TRUE)
  sub("^# ---- ([^ ]+) ----.*", "\\1", hits)
}

errors <- character()
for (vig in names(vignette_map)) {
  entry <- vignette_map[[vig]]
  rmd <- file.path(repo_root, "vignettes", vig)
  script <- file.path(repo_root, "inst", "tutorials", entry$script)
  if (!file.exists(rmd) || !file.exists(script)) {
    errors <- c(errors, paste("Missing file for", vig))
    next
  }
  chunk_labels <- extract_chunk_labels(rmd)
  chunk_labels <- chunk_labels[chunk_labels != "setup"]
  if (!is.null(entry$chunks)) {
    chunk_labels <- entry$chunks
  }
  script_labels <- extract_script_labels(script)
  missing <- setdiff(chunk_labels, script_labels)
  if (length(missing)) {
    errors <- c(
      errors,
      paste0(vig, ": missing labels in ", entry$script, ": ",
             paste(missing, collapse = ", "))
    )
  }
}

if (length(errors)) {
  cat(paste(errors, collapse = "\n"), "\n")
  quit(status = 1L)
}

cat("All vignette chunk labels found in tutorial scripts.\n")
quit(status = 0L)
