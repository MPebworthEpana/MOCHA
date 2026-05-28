#!/usr/bin/env Rscript
# Validate that vignette chunk labels exist in inst/tutorials/*.R scripts.

repo_root <- normalizePath(
  getwd(),
  mustWork = FALSE
)
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg)) {
  repo_root <- normalizePath(
    file.path(dirname(sub("^--file=", "", file_arg[1])), ".."),
    mustWork = TRUE
  )
}

vignette_map <- c(
  "MOCHA-workflow-tutorial.Rmd" = "01-workflow.R",
  "Data-Import-Tutorial.Rmd" = "02-import-bundled.R",
  "MOCHA-downstream-workflows.Rmd" = "03-downstream.R",
  "MOCHA-export-and-sharing.Rmd" = "04-export.R",
  "MOCHA-advanced-modeling.Rmd" = "05-advanced-modeling.R",
  "Alternative-TSS-TF-regulation.Rmd" = "06-alt-tss-motifs.R"
)

extract_chunk_labels <- function(rmd_path) {
  lines <- readLines(rmd_path, warn = FALSE)
  grep("^```\\{r ([^,}]+)", lines, value = TRUE)
}

extract_script_labels <- function(r_path) {
  lines <- readLines(r_path, warn = FALSE)
  hits <- grep("^# ---- ([^ ]+) ----", lines, value = TRUE)
  sub("^# ---- ([^ ]+) ----.*", "\\1", hits)
}

errors <- character()
for (vig in names(vignette_map)) {
  rmd <- file.path(repo_root, "vignettes", vig)
  script <- file.path(repo_root, "inst", "tutorials", vignette_map[[vig]])
  if (!file.exists(rmd) || !file.exists(script)) {
    errors <- c(errors, paste("Missing file for", vig))
    next
  }
  chunk_lines <- extract_chunk_labels(rmd)
  chunk_labels <- unique(sub("^```\\{r ([^,}]+).*$", "\\1", chunk_lines))
  chunk_labels <- chunk_labels[chunk_labels != "setup"]
  script_labels <- extract_script_labels(script)
  missing <- setdiff(chunk_labels, script_labels)
  if (length(missing)) {
    errors <- c(
      errors,
      paste0(vig, ": missing labels in ", vignette_map[[vig]], ": ",
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
