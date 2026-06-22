#!/usr/bin/env Rscript
# Export GLMM results (RDS) to stats CSVs for raincloud brackets and full model tables.
# Run after AOC_stats_glmm_nback.R and AOC_stats_glmm_sternberg.R.

options(scipen = 999)

script_dir <- {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("--file=", "", args[grep("^--file=", args)])
  if (length(file_arg) > 0) {
    dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE))
  } else {
    getwd()
  }
}

source(file.path(script_dir, "AOC_stats_glmm_export_helpers.R"))

message("=== AOC GLMM export ===")
export_all_tasks()

message("\n=== Full model Word tables ===")
source(file.path(script_dir, "AOC_stats_fullModelTables.R"), local = TRUE)
