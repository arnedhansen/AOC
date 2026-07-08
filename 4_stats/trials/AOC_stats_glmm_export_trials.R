#!/usr/bin/env Rscript
# Export trial-level GLMM results (RDS) to stats CSVs for full model tables.
# Run after AOC_stats_glmm_nback_trials.R and AOC_stats_glmm_sternberg_trials.R.
# Outputs: data/stats/trials (separate from condition-level data/stats).

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

source(file.path(script_dir, "AOC_stats_glmm_export_helpers_trials.R"))

message("=== AOC GLMM export (trials) ===")
export_all_tasks_trials()

message("\n=== Full model Word tables (trials) ===")
source(file.path(script_dir, "AOC_stats_fullModelTables_trials.R"), local = TRUE)
