#!/usr/bin/env Rscript
# AOC GLMMs — Sternberg
# Loads merged CSV and fits preregistered mixed models (lme4).
# Export tables: run AOC_stats_glmm_export.R afterward.

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

source(file.path(script_dir, "AOC_stats_glmm_fit.R"))

message("=== AOC GLMMs: Sternberg ===")
results <- run_paper_glmm(STERNBERG_TASK_CONFIG)
save_glmm_results(results)
message("Done. Run AOC_stats_glmm_export.R to write CSV / Word table inputs.")
