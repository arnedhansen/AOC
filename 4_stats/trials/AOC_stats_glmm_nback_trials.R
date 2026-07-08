#!/usr/bin/env Rscript
# AOC GLMMs — N-back (trial-level)
# Loads AOC_merged_data_nback_trials.csv and fits mixed models.
# Export tables: run AOC_stats_glmm_export_trials.R afterward.

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

source(file.path(script_dir, "AOC_stats_glmm_fit_trials.R"))

message("=== AOC GLMMs (trials): n-back ===")
results <- run_paper_glmm_trials(NBACK_TASK_CONFIG_TRIALS)
save_glmm_results_trials(results)
message("Done. Run AOC_stats_glmm_export_trials.R to write CSV / Word table inputs.")
