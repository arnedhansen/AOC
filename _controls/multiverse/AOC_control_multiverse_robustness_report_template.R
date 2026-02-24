# AOC Control — Robustness Report Template Builder
# Creates report-ready comparison tables across:
#   - trial-level (primary) estimates
#   - subject-level (robustness) estimates
#   - FOOOF R2 threshold retention summaries
#
# Output:
#   multiverse_{task}_robustness_report_template.csv

library(tidyverse)

csv_dir <- Sys.getenv(
  "AOC_MULTIVERSE_DIR",
  unset = "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/multiverse"
)

safe_read <- function(path) {
  if (!file.exists(path)) return(tibble())
  read.csv(path, stringsAsFactors = FALSE, na.strings = c("NA", "NaN", ""))
}

extract_hypothesis <- function(df, label, level_label) {
  if (nrow(df) == 0) return(tibble())
  df %>%
    transmute(
      task,
      hypothesis = label,
      analysis_level = level_label,
      estimate,
      conf.low,
      conf.high,
      p.value,
      std.error,
      n_rows = NA_real_,
      n_subjects = NA_real_,
      n_trials = NA_real_,
      n_universes = if (".universe" %in% names(df)) n_distinct(.universe) else NA_real_,
      r2_threshold = NA_real_,
      notes = NA_character_
    )
}

for (task in c("sternberg", "nback")) {
  trial_alpha_path <- file.path(csv_dir, paste0("multiverse_", task, "_condition_results.csv"))
  trial_ap_path <- file.path(csv_dir, paste0("multiverse_", task, "_aperiodic_condition_results.csv"))
  subj_alpha_path <- file.path(csv_dir, paste0("multiverse_", task, "_subject_condition_results.csv"))
  subj_ap_path <- file.path(csv_dir, paste0("multiverse_", task, "_subject_aperiodic_condition_results.csv"))
  r2_summary_path <- file.path(csv_dir, paste0("multiverse_", task, "_fooof_r2_threshold_sensitivity.csv"))

  trial_alpha <- safe_read(trial_alpha_path)
  trial_ap <- safe_read(trial_ap_path)
  subj_alpha <- safe_read(subj_alpha_path)
  subj_ap <- safe_read(subj_ap_path)
  r2_summary <- safe_read(r2_summary_path)

  if (nrow(trial_alpha) > 0) trial_alpha$task <- task
  if (nrow(trial_ap) > 0) trial_ap$task <- task
  if (nrow(subj_alpha) > 0) subj_alpha$task <- task
  if (nrow(subj_ap) > 0) subj_ap$task <- task

  trial_cond_alpha <- extract_hypothesis(trial_alpha, "Condition -> Alpha", "trial_primary")
  subject_cond_alpha <- extract_hypothesis(subj_alpha, "Condition -> Alpha", "subject_robustness")

  trial_cond_ap <- extract_hypothesis(
    trial_ap %>% filter(aperiodic_measure %in% c("Exponent", "Offset")),
    "Condition -> Aperiodic",
    "trial_primary"
  )
  subject_cond_ap <- extract_hypothesis(
    subj_ap %>% filter(aperiodic_measure %in% c("Exponent", "Offset")),
    "Condition -> Aperiodic",
    "subject_robustness"
  )

  threshold_rows <- if (nrow(r2_summary) > 0) {
    r2_summary %>%
      transmute(
        task,
        hypothesis = "FOOOF Fit Retention",
        analysis_level = "trial_primary_qc",
        estimate = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        p.value = NA_real_,
        std.error = NA_real_,
        n_rows = rows_retained,
        n_subjects = subjects_retained,
        n_trials = trials_retained,
        n_universes = universes_retained,
        r2_threshold,
        notes = paste0(
          "Rows removed: ", rows_removed, " (", round(rows_removed_pct, 1), "%); ",
          "FOOOFed removed: ", fooof_rows_removed, " (", round(fooof_rows_removed_pct, 1), "%)"
        )
      )
  } else {
    tibble()
  }

  template <- bind_rows(
    trial_cond_alpha,
    subject_cond_alpha,
    trial_cond_ap,
    subject_cond_ap,
    threshold_rows
  ) %>%
    mutate(
      estimate_direction = case_when(
        is.na(estimate) ~ NA_character_,
        estimate > 0 ~ "increase",
        estimate < 0 ~ "decrease",
        TRUE ~ "null"
      )
    )

  out_path <- file.path(csv_dir, paste0("multiverse_", task, "_robustness_report_template.csv"))
  write.csv(template, out_path, row.names = FALSE)
  message("Saved: ", out_path)
}

message("=== Robustness report templates complete ===")
