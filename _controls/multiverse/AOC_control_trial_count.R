#!/usr/bin/env Rscript
# AOC Multiverse Control — Trial Count Report
#
# Per subject × condition: how many valid (non-NaN) trials exist for
# alpha and gaze_value in a representative universe (raw baseline,
# canonical alpha, posterior electrodes, non-FOOOFed, 0-2000ms, SPL)?
# Also reports total trials per subject × condition across all universes.
#
# Reads:  multiverse_{task}.csv (trial-level)
# Writes: figures  → .../figures/controls/multiverse/trial_count/
#         data     → .../data/controls/multiverse/trial_count/

library(tidyverse)

# ========== PATHS ==========
base      <- "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC"
csv_dir   <- file.path(base, "data", "features")
fig_dir   <- file.path(base, "figures", "controls", "multiverse", "trial_count")
data_dir  <- file.path(base, "data", "controls", "multiverse", "trial_count")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

tasks <- c("sternberg", "nback")

for (task in tasks) {
  csv_path <- file.path(csv_dir, paste0("multiverse_", task, ".csv"))
  if (!file.exists(csv_path)) { message("Skip: ", csv_path); next }

  dat <- read.csv(csv_path, stringsAsFactors = FALSE, na.strings = c("NA", "NaN", ""))
  message(sprintf("[%s] Loaded %d rows.", task, nrow(dat)))

  # Pick one representative universe to count unique trials per subject × condition
  ref <- dat %>%
    filter(electrodes == "posterior", fooof == "nonFOOOFed",
           latency_ms == "0_2000ms", alpha_type == "canonical",
           gaze_measure == "scan_path_length",
           baseline_eeg == "raw", baseline_gaze == "raw")

  trial_counts <- ref %>%
    group_by(subjectID, Condition) %>%
    summarise(
      n_total      = n(),
      n_alpha_ok   = sum(!is.na(alpha)),
      n_gaze_ok    = sum(!is.na(gaze_value)),
      n_both_ok    = sum(!is.na(alpha) & !is.na(gaze_value)),
      .groups = "drop"
    ) %>%
    mutate(task = task)

  write.csv(trial_counts,
            file.path(data_dir, paste0("trial_counts_", task, ".csv")),
            row.names = FALSE)

  # --- Bar chart: valid trials per subject × condition ---
  trial_counts$subjectID <- factor(trial_counts$subjectID)
  trial_counts$Condition <- factor(trial_counts$Condition)

  p <- ggplot(trial_counts, aes(x = subjectID, y = n_both_ok, fill = Condition)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "red", linewidth = 0.5) +
    annotate("text", x = Inf, y = 10, label = "min = 10", hjust = 1.1, vjust = -0.5,
             color = "red", size = 3.5) +
    labs(title = paste0("Valid trials (alpha & gaze non-NaN) — ", task),
         subtitle = "Reference universe: posterior, no SpecParam, 0–2000 ms, canonical, SPL, raw baseline",
         x = "Subject", y = "Valid trials") +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(fig_dir, paste0("AOC_trial_count_", task, ".png")),
         plot = p, width = 12, height = 6, dpi = 300)

  # --- Summary: overall valid-trial counts across all gaze measures ---
  all_gaze <- dat %>%
    filter(electrodes == "posterior", fooof == "nonFOOOFed",
           latency_ms == "0_2000ms", alpha_type == "canonical",
           baseline_eeg == "raw", baseline_gaze == "raw") %>%
    group_by(subjectID, Condition, gaze_measure) %>%
    summarise(n_both_ok = sum(!is.na(alpha) & !is.na(gaze_value)), .groups = "drop") %>%
    mutate(task = task)

  p2 <- ggplot(all_gaze, aes(x = factor(subjectID), y = n_both_ok, fill = Condition)) +
    geom_col(position = "dodge", width = 0.7) +
    facet_wrap(~gaze_measure, ncol = 2) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "red", linewidth = 0.4) +
    labs(title = paste0("Valid trials by gaze measure — ", task),
         x = "Subject", y = "Valid trials") +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  ggsave(file.path(fig_dir, paste0("AOC_trial_count_by_gaze_", task, ".png")),
         plot = p2, width = 14, height = 10, dpi = 300)

  message(sprintf("[%s] Saved trial count figures and data.", task))
}

message("=== Trial Count Control DONE ===")
