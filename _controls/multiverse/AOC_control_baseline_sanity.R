#!/usr/bin/env Rscript
# AOC Multiverse Control — Baseline Correction Sanity
#
# Checks for Inf, -Inf, and extreme values in dB and pct_change baselines.
# The R analysis drops dB rows, but pct_change rows are kept — if near-zero
# baselines produce extreme values, the robust z-scoring should handle them,
# but it's useful to know the rate and magnitude.
#
# Reads:  multiverse_{task}.csv (trial-level, includes ALL baselines)
# Writes: figures  → .../figures/controls/multiverse/baseline_sanity/
#         data     → .../data/controls/multiverse/baseline_sanity/

library(tidyverse)

# ========== PATHS ==========
base      <- "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC"
csv_dir   <- file.path(base, "data", "features")
fig_dir   <- file.path(base, "figures", "controls", "multiverse", "baseline_sanity")
data_dir  <- file.path(base, "data", "controls", "multiverse", "baseline_sanity")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

tasks <- c("sternberg", "nback")

for (task in tasks) {
  csv_path <- file.path(csv_dir, paste0("multiverse_", task, ".csv"))
  if (!file.exists(csv_path)) { message("Skip: ", csv_path); next }

  dat <- read.csv(csv_path, stringsAsFactors = FALSE, na.strings = c("NA", "NaN", ""))
  dat <- dat %>% filter(!electrodes %in% c("all", "parietal"))
  message(sprintf("[%s] Loaded: %d rows.", task, nrow(dat)))

  # --- EEG baseline comparison: raw vs dB vs pct_change ---
  eeg_bl <- dat %>%
    filter(!is.na(alpha)) %>%
    group_by(baseline_eeg) %>%
    summarise(
      n            = n(),
      n_inf        = sum(is.infinite(alpha)),
      n_nan        = sum(is.nan(alpha)),
      mean_alpha   = mean(alpha, na.rm = TRUE),
      sd_alpha     = sd(alpha, na.rm = TRUE),
      min_alpha    = min(alpha, na.rm = TRUE),
      max_alpha    = max(alpha, na.rm = TRUE),
      q01          = quantile(alpha, 0.01, na.rm = TRUE),
      q99          = quantile(alpha, 0.99, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(task = task, variable = "alpha")

  # --- Gaze baseline comparison: raw vs dB vs pct_change ---
  gaze_bl <- dat %>%
    filter(!is.na(gaze_value)) %>%
    group_by(baseline_gaze) %>%
    summarise(
      n            = n(),
      n_inf        = sum(is.infinite(gaze_value)),
      n_nan        = sum(is.nan(gaze_value)),
      mean_val     = mean(gaze_value, na.rm = TRUE),
      sd_val       = sd(gaze_value, na.rm = TRUE),
      min_val      = min(gaze_value, na.rm = TRUE),
      max_val      = max(gaze_value, na.rm = TRUE),
      q01          = quantile(gaze_value, 0.01, na.rm = TRUE),
      q99          = quantile(gaze_value, 0.99, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(task = task, variable = "gaze_value")

  write.csv(eeg_bl,
            file.path(data_dir, paste0("baseline_stats_eeg_", task, ".csv")),
            row.names = FALSE)
  write.csv(gaze_bl,
            file.path(data_dir, paste0("baseline_stats_gaze_", task, ".csv")),
            row.names = FALSE)

  # --- Density plots: alpha distribution by EEG baseline ---
  alpha_for_plot <- dat %>%
    filter(!is.na(alpha), is.finite(alpha)) %>%
    select(baseline_eeg, alpha)

  p1 <- ggplot(alpha_for_plot, aes(x = alpha, fill = baseline_eeg)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~baseline_eeg, scales = "free") +
    labs(title = paste0("Alpha distribution by EEG baseline — ", task),
         x = "Alpha power", y = "Density") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"), legend.position = "none")

  ggsave(file.path(fig_dir, paste0("AOC_baseline_alpha_density_", task, ".png")),
         plot = p1, width = 14, height = 5, dpi = 300)

  # --- Density plots: gaze distribution by gaze baseline ---
  gaze_for_plot <- dat %>%
    filter(!is.na(gaze_value), is.finite(gaze_value)) %>%
    select(baseline_gaze, gaze_value)

  p2 <- ggplot(gaze_for_plot, aes(x = gaze_value, fill = baseline_gaze)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~baseline_gaze, scales = "free") +
    labs(title = paste0("Gaze distribution by gaze baseline — ", task),
         x = "Gaze value", y = "Density") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"), legend.position = "none")

  ggsave(file.path(fig_dir, paste0("AOC_baseline_gaze_density_", task, ".png")),
         plot = p2, width = 14, height = 5, dpi = 300)

  # --- Extreme value check for pct_change (the baseline actually used in analysis) ---
  pct_alpha <- dat %>%
    filter(baseline_eeg == "pct_change", !is.na(alpha), is.finite(alpha))
  pct_gaze <- dat %>%
    filter(baseline_gaze == "pct_change", !is.na(gaze_value), is.finite(gaze_value))

  extreme_thresh <- 1000  # flag values > 1000% change as extreme

  extreme_report <- tibble(
    task = task,
    alpha_pct_n_total    = nrow(pct_alpha),
    alpha_pct_n_extreme  = sum(abs(pct_alpha$alpha) > extreme_thresh),
    alpha_pct_extreme_rate = 100 * sum(abs(pct_alpha$alpha) > extreme_thresh) / max(nrow(pct_alpha), 1),
    gaze_pct_n_total     = nrow(pct_gaze),
    gaze_pct_n_extreme   = sum(abs(pct_gaze$gaze_value) > extreme_thresh),
    gaze_pct_extreme_rate = 100 * sum(abs(pct_gaze$gaze_value) > extreme_thresh) / max(nrow(pct_gaze), 1)
  )

  write.csv(extreme_report,
            file.path(data_dir, paste0("baseline_extreme_values_", task, ".csv")),
            row.names = FALSE)

  message(sprintf("[%s] pct_change extreme (>%d): alpha=%d (%.2f%%), gaze=%d (%.2f%%)",
                  task, extreme_thresh,
                  extreme_report$alpha_pct_n_extreme, extreme_report$alpha_pct_extreme_rate,
                  extreme_report$gaze_pct_n_extreme, extreme_report$gaze_pct_extreme_rate))
}

message("=== Baseline Sanity Control DONE ===")
