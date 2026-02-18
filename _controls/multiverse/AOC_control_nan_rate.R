#!/usr/bin/env Rscript
# AOC Multiverse Control — NaN Rate Heatmap
#
# For each analysis dimension, computes the percentage of NaN values in
# alpha and gaze_value. Produces heatmaps faceted by task showing which
# dimension combinations suffer from data loss.
#
# Reads:  multiverse_{task}.csv (trial-level)
# Writes: figures  → .../figures/controls/multiverse/NaN_rate/
#         data     → .../data/controls/multiverse/NaN_rate/

library(tidyverse)
library(scales)

# ========== PATHS ==========
base      <- "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC"
csv_dir   <- file.path(base, "data", "features")
fig_dir   <- file.path(base, "figures", "controls", "multiverse", "NaN_rate")
data_dir  <- file.path(base, "data", "controls", "multiverse", "NaN_rate")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

tasks <- c("sternberg", "nback")
dims  <- c("electrodes", "fooof", "latency_ms", "alpha_type",
           "gaze_measure", "baseline_eeg", "baseline_gaze")

for (task in tasks) {
  csv_path <- file.path(csv_dir, paste0("multiverse_", task, ".csv"))
  if (!file.exists(csv_path)) { message("Skip: ", csv_path); next }

  dat <- read.csv(csv_path, stringsAsFactors = FALSE, na.strings = c("NA", "NaN", ""))
  message(sprintf("[%s] Loaded %d rows.", task, nrow(dat)))

  # --- Per-dimension NaN rates ---
  nan_rows <- list()
  for (dim in dims) {
    if (!dim %in% names(dat)) next
    grp <- dat %>%
      group_by(across(all_of(dim))) %>%
      summarise(
        n_total    = n(),
        alpha_nan_pct      = 100 * mean(is.na(alpha)),
        gaze_value_nan_pct = 100 * mean(is.na(gaze_value)),
        .groups = "drop"
      ) %>%
      mutate(dimension = dim, task = task)
    names(grp)[1] <- "value"
    grp$value <- as.character(grp$value)
    nan_rows[[dim]] <- grp
  }
  nan_df <- bind_rows(nan_rows)

  write.csv(nan_df, file.path(data_dir, paste0("nan_rate_per_dimension_", task, ".csv")),
            row.names = FALSE)

  # --- Cross-tabulation: gaze_measure × latency_ms for gaze NaN ---
  cross <- dat %>%
    group_by(gaze_measure, latency_ms) %>%
    summarise(gaze_nan_pct = 100 * mean(is.na(gaze_value)), .groups = "drop")

  write.csv(cross, file.path(data_dir, paste0("nan_rate_gaze_x_latency_", task, ".csv")),
            row.names = FALSE)

  # --- Heatmap: gaze NaN by gaze_measure × latency ---
  p_cross <- ggplot(cross, aes(x = latency_ms, y = gaze_measure, fill = gaze_nan_pct)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.1f%%", gaze_nan_pct)), size = 4) +
    scale_fill_gradient(low = "white", high = "firebrick", name = "% NaN",
                        limits = c(0, 100)) +
    labs(title = paste0("Gaze NaN rate — ", task),
         subtitle = "gaze_measure × latency_ms",
         x = "Latency", y = "Gaze measure") +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title = element_text(face = "bold"))

  ggsave(file.path(fig_dir, paste0("AOC_nan_rate_gaze_heatmap_", task, ".png")),
         plot = p_cross, width = 10, height = 6, dpi = 300)

  # --- Bar chart: per-dimension NaN rates for alpha and gaze ---
  nan_long <- nan_df %>%
    pivot_longer(cols = c(alpha_nan_pct, gaze_value_nan_pct),
                 names_to = "variable", values_to = "nan_pct") %>%
    mutate(variable = recode(variable,
                             "alpha_nan_pct" = "alpha",
                             "gaze_value_nan_pct" = "gaze_value"))

  p_bar <- ggplot(nan_long, aes(x = value, y = nan_pct, fill = variable)) +
    geom_col(position = "dodge", width = 0.7) +
    facet_wrap(~dimension, scales = "free_x", ncol = 2) +
    scale_fill_manual(values = c("alpha" = "#4A90D9", "gaze_value" = "#D94A4A"),
                      name = "Variable") +
    labs(title = paste0("NaN rate per dimension — ", task),
         x = "", y = "% NaN") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold"))

  ggsave(file.path(fig_dir, paste0("AOC_nan_rate_per_dimension_", task, ".png")),
         plot = p_bar, width = 14, height = 12, dpi = 300)

  message(sprintf("[%s] Saved NaN rate figures and data.", task))
}

message("=== NaN Rate Control DONE ===")
