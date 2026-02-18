#!/usr/bin/env Rscript
# AOC Multiverse Control — Robust Z-Score Winsorization Rate
#
# Applies the same robust_z(x, winsor=3) used in the analysis scripts to
# alpha and gaze_value per universe, and reports what fraction of values
# get clipped at ±3. High rates (>10-15%) suggest the underlying
# distribution may be problematic for that universe.
#
# Reads:  multiverse_{task}.csv (trial-level)
# Writes: figures  → .../figures/controls/multiverse/winsorization/
#         data     → .../data/controls/multiverse/winsorization/

library(tidyverse)

# ========== PATHS ==========
base      <- "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC"
csv_dir   <- file.path(base, "data", "features")
fig_dir   <- file.path(base, "figures", "controls", "multiverse", "winsorization")
data_dir  <- file.path(base, "data", "controls", "multiverse", "winsorization")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

robust_z <- function(x, winsor = 3) {
  med <- median(x, na.rm = TRUE)
  ma  <- mad(x, na.rm = TRUE)
  if (is.na(ma) || ma == 0) return(rep(NaN, length(x)))
  z <- (x - med) / ma
  z[z >  winsor] <-  winsor
  z[z < -winsor] <- -winsor
  z
}

winsor_rate <- function(x, winsor = 3) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(NaN)
  med <- median(x)
  ma  <- mad(x)
  if (is.na(ma) || ma == 0) return(NaN)
  z <- (x - med) / ma
  mean(abs(z) >= winsor)
}

tasks <- c("sternberg", "nback")

# After R analysis filtering, only raw and pct_change baselines are used
bl_eeg_keep  <- c("raw", "pct_change")
bl_gaze_keep <- c("raw", "pct_change")

for (task in tasks) {
  csv_path <- file.path(csv_dir, paste0("multiverse_", task, ".csv"))
  if (!file.exists(csv_path)) { message("Skip: ", csv_path); next }

  dat <- read.csv(csv_path, stringsAsFactors = FALSE, na.strings = c("NA", "NaN", ""))
  dat <- dat %>%
    filter(!electrodes %in% c("all", "parietal"),
           baseline_eeg %in% bl_eeg_keep,
           baseline_gaze %in% bl_gaze_keep)
  message(sprintf("[%s] Loaded and filtered: %d rows.", task, nrow(dat)))

  # Compute winsorization rate per universe
  dims <- c("electrodes", "fooof", "latency_ms", "alpha_type",
            "gaze_measure", "baseline_eeg", "baseline_gaze")

  win_rates <- dat %>%
    group_by(across(all_of(dims))) %>%
    summarise(
      n_valid_alpha = sum(!is.na(alpha)),
      n_valid_gaze  = sum(!is.na(gaze_value)),
      alpha_winsor_pct = 100 * winsor_rate(alpha),
      gaze_winsor_pct  = 100 * winsor_rate(gaze_value),
      .groups = "drop"
    ) %>%
    mutate(task = task)

  write.csv(win_rates,
            file.path(data_dir, paste0("winsorization_rates_", task, ".csv")),
            row.names = FALSE)

  # --- Histogram: distribution of winsorization rates ---
  win_long <- win_rates %>%
    select(all_of(dims), alpha_winsor_pct, gaze_winsor_pct, task) %>%
    pivot_longer(cols = c(alpha_winsor_pct, gaze_winsor_pct),
                 names_to = "variable", values_to = "winsor_pct") %>%
    mutate(variable = recode(variable,
                             "alpha_winsor_pct" = "alpha",
                             "gaze_winsor_pct" = "gaze_value")) %>%
    filter(is.finite(winsor_pct))

  p1 <- ggplot(win_long, aes(x = winsor_pct, fill = variable)) +
    geom_histogram(binwidth = 1, position = "dodge", alpha = 0.8) +
    geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
    geom_vline(xintercept = 15, linetype = "dotted", color = "darkred") +
    annotate("text", x = 10, y = Inf, label = "10%", vjust = 1.5,
             color = "red", size = 3.5) +
    annotate("text", x = 15, y = Inf, label = "15%", vjust = 1.5,
             color = "darkred", size = 3.5) +
    scale_fill_manual(values = c("alpha" = "#4A90D9", "gaze_value" = "#D94A4A"),
                      name = "Variable") +
    labs(title = paste0("Winsorization rate distribution — ", task),
         subtitle = "% of values clipped at ±3 MADs per universe",
         x = "Winsorization rate (%)", y = "Number of universes") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(fig_dir, paste0("AOC_winsorization_distribution_", task, ".png")),
         plot = p1, width = 10, height = 6, dpi = 300)

  # --- Breakdown by dimension: which dimensions have highest winsorization? ---
  win_by_dim <- list()
  for (dim in dims) {
    agg <- win_long %>%
      group_by(across(all_of(dim)), variable) %>%
      summarise(mean_winsor_pct = mean(winsor_pct, na.rm = TRUE), .groups = "drop") %>%
      mutate(dimension = dim)
    names(agg)[1] <- "value"
    agg$value <- as.character(agg$value)
    win_by_dim[[dim]] <- agg
  }
  win_by_dim_df <- bind_rows(win_by_dim)

  p2 <- ggplot(win_by_dim_df, aes(x = value, y = mean_winsor_pct, fill = variable)) +
    geom_col(position = "dodge", width = 0.7) +
    facet_wrap(~dimension, scales = "free_x", ncol = 2) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "red", linewidth = 0.4) +
    scale_fill_manual(values = c("alpha" = "#4A90D9", "gaze_value" = "#D94A4A"),
                      name = "Variable") +
    labs(title = paste0("Mean winsorization rate per dimension — ", task),
         x = "", y = "Mean % clipped at ±3") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold"))

  ggsave(file.path(fig_dir, paste0("AOC_winsorization_by_dimension_", task, ".png")),
         plot = p2, width = 14, height = 12, dpi = 300)

  # Report worst universes
  worst <- win_rates %>%
    filter(is.finite(alpha_winsor_pct) | is.finite(gaze_winsor_pct)) %>%
    mutate(max_winsor = pmax(alpha_winsor_pct, gaze_winsor_pct, na.rm = TRUE)) %>%
    arrange(desc(max_winsor)) %>%
    head(20)
  write.csv(worst,
            file.path(data_dir, paste0("worst_winsorization_universes_", task, ".csv")),
            row.names = FALSE)

  n_high <- sum(win_rates$alpha_winsor_pct > 15 | win_rates$gaze_winsor_pct > 15,
                na.rm = TRUE)
  message(sprintf("[%s] %d universes with >15%% winsorization. Saved.", task, n_high))
}

message("=== Winsorization Control DONE ===")
