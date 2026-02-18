# AOC Multiverse — Sternberg: Visualization
# Loads result CSVs from the analysis script and generates specification curve figures.
#
# Figures:
#   1 (_estimate):        alpha ~ gaze [sorted by estimate]
#   2 (_grouped):         alpha ~ gaze [grouped by processing dimension]
#   3 (_condition_alpha): alpha ~ condition [EEG-only]
#   4 (_interaction):     alpha ~ gaze × condition [interaction term]
#   5 (_condition_gaze):  gaze ~ condition [gaze-only]
#   6 (_aperiodic):       aperiodic component [combined forest plot]
#
# Requires: AOC_multiverse_sternberg_analysis.R to have been run first.
# Figure output: '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/multiverse'

library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)

# ========== PATHS ==========
csv_dir      <- Sys.getenv("AOC_MULTIVERSE_DIR",
                            unset = "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features")
storage_plot <- Sys.getenv("AOC_MULTIVERSE_FIGURES",
                            unset = "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/multiverse")
if (!dir.exists(storage_plot)) dir.create(storage_plot, recursive = TRUE)

# ========== THEME & AESTHETICS ==========
v_common_theme <- theme(
  plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
  axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
  axis.text.y = element_text(size = 14),
  axis.title.x = element_text(size = 14, face = "bold"),
  axis.title.y = element_text(size = 14, face = "bold"),
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 12, face = "bold"),
  plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
  plot.subtitle = element_text(size = 12, hjust = 0.5)
)

sig_colors <- c("Positive" = "#33CC66", "Negative" = "#fe0000",
                "Non-significant" = "#d1d1d1")
sig_levels <- c("Positive", "Negative", "Non-significant")

# ========== LABEL MAPPINGS ==========
group_labels <- c(
  "electrodes" = "Electrodes", "fooof" = "Spectral\nParameterization", "latency_ms" = "Latency",
  "alpha_type" = "Alpha", "gaze_measure" = "Gaze Measure",
  "baseline_eeg" = "EEG Baseline", "baseline_gaze" = "Gaze Baseline"
)

rename_opts <- function(x) {
  recode(x,
    "posterior" = "Posterior", "occipital" = "Occipital",
    "FOOOFed" = "SpecParam", "nonFOOOFed" = "No SpecParam",
    "0_500ms" = "0\u2013500 ms", "0_1000ms" = "0\u20131000 ms",
    "0_2000ms" = "0\u20132000 ms", "1000_2000ms" = "1000\u20132000 ms",
    "canonical" = "Canonical", "IAF" = "IAF",
    "scan_path_length" = "Scan Path Length", "gaze_velocity" = "Gaze Velocity",
    "microsaccades" = "Microsaccades", "BCEA" = "BCEA", "gaze_deviation" = "Gaze Deviation",
    "raw" = "Raw", "dB" = "dB", "pct_change" = "% Change"
  )
}

value_levels <- c(
  "% Change", "Raw",
  "Gaze Deviation", "Gaze Velocity",
  "IAF", "Canonical",
  "No SpecParam", "SpecParam",
  "Occipital", "Posterior",
  "1000\u20132000 ms", "0\u20132000 ms", "0\u20131000 ms", "0\u2013500 ms"
)

v_p2_group_order <- c("Latency", "Electrodes", "Spectral\nParameterization",
                       "EEG Baseline", "Alpha", "Gaze Measure", "Gaze Baseline")

elec_order <- c("posterior", "occipital")
lat_order  <- c("0_500ms", "0_1000ms", "0_2000ms", "1000_2000ms")

# ========== PANEL HELPERS ==========
make_panel_long <- function(df, x_col) {
  df %>%
    pivot_longer(cols = c(electrodes, fooof, latency_ms, alpha_type, gaze_measure,
                          baseline_eeg, baseline_gaze),
                 names_to = "Variable", values_to = "value") %>%
    mutate(
      Variable = recode(Variable, !!!group_labels),
      Variable = factor(Variable, levels = v_p2_group_order),
      value = rename_opts(value),
      value = factor(value, levels = value_levels)
    )
}

eeg_group_labels <- c(
  "electrodes" = "Electrodes", "fooof" = "Spectral\nParameterization", "latency_ms" = "Latency",
  "alpha_type" = "Alpha", "baseline_eeg" = "EEG Baseline"
)
eeg_group_order <- c("Latency", "Electrodes", "Spectral\nParameterization", "EEG Baseline", "Alpha")

make_eeg_panel_long <- function(df, x_col) {
  df %>%
    pivot_longer(cols = c(electrodes, fooof, latency_ms, alpha_type, baseline_eeg),
                 names_to = "Variable", values_to = "value") %>%
    mutate(
      Variable = recode(Variable, !!!eeg_group_labels),
      Variable = factor(Variable, levels = eeg_group_order),
      value = rename_opts(value),
      value = factor(value, levels = value_levels)
    )
}

gaze_group_labels <- c(
  "latency_ms" = "Latency", "gaze_measure" = "Gaze Measure",
  "baseline_gaze" = "Gaze Baseline"
)
gaze_group_order <- c("Latency", "Gaze Measure", "Gaze Baseline")

make_gaze_panel_long <- function(df, x_col) {
  df %>%
    pivot_longer(cols = c(latency_ms, gaze_measure, baseline_gaze),
                 names_to = "Variable", values_to = "value") %>%
    mutate(
      Variable = recode(Variable, !!!gaze_group_labels),
      Variable = factor(Variable, levels = gaze_group_order),
      value = rename_opts(value),
      value = factor(value, levels = value_levels)
    )
}

# ========== LOAD RESULT CSVs ==========
message("Loading result CSVs...")

cond_path <- file.path(csv_dir, "multiverse_sternberg_conditions_results.csv")
if (!file.exists(cond_path)) stop("Conditions CSV not found: ", cond_path, "\nRun AOC_multiverse_sternberg_analysis.R first.")
M_cond <- read.csv(cond_path, stringsAsFactors = FALSE)
M_cond$condition <- factor(M_cond$condition, levels = sig_levels)
M_cond <- M_cond %>% filter(gaze_measure %in% c("gaze_deviation", "gaze_velocity"))

# Determine highest condition from data
cond_labels_in_data <- unique(M_cond$cond_label)
cond_nums <- as.numeric(gsub("[^0-9]", "", cond_labels_in_data))
highest_label <- cond_labels_in_data[which.max(cond_nums)]
message(sprintf("Highest condition label: %s", highest_label))

# ========== FIGURES 1 & 2: GAZE -> ALPHA (highest condition) ==========
M_high <- M_cond %>% filter(cond_label == highest_label)

# --- Figure 1: sorted by estimate ---
ord_high <- M_high %>% arrange(fooof, estimate) %>% pull(.universe)
ord_df   <- data.frame(.universe = ord_high, ordered_universe = seq_along(ord_high))
M_high_est <- M_high %>% left_join(ord_df, by = ".universe")

ymax_est <- max(abs(c(M_high_est$conf.low, M_high_est$conf.high)), na.rm = TRUE) * 1.05
ylim_est <- c(-ymax_est, ymax_est)

p_curve <- ggplot(M_high_est, aes(x = ordered_universe, y = estimate, color = condition)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.3, alpha = 0.7) +
  geom_point(size = 1.2, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  scale_color_manual(values = sig_colors, name = "Significance") +
  guides(alpha = "none") +
  labs(title = expression(bold(alpha ~ "~" ~ gaze)),
       subtitle = "alpha ~ gaze_value * Condition + (1|subjectID)",
       x = "Universe", y = expression(bold("Standardized " * beta))) +
  theme_minimal() + theme(legend.position = "none") + v_common_theme +
  coord_cartesian(ylim = ylim_est)
legend_spec <- get_legend(p_curve + theme(legend.position = "bottom") + guides(alpha = "none"))

df_specs <- M_high_est %>%
  select(ordered_universe, .universe, electrodes, fooof, latency_ms, alpha_type,
         gaze_measure, baseline_eeg, baseline_gaze, condition)
df_long <- make_panel_long(df_specs, "ordered_universe")

p_panel <- ggplot(df_long, aes(x = ordered_universe, y = value, fill = condition)) +
  geom_tile() +
  scale_fill_manual(values = sig_colors, name = "Significance") +
  facet_grid(Variable ~ ., scales = "free_y", space = "free") +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0, size = 13, face = "bold", hjust = 0),
        legend.position = "bottom") +
  labs(x = "Universe", y = "Analysis Decision") + v_common_theme

p_combined <- p_curve / legend_spec / p_panel + plot_layout(heights = c(0.8, 0.1, 1.5))
ggsave(file.path(storage_plot, "AOC_multiverse_sternberg_estimate.png"),
       plot = p_combined, width = 14, height = 12, dpi = 600)
message("Saved: AOC_multiverse_sternberg_estimate.png")

# --- Figure 2: grouped by processing hierarchy ---
df_grouped <- M_high %>%
  mutate(
    .elec_ord = match(electrodes, elec_order),
    .lat_ord  = match(latency_ms, lat_order)
  ) %>%
  arrange(fooof, .lat_ord, .elec_ord, baseline_eeg, alpha_type, gaze_measure, baseline_gaze) %>%
  select(-.elec_ord, -.lat_ord)
df_grouped$grouped_universe <- seq_len(nrow(df_grouped))

ymax_grp <- max(abs(c(df_grouped$conf.low, df_grouped$conf.high)), na.rm = TRUE) * 1.05
ylim_grp <- c(-ymax_grp, ymax_grp)

p_grp_curve <- ggplot(df_grouped, aes(x = grouped_universe, y = estimate, color = condition)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.3, alpha = 0.7) +
  geom_point(size = 0.8, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  scale_color_manual(values = sig_colors, name = "Significance") +
  guides(alpha = "none") +
  labs(title = expression(bold(alpha ~ "~" ~ gaze ~ "(grouped)")),
       subtitle = "alpha ~ gaze_value * Condition + (1|subjectID)",
       x = "Universe", y = expression(bold("Standardized " * beta))) +
  theme_minimal() + theme(legend.position = "none") + v_common_theme +
  coord_cartesian(ylim = ylim_grp)
legend_grp <- get_legend(p_grp_curve + theme(legend.position = "bottom"))

df_grp_long <- make_panel_long(df_grouped, "grouped_universe")

p_grp_panel <- ggplot(df_grp_long, aes(x = grouped_universe, y = value, fill = condition)) +
  geom_tile() +
  scale_fill_manual(values = sig_colors, name = "Significance") +
  facet_grid(Variable ~ ., scales = "free_y", space = "free") +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0, size = 13, face = "bold", hjust = 0),
        legend.position = "bottom") +
  labs(x = "Universe", y = "Analysis Decision") + v_common_theme

p_grp_combined <- p_grp_curve / legend_grp / p_grp_panel +
  plot_layout(heights = c(0.8, 0.1, 1.5))
ggsave(file.path(storage_plot, "AOC_multiverse_sternberg_grouped.png"),
       plot = p_grp_combined, width = 14, height = 12, dpi = 600)
message("Saved: AOC_multiverse_sternberg_grouped.png")

# ========== FIGURE 3: CONDITION -> ALPHA (EEG-only) ==========
ca_path <- file.path(csv_dir, "multiverse_sternberg_condition_results.csv")
if (file.exists(ca_path)) {
  M_ca <- read.csv(ca_path, stringsAsFactors = FALSE)
  M_ca$condition <- factor(M_ca$condition, levels = sig_levels)

  M_ca <- M_ca %>%
    mutate(.lat_ord = match(latency_ms, lat_order),
           .elec_ord = match(electrodes, elec_order)) %>%
    arrange(.lat_ord, .elec_ord, fooof, baseline_eeg, alpha_type) %>%
    select(-.lat_ord, -.elec_ord) %>%
    mutate(ordered_universe = row_number())

  ymax_ca <- max(abs(c(M_ca$conf.low, M_ca$conf.high)), na.rm = TRUE) * 1.05
  ylim_ca <- c(-ymax_ca, ymax_ca)

  p_ca_curve <- ggplot(M_ca, aes(x = ordered_universe, y = estimate, color = condition)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.3, alpha = 0.7) +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    scale_color_manual(values = sig_colors, name = "Significance") +
    guides(alpha = "none") +
    labs(title = expression(bold(alpha ~ "~" ~ condition)),
         subtitle = "alpha ~ Condition + (1|subjectID)",
         x = "Universe", y = expression(bold("Standardized " * beta))) +
    theme_minimal() + theme(legend.position = "none") + v_common_theme +
    coord_cartesian(ylim = ylim_ca)
  legend_ca <- get_legend(p_ca_curve + theme(legend.position = "bottom"))

  df_ca_specs <- M_ca %>%
    select(ordered_universe, electrodes, fooof, latency_ms, alpha_type, baseline_eeg, condition)
  df_ca_long <- make_eeg_panel_long(df_ca_specs, "ordered_universe")

  p_ca_panel <- ggplot(df_ca_long, aes(x = ordered_universe, y = value, fill = condition)) +
    geom_tile() +
    scale_fill_manual(values = sig_colors, name = "Significance") +
    facet_grid(Variable ~ ., scales = "free_y", space = "free") +
    theme_minimal() +
    theme(strip.text.y = element_text(angle = 0, size = 13, face = "bold", hjust = 0),
          legend.position = "bottom") +
    labs(x = "Universe", y = "Analysis Decision") + v_common_theme

  p_ca_combined <- p_ca_curve / legend_ca / p_ca_panel +
    plot_layout(heights = c(0.8, 0.1, 1.2))
  ggsave(file.path(storage_plot, "AOC_multiverse_sternberg_condition_alpha.png"),
         plot = p_ca_combined, width = 14, height = 10, dpi = 600)
  message("Saved: AOC_multiverse_sternberg_condition_alpha.png")
} else {
  message("Skipping Figure 3: condition_results.csv not found.")
}

# ========== FIGURE 4: INTERACTION (gaze x condition -> alpha) ==========
int_path <- file.path(csv_dir, "multiverse_sternberg_interaction_results.csv")
if (file.exists(int_path)) {
  M_interaction <- read.csv(int_path, stringsAsFactors = FALSE)
  M_interaction$condition <- factor(M_interaction$condition, levels = sig_levels)
  M_interaction <- M_interaction %>% filter(gaze_measure %in% c("gaze_deviation", "gaze_velocity"))

  M_interaction <- M_interaction %>%
    mutate(.lat_ord = match(latency_ms, lat_order),
           .elec_ord = match(electrodes, elec_order)) %>%
    arrange(.lat_ord, .elec_ord, fooof, baseline_eeg, alpha_type, gaze_measure, baseline_gaze) %>%
    select(-.lat_ord, -.elec_ord) %>%
    mutate(ordered_universe = row_number())

  ymax_int <- max(abs(c(M_interaction$conf.low, M_interaction$conf.high)), na.rm = TRUE) * 1.05
  ylim_int <- c(-ymax_int, ymax_int)

  p_int_curve <- ggplot(M_interaction, aes(x = ordered_universe, y = estimate, color = condition)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.3, alpha = 0.7) +
    geom_point(size = 1.2, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    scale_color_manual(values = sig_colors, name = "Significance") +
    guides(alpha = "none") +
    labs(title = expression(bold(alpha ~ "~" ~ gaze ~ "\u00D7" ~ condition)),
         subtitle = "alpha ~ gaze_value * Condition + (1|subjectID)",
         x = "Universe", y = expression(bold("Standardized " * beta))) +
    theme_minimal() + theme(legend.position = "none") + v_common_theme +
    coord_cartesian(ylim = ylim_int)
  legend_int <- get_legend(p_int_curve + theme(legend.position = "bottom"))

  df_int_specs <- M_interaction %>%
    select(ordered_universe, .universe, electrodes, fooof, latency_ms, alpha_type,
           gaze_measure, baseline_eeg, baseline_gaze, condition)
  df_int_long <- make_panel_long(df_int_specs, "ordered_universe")

  p_int_panel <- ggplot(df_int_long, aes(x = ordered_universe, y = value, fill = condition)) +
    geom_tile() +
    scale_fill_manual(values = sig_colors, name = "Significance") +
    facet_grid(Variable ~ ., scales = "free_y", space = "free") +
    theme_minimal() +
    theme(strip.text.y = element_text(angle = 0, size = 13, face = "bold", hjust = 0),
          legend.position = "bottom") +
    labs(x = "Universe", y = "Analysis Decision") + v_common_theme

  p_int_combined <- p_int_curve / legend_int / p_int_panel +
    plot_layout(heights = c(0.8, 0.1, 1.5))
  ggsave(file.path(storage_plot, "AOC_multiverse_sternberg_interaction.png"),
         plot = p_int_combined, width = 14, height = 12, dpi = 600)
  message("Saved: AOC_multiverse_sternberg_interaction.png")
} else {
  message("Skipping Figure 4: interaction_results.csv not found.")
}

# ========== FIGURE 5: CONDITION -> GAZE (gaze-only) ==========
cg_path <- file.path(csv_dir, "multiverse_sternberg_condition_gaze_results.csv")
if (file.exists(cg_path)) {
  M_cg <- read.csv(cg_path, stringsAsFactors = FALSE)
  M_cg$condition <- factor(M_cg$condition, levels = sig_levels)
  M_cg <- M_cg %>% filter(gaze_measure %in% c("gaze_deviation", "gaze_velocity"))

  M_cg <- M_cg %>%
    mutate(.lat_ord = match(latency_ms, lat_order)) %>%
    arrange(.lat_ord, gaze_measure, baseline_gaze) %>%
    select(-.lat_ord) %>%
    mutate(ordered_universe = row_number())

  ymax_cg <- max(abs(c(M_cg$conf.low, M_cg$conf.high)), na.rm = TRUE) * 1.05
  ylim_cg <- c(-ymax_cg, ymax_cg)

  p_cg_curve <- ggplot(M_cg, aes(x = ordered_universe, y = estimate, color = condition)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.3, alpha = 0.7) +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    scale_color_manual(values = sig_colors, name = "Significance") +
    guides(alpha = "none") +
    labs(title = expression(bold(gaze ~ "~" ~ condition)),
         subtitle = "gaze ~ Condition + (1|subjectID)",
         x = "Universe", y = expression(bold("Standardized " * beta))) +
    theme_minimal() + theme(legend.position = "none") + v_common_theme +
    coord_cartesian(ylim = ylim_cg)
  legend_cg <- get_legend(p_cg_curve + theme(legend.position = "bottom"))

  df_cg_specs <- M_cg %>%
    select(ordered_universe, latency_ms, gaze_measure, baseline_gaze, condition)
  df_cg_long <- make_gaze_panel_long(df_cg_specs, "ordered_universe")

  p_cg_panel <- ggplot(df_cg_long, aes(x = ordered_universe, y = value, fill = condition)) +
    geom_tile() +
    scale_fill_manual(values = sig_colors, name = "Significance") +
    facet_grid(Variable ~ ., scales = "free_y", space = "free") +
    theme_minimal() +
    theme(strip.text.y = element_text(angle = 0, size = 13, face = "bold", hjust = 0),
          legend.position = "bottom") +
    labs(x = "Universe", y = "Analysis Decision") + v_common_theme

  p_cg_combined <- p_cg_curve / legend_cg / p_cg_panel +
    plot_layout(heights = c(0.8, 0.1, 0.8))
  ggsave(file.path(storage_plot, "AOC_multiverse_sternberg_condition_gaze.png"),
         plot = p_cg_combined, width = 14, height = 8, dpi = 600)
  message("Saved: AOC_multiverse_sternberg_condition_gaze.png")
} else {
  message("Skipping Figure 5: condition_gaze_results.csv not found.")
}

# ========== FIGURE 6: APERIODIC COMPONENT — COMBINED FOREST PLOT ==========
# Top: aperiodic ~ gaze — forest plot grouped by gaze measure (5 × 16 universes)
# Bottom: aperiodic ~ condition — labeled coefficient plot (8 universes)
# Columns: Exponent | Offset

ap_gaze_path <- file.path(csv_dir, "multiverse_sternberg_aperiodic_gaze_results.csv")
ap_cond_path <- file.path(csv_dir, "multiverse_sternberg_aperiodic_condition_results.csv")

has_ap_gaze <- file.exists(ap_gaze_path)
has_ap_cond <- file.exists(ap_cond_path)

if (has_ap_gaze || has_ap_cond) {

  p_forest_gaze <- p_forest_cond <- NULL

  if (has_ap_gaze) {
    M_ap_gaze <- read.csv(ap_gaze_path, stringsAsFactors = FALSE)
    M_ap_gaze$condition <- factor(M_ap_gaze$condition, levels = sig_levels)

    M_ap_gaze <- M_ap_gaze %>%
      filter(gaze_measure %in% c("gaze_deviation", "gaze_velocity")) %>%
      mutate(
        gaze_label = rename_opts(gaze_measure),
        gaze_label = factor(gaze_label, levels = c("Gaze Deviation", "Gaze Velocity")),
        aperiodic_measure = factor(aperiodic_measure, levels = c("Exponent", "Offset"))
      ) %>%
      group_by(aperiodic_measure, gaze_label) %>%
      arrange(estimate) %>%
      mutate(rank = row_number()) %>%
      ungroup()

    p_forest_gaze <- ggplot(M_ap_gaze, aes(x = estimate, y = rank, color = condition)) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
      geom_pointrange(aes(xmin = conf.low, xmax = conf.high), size = 0.3, fatten = 2) +
      scale_color_manual(values = sig_colors, name = "Significance") +
      facet_grid(gaze_label ~ aperiodic_measure) +
      labs(title = expression(bold("Aperiodic ~ gaze")),
           subtitle = "aperiodic ~ gaze_value + (1|subjectID)",
           x = expression(bold("Standardized " * beta)), y = "") +
      theme_minimal() + v_common_theme +
      theme(
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0, size = 12, face = "bold", hjust = 0),
        strip.text.x = element_text(size = 13, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
      )
  }

  if (has_ap_cond) {
    M_ap_cond <- read.csv(ap_cond_path, stringsAsFactors = FALSE)
    M_ap_cond$condition <- factor(M_ap_cond$condition, levels = sig_levels)

    M_ap_cond <- M_ap_cond %>%
      mutate(
        universe_label = paste0(rename_opts(electrodes), ", ", rename_opts(latency_ms)),
        aperiodic_measure = factor(aperiodic_measure, levels = c("Exponent", "Offset"))
      )

    label_order <- M_ap_cond %>%
      group_by(universe_label) %>%
      summarise(mean_est = mean(estimate, na.rm = TRUE), .groups = "drop") %>%
      arrange(mean_est) %>%
      pull(universe_label)
    M_ap_cond$universe_label <- factor(M_ap_cond$universe_label, levels = label_order)

    p_forest_cond <- ggplot(M_ap_cond, aes(x = estimate, y = universe_label, color = condition)) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
      geom_pointrange(aes(xmin = conf.low, xmax = conf.high), size = 0.5, fatten = 3) +
      scale_color_manual(values = sig_colors, name = "Significance") +
      facet_wrap(~aperiodic_measure, nrow = 1) +
      labs(title = expression(bold("Aperiodic ~ condition")),
           subtitle = "aperiodic ~ Condition + (1|subjectID)",
           x = expression(bold("Standardized " * beta)), y = "") +
      theme_minimal() + v_common_theme +
      theme(
        legend.position = "none",
        strip.text = element_text(size = 13, face = "bold")
      )
  }

  legend_ap <- get_legend(
    (if (!is.null(p_forest_gaze)) p_forest_gaze else p_forest_cond) +
    theme(legend.position = "bottom")
  )

  if (!is.null(p_forest_gaze) && !is.null(p_forest_cond)) {
    p_ap_combined <- p_forest_gaze / legend_ap / p_forest_cond +
      plot_layout(heights = c(1.0, 0.05, 0.35))
  } else if (!is.null(p_forest_gaze)) {
    p_ap_combined <- p_forest_gaze / legend_ap +
      plot_layout(heights = c(1.0, 0.05))
  } else {
    p_ap_combined <- p_forest_cond / legend_ap +
      plot_layout(heights = c(1.0, 0.05))
  }

  ggsave(file.path(storage_plot, "AOC_multiverse_sternberg_aperiodic.png"),
         plot = p_ap_combined, width = 12, height = 14, dpi = 600)
  message("Saved: AOC_multiverse_sternberg_aperiodic.png")

} else {
  message("Skipping Figure 6: no aperiodic result CSVs found.")
}

message("=== Sternberg multiverse VISUALIZATION complete ===")
