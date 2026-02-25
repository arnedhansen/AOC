#!/usr/bin/env Rscript

# AOC Control — Subject-level robustness compact report
# Builds a compact cross-task table and visualization from
# multiverse_*_subject_robustness_summary.csv outputs.

library(tidyverse)

csv_dir <- Sys.getenv(
  "AOC_MULTIVERSE_DIR",
  unset = "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/multiverse"
)
fig_dir <- Sys.getenv(
  "AOC_MULTIVERSE_FIGURES",
  unset = "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/multiverse/subject-level"
)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

paths <- tibble(
  task = c("sternberg", "nback"),
  path = c(
    file.path(csv_dir, "multiverse_sternberg_subject_robustness_summary.csv"),
    file.path(csv_dir, "multiverse_nback_subject_robustness_summary.csv")
  )
)

missing <- paths %>% filter(!file.exists(path))
if (nrow(missing) > 0) {
  stop("Missing robustness files: ", paste(missing$path, collapse = ", "))
}

raw <- paths %>%
  mutate(dat = purrr::map(path, ~ read.csv(.x, stringsAsFactors = FALSE, na.strings = c("NA", "NaN", "")))) %>%
  select(task, dat) %>%
  tidyr::unnest(dat)

main_int_high <- raw %>%
  filter(model_family == "main_model_terms", grepl("gaze_value(_z)?:Condition", term)) %>%
  mutate(cond_num = suppressWarnings(as.numeric(gsub("[^0-9]", "", term)))) %>%
  group_by(task) %>%
  filter(cond_num == max(cond_num, na.rm = TRUE)) %>%
  ungroup() %>%
  select(task, term)

compact <- raw %>%
  left_join(main_int_high %>% mutate(keep_highest_interaction = TRUE), by = c("task", "term")) %>%
  mutate(
    effect_label = case_when(
      model_family == "main_model_terms" & term %in% c("gaze_value_z", "gaze_value") ~ "Main: gaze slope",
      model_family == "main_model_terms" & grepl("gaze_value(_z)?:(Condition|load_trend)", term) &
        !is.na(keep_highest_interaction) & keep_highest_interaction ~ "Main: gaze x load",
      model_family == "condition_to_alpha" & grepl("Condition|load_trend", term) ~ "Load -> alpha",
      model_family == "condition_to_gaze" & grepl("Condition|load_trend", term) ~ "Load -> gaze",
      model_family == "aperiodic_to_gaze" & term %in% c("gaze_value_z", "gaze_value") & aperiodic_measure == "Exponent" ~ "Aperiodic exponent <- gaze",
      model_family == "aperiodic_to_gaze" & term %in% c("gaze_value_z", "gaze_value") & aperiodic_measure == "Offset" ~ "Aperiodic offset <- gaze",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(effect_label)) %>%
  select(
    task, effect_label, model_family, term, aperiodic_measure,
    n_universes, median_beta, dominant_sign, pct_cis_crossing_zero
  ) %>%
  arrange(task, effect_label)

table_out <- file.path(csv_dir, "multiverse_subject_robustness_compact_table.csv")
write.csv(compact, table_out, row.names = FALSE)

plot_dat <- compact %>%
  mutate(
    task = factor(task, levels = c("nback", "sternberg"), labels = c("N-back", "Sternberg")),
    effect_label = factor(effect_label, levels = c(
      "Main: gaze slope",
      "Main: gaze x load",
      "Load -> alpha",
      "Load -> gaze",
      "Aperiodic exponent <- gaze",
      "Aperiodic offset <- gaze"
    ))
  )

p <- ggplot(plot_dat, aes(x = median_beta, y = effect_label, color = pct_cis_crossing_zero)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_point(size = 3) +
  facet_wrap(~ task, ncol = 1) +
  scale_color_gradient(low = "#1a9641", high = "#d7191c", name = "% CIs crossing 0") +
  labs(
    title = "Subject-level robustness: compact effect profile",
    subtitle = "Median beta (x-axis) and uncertainty burden (color) across key model families",
    x = "Median standardized beta",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

fig_out <- file.path(fig_dir, "AOC_subject_robustness_compact_profile.png")
ggsave(fig_out, p, width = 11, height = 8, dpi = 300, bg = "white")

message("Saved table: ", table_out, " (rows=", nrow(compact), ")")
message("Saved figure: ", fig_out)

