#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(glue)
})

if (isTRUE(.Platform$OS.type == "windows")) {
  base_data <- "W:/Students/Arne/AOC"
} else {
  base_data <- "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC"
}

controls_dir <- file.path(base_data, "data", "controls", "multiverse")
multiverse_dir <- file.path(base_data, "data", "multiverse")
mode_a <- "singleFFT"
mode_b <- "welch500_50"
tasks <- c("nback", "sternberg")

sanitize_mode <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
mode_a_tag <- sanitize_mode(mode_a)
mode_b_tag <- sanitize_mode(mode_b)

read_mode_file <- function(task, mode_tag) {
  path <- file.path(controls_dir, paste0("fooof_r2_", task, "_", mode_tag, ".csv"))
  if (!file.exists(path)) return(NULL)
  dat <- read_csv(path, show_col_types = FALSE)
  dat$task <- task
  dat
}

summarize_retention <- function(dat, mode_label) {
  dat %>%
    mutate(mode = mode_label) %>%
    group_by(task, mode, latency_ms, electrodes) %>%
    summarise(
      n_total = n(),
      n_low_r2 = sum(!is.na(r_squared) & r_squared < 0.60),
      pct_low_r2 = 100 * n_low_r2 / pmax(n_total, 1),
      median_r2 = median(r_squared, na.rm = TRUE),
      median_error = median(fooof_error, na.rm = TRUE),
      median_exp = median(aperiodic_exponent, na.rm = TRUE),
      median_off = median(aperiodic_offset, na.rm = TRUE),
      .groups = "drop"
    )
}

pairwise_delta <- function(dat_a, dat_b, task_name) {
  key_cols <- c("subjectID", "Condition", "Trial", "latency_ms", "electrodes")
  keep_cols <- c(key_cols, "r_squared", "fooof_error", "aperiodic_exponent", "aperiodic_offset")
  da <- dat_a %>% select(all_of(keep_cols))
  db <- dat_b %>% select(all_of(keep_cols))
  merged <- inner_join(da, db, by = key_cols, suffix = c("_a", "_b"))
  if (nrow(merged) == 0) return(tibble())
  merged %>%
    mutate(
      task = task_name,
      delta_r2 = r_squared_b - r_squared_a,
      delta_error = fooof_error_b - fooof_error_a,
      delta_exp = aperiodic_exponent_b - aperiodic_exponent_a,
      delta_off = aperiodic_offset_b - aperiodic_offset_a
    ) %>%
    group_by(task, latency_ms, electrodes) %>%
    summarise(
      n_paired = n(),
      median_delta_r2 = median(delta_r2, na.rm = TRUE),
      median_delta_error = median(delta_error, na.rm = TRUE),
      median_delta_exp = median(delta_exp, na.rm = TRUE),
      median_delta_off = median(delta_off, na.rm = TRUE),
      pct_delta_r2_positive = 100 * mean(delta_r2 > 0, na.rm = TRUE),
      .groups = "drop"
    )
}

read_multiverse_mode <- function(task, mode_tag) {
  path <- file.path(multiverse_dir, paste0("multiverse_", task, "_", mode_tag, ".csv"))
  if (!file.exists(path)) return(NULL)
  dat <- read_csv(path, show_col_types = FALSE)
  dat$task <- task
  dat
}

alpha_shift_summary <- function(dat_a, dat_b, task_name) {
  key_cols <- c(
    "subjectID", "Trial", "Condition", "electrodes", "latency_ms",
    "alpha_type", "baseline_eeg", "baseline_gaze", "gaze_measure", "fooof"
  )
  keep_cols <- c(key_cols, "alpha", "aperiodic_exponent", "aperiodic_offset")
  da <- dat_a %>% filter(fooof == "FOOOFed") %>% select(all_of(keep_cols))
  db <- dat_b %>% filter(fooof == "FOOOFed") %>% select(all_of(keep_cols))
  merged <- inner_join(da, db, by = key_cols, suffix = c("_a", "_b"))
  if (nrow(merged) == 0) return(tibble())
  merged %>%
    mutate(
      task = task_name,
      delta_alpha = alpha_b - alpha_a,
      delta_exp = aperiodic_exponent_b - aperiodic_exponent_a,
      delta_off = aperiodic_offset_b - aperiodic_offset_a
    ) %>%
    group_by(task, latency_ms, electrodes, alpha_type, baseline_eeg) %>%
    summarise(
      n_paired = n(),
      median_delta_alpha = median(delta_alpha, na.rm = TRUE),
      median_delta_exp = median(delta_exp, na.rm = TRUE),
      median_delta_off = median(delta_off, na.rm = TRUE),
      pct_alpha_delta_positive = 100 * mean(delta_alpha > 0, na.rm = TRUE),
      .groups = "drop"
    )
}

all_retention <- list()
all_delta <- list()
all_alpha_delta <- list()
missing <- c()
missing_alpha <- c()

for (task in tasks) {
  dat_a <- read_mode_file(task, mode_a_tag)
  dat_b <- read_mode_file(task, mode_b_tag)
  if (is.null(dat_a) || is.null(dat_b)) {
    missing <- c(missing, task)
    next
  }
  all_retention[[paste0(task, "_a")]] <- summarize_retention(dat_a, mode_a)
  all_retention[[paste0(task, "_b")]] <- summarize_retention(dat_b, mode_b)
  all_delta[[task]] <- pairwise_delta(dat_a, dat_b, task)

  mv_a <- read_multiverse_mode(task, mode_a_tag)
  mv_b <- read_multiverse_mode(task, mode_b_tag)
  if (is.null(mv_a) || is.null(mv_b)) {
    missing_alpha <- c(missing_alpha, task)
  } else {
    all_alpha_delta[[task]] <- alpha_shift_summary(mv_a, mv_b, task)
  }
}

if (length(all_retention) == 0) {
  msg <- glue(
    "Missing mode-specific R2 files in {controls_dir}.\n",
    "Expected: fooof_r2_<task>_{mode_a_tag}.csv and fooof_r2_<task>_{mode_b_tag}.csv\n",
    "Run AOC_multiverse_prep.m once per mode and re-run this script."
  )
  stop(msg)
}

retention_tbl <- bind_rows(all_retention)
delta_tbl <- bind_rows(all_delta)
alpha_delta_tbl <- bind_rows(all_alpha_delta)

out_retention <- file.path(controls_dir, paste0("fooof_mode_retention_compare_", mode_a_tag, "_vs_", mode_b_tag, ".csv"))
out_delta <- file.path(controls_dir, paste0("fooof_mode_delta_compare_", mode_a_tag, "_vs_", mode_b_tag, ".csv"))
out_alpha_delta <- file.path(controls_dir, paste0("fooof_mode_alpha_delta_compare_", mode_a_tag, "_vs_", mode_b_tag, ".csv"))
write_csv(retention_tbl, out_retention)
write_csv(delta_tbl, out_delta)
if (nrow(alpha_delta_tbl) > 0) {
  write_csv(alpha_delta_tbl, out_alpha_delta)
}

report_path <- file.path(controls_dir, paste0("fooof_mode_compare_", mode_a_tag, "_vs_", mode_b_tag, ".md"))
lines <- c(
  glue("# FOOOF Mode Compare: {mode_a} vs {mode_b}"),
  "",
  glue("- Controls dir: `{controls_dir}`"),
  glue("- Retention table: `{out_retention}`"),
  glue("- Paired delta table: `{out_delta}`"),
  glue("- Alpha/aperiodic paired table: `{out_alpha_delta}` (written only if mode-specific multiverse CSVs exist)"),
  ""
)

if (length(missing) > 0) {
  lines <- c(lines, "- Missing tasks:", paste0("  - ", missing), "")
}
if (length(missing_alpha) > 0) {
  lines <- c(lines, "- Missing alpha-compare multiverse files:", paste0("  - ", missing_alpha), "")
}

writeLines(lines, report_path)
message(glue("Wrote: {out_retention}"))
message(glue("Wrote: {out_delta}"))
if (nrow(alpha_delta_tbl) > 0) {
  message(glue("Wrote: {out_alpha_delta}"))
}
message(glue("Wrote: {report_path}"))
