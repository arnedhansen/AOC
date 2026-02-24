#!/usr/bin/env Rscript
# AOC Multiverse Control — Aperiodic/Gaze Fit-Quality Sensitivity
#
# Purpose:
#   Tests whether aperiodic-gaze associations persist after controlling FOOOF fit
#   quality (R^2), and whether gaze predicts fit quality itself.
#
# Checks (trial-level, per analysis universe):
#   1) aperiodic ~ gaze                        (baseline)
#   2) aperiodic ~ gaze + r2                   (R^2 covariate control)
#   3) aperiodic ~ gaze * r2                   (interaction / moderation)
#   4) aperiodic ~ gaze on strict R^2 filters  (threshold sensitivity)
#   5) r2 ~ gaze                               (artifact pathway check)
#
# Reads:
#   multiverse_{task}.csv      from AOC_MULTIVERSE_DIR
#   fooof_r2_{task}.csv        from AOC_R2_DIR
#
# Writes:
#   data/controls/multiverse/aperiodic_fit_quality/
#     - aperiodic_fit_quality_{task}_detail.csv
#     - aperiodic_fit_quality_{task}_summary.csv
#     - aperiodic_fit_quality_{task}_join_qc.csv

library(tidyverse)
library(lme4)
library(lmerTest)
library(broom.mixed)

# ========== PATHS ==========
base <- Sys.getenv("AOC_BASE_DIR", unset = "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC")
csv_dir <- Sys.getenv("AOC_MULTIVERSE_DIR", unset = file.path(base, "data", "multiverse"))
r2_dir <- Sys.getenv("AOC_R2_DIR", unset = file.path(base, "data", "controls", "multiverse"))
data_dir <- file.path(base, "data", "controls", "multiverse", "aperiodic_fit_quality")
fig_dir <- file.path(base, "figures", "controls", "multiverse", "aperiodic_fit_quality")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ========== SETTINGS ==========
tasks <- c("sternberg", "nback")
strict_r2_thresholds <- suppressWarnings(
  as.numeric(trimws(unlist(strsplit(Sys.getenv("AOC_STRICT_R2_THRESHOLDS", unset = "0.90,0.95"), ","))))
)
strict_r2_thresholds <- strict_r2_thresholds[is.finite(strict_r2_thresholds) & strict_r2_thresholds >= 0 & strict_r2_thresholds <= 1]
if (length(strict_r2_thresholds) == 0) strict_r2_thresholds <- c(0.90, 0.95)
strict_r2_thresholds <- sort(unique(strict_r2_thresholds))

# Match main plotting order where possible.
elec_order <- c("posterior", "occipital")
lat_order <- c("0_500ms", "0_1000ms", "0_2000ms", "1000_2000ms")
gaze_pref_order <- c("gaze_deviation", "gaze_velocity", "scan_path_length", "BCEA")

# ========== HELPERS ==========
robust_z <- function(x, winsor = 3) {
  med <- median(x, na.rm = TRUE)
  ma <- mad(x, na.rm = TRUE)
  if (is.na(ma) || ma == 0) return(rep(NaN, length(x)))
  z <- (x - med) / ma
  z[z > winsor] <- winsor
  z[z < -winsor] <- -winsor
  z
}

run_lmer <- function(formula_obj, data_obj) {
  tryCatch(
    suppressWarnings(
      suppressMessages(
        lmer(formula_obj, data = data_obj, control = lmerControl(optimizer = "bobyqa"))
      )
    ),
    error = function(e) NULL
  )
}

safe_term <- function(fit, term_name) {
  if (is.null(fit)) return(tibble())
  out <- broom.mixed::tidy(fit, conf.int = TRUE) %>% filter(term == term_name)
  if (nrow(out) == 0) return(tibble())
  out
}

sig_label <- function(estimate, p_value) {
  case_when(
    p_value < 0.05 & estimate > 0 ~ "Positive",
    p_value < 0.05 & estimate < 0 ~ "Negative",
    TRUE ~ "Non-significant"
  )
}

for (task in tasks) {
  message(sprintf("\n=== Aperiodic Fit-Quality Control: %s ===", toupper(task)))

  csv_path <- file.path(csv_dir, paste0("multiverse_", task, ".csv"))
  r2_path <- file.path(r2_dir, paste0("fooof_r2_", task, ".csv"))
  if (!file.exists(csv_path)) {
    message("Skip: missing multiverse CSV: ", csv_path)
    next
  }
  if (!file.exists(r2_path)) {
    message("Skip: missing R^2 CSV: ", r2_path)
    next
  }

  dat <- read.csv(csv_path, stringsAsFactors = FALSE, na.strings = c("NA", "NaN", ""))
  r2 <- read.csv(r2_path, stringsAsFactors = FALSE, na.strings = c("NA", "NaN", ""))

  # Standard filters aligned with primary pipeline.
  if ("gaze_measure" %in% names(dat)) {
    dat <- dat %>% filter(gaze_measure != "gaze_density" | is.na(gaze_measure))
  }
  dat <- dat %>%
    filter(!electrodes %in% c("all", "parietal")) %>%
    filter(baseline_gaze != "dB", baseline_eeg != "dB") %>%
    filter(fooof == "FOOOFed") %>%
    filter(!is.na(aperiodic_offset), !is.na(aperiodic_exponent))

  if (nrow(dat) == 0) {
    message("Skip: no usable FOOOFed aperiodic rows after filtering.")
    next
  }

  join_cols <- c("subjectID", "Condition", "Trial", "latency_ms", "electrodes")
  if (!all(join_cols %in% names(dat)) || !all(c(join_cols, "r_squared") %in% names(r2))) {
    message("Skip: join columns missing in multiverse or R^2 CSV.")
    next
  }

  # Join robustly on numeric subjectID plus task keys.
  dat$subjectID_num <- suppressWarnings(as.numeric(as.character(dat$subjectID)))
  r2$subjectID_num <- suppressWarnings(as.numeric(as.character(r2$subjectID)))
  dat$Condition <- as.factor(dat$Condition)
  r2$Condition <- as.factor(r2$Condition)

  dat_qc <- dat %>%
    left_join(
      r2 %>% select(subjectID_num, Condition, Trial, latency_ms, electrodes, r_squared),
      by = c("subjectID_num", "Condition", "Trial", "latency_ms", "electrodes")
    )

  join_qc <- tibble(
    task = task,
    rows_total = nrow(dat_qc),
    rows_with_r2 = sum(is.finite(dat_qc$r_squared)),
    rows_missing_r2 = sum(!is.finite(dat_qc$r_squared)),
    pct_missing_r2 = 100 * sum(!is.finite(dat_qc$r_squared)) / max(nrow(dat_qc), 1),
    subjects = n_distinct(dat_qc$subjectID),
    universes = n_distinct(dat_qc$universe_id)
  )
  write.csv(join_qc, file.path(data_dir, paste0("aperiodic_fit_quality_", task, "_join_qc.csv")), row.names = FALSE)

  dat_model <- dat_qc %>%
    filter(is.finite(r_squared)) %>%
    select(
      subjectID, Condition, Trial, electrodes, latency_ms, gaze_measure, baseline_gaze,
      aperiodic_offset, aperiodic_exponent, gaze_value, r_squared
    ) %>%
    distinct()

  if (nrow(dat_model) == 0) {
    message("Skip: all rows missing matched R^2 after join.")
    next
  }

  dat_model$subjectID <- as.factor(dat_model$subjectID)

  specs <- dat_model %>%
    distinct(electrodes, latency_ms, gaze_measure, baseline_gaze) %>%
    mutate(
      .elec_ord = match(electrodes, elec_order),
      .lat_ord = match(latency_ms, lat_order),
      .gaze_ord = match(gaze_measure, gaze_pref_order)
    ) %>%
    mutate(
      .elec_ord = ifelse(is.na(.elec_ord), 99, .elec_ord),
      .lat_ord = ifelse(is.na(.lat_ord), 99, .lat_ord),
      .gaze_ord = ifelse(is.na(.gaze_ord), 99, .gaze_ord)
    ) %>%
    arrange(.lat_ord, .elec_ord, .gaze_ord, baseline_gaze, gaze_measure) %>%
    select(-.elec_ord, -.lat_ord, -.gaze_ord) %>%
    mutate(spec_id = row_number())

  detail_rows <- vector("list", length = 0)
  k <- 1

  for (i in seq_len(nrow(specs))) {
    s <- specs[i, ]
    ds <- dat_model %>%
      filter(
        electrodes == s$electrodes,
        latency_ms == s$latency_ms,
        gaze_measure == s$gaze_measure,
        baseline_gaze == s$baseline_gaze
      ) %>%
      filter(complete.cases(subjectID, gaze_value, r_squared, aperiodic_exponent, aperiodic_offset))

    if (nrow(ds) < 10) next

    ds <- ds %>%
      mutate(
        gaze_z = robust_z(gaze_value),
        r2_z = robust_z(r_squared),
        exponent_z = robust_z(aperiodic_exponent),
        offset_z = robust_z(aperiodic_offset)
      ) %>%
      filter(is.finite(gaze_z), is.finite(r2_z), is.finite(exponent_z), is.finite(offset_z))

    if (nrow(ds) < 10) next

    # 1) Artifact pathway check: does gaze predict fit quality?
    fit_r2 <- run_lmer(r2_z ~ gaze_z + (1 | subjectID), ds)
    t_r2 <- safe_term(fit_r2, "gaze_z")
    if (nrow(t_r2) > 0) {
      detail_rows[[k]] <- t_r2 %>%
        mutate(
          task = task,
          spec_id = s$spec_id,
          electrodes = s$electrodes,
          latency_ms = s$latency_ms,
          gaze_measure = s$gaze_measure,
          baseline_gaze = s$baseline_gaze,
          outcome = "r_squared",
          check = "r2_from_gaze",
          threshold = NA_real_,
          n_rows = nrow(ds),
          n_subjects = n_distinct(ds$subjectID)
        )
      k <- k + 1
    }

    for (outcome_name in c("exponent_z", "offset_z")) {
      measure_label <- ifelse(outcome_name == "exponent_z", "Exponent", "Offset")

      # 2) Baseline: aperiodic ~ gaze
      fit_base <- run_lmer(stats::as.formula(paste0(outcome_name, " ~ gaze_z + (1 | subjectID)")), ds)
      t_base <- safe_term(fit_base, "gaze_z")
      if (nrow(t_base) > 0) {
        detail_rows[[k]] <- t_base %>%
          mutate(
            task = task,
            spec_id = s$spec_id,
            electrodes = s$electrodes,
            latency_ms = s$latency_ms,
            gaze_measure = s$gaze_measure,
            baseline_gaze = s$baseline_gaze,
            outcome = measure_label,
            check = "baseline_gaze",
            threshold = NA_real_,
            n_rows = nrow(ds),
            n_subjects = n_distinct(ds$subjectID)
          )
        k <- k + 1
      }

      # 3) Add R^2 covariate: aperiodic ~ gaze + r2
      fit_cov <- run_lmer(stats::as.formula(paste0(outcome_name, " ~ gaze_z + r2_z + (1 | subjectID)")), ds)
      t_cov <- safe_term(fit_cov, "gaze_z")
      if (nrow(t_cov) > 0) {
        detail_rows[[k]] <- t_cov %>%
          mutate(
            task = task,
            spec_id = s$spec_id,
            electrodes = s$electrodes,
            latency_ms = s$latency_ms,
            gaze_measure = s$gaze_measure,
            baseline_gaze = s$baseline_gaze,
            outcome = measure_label,
            check = "gaze_plus_r2",
            threshold = NA_real_,
            n_rows = nrow(ds),
            n_subjects = n_distinct(ds$subjectID)
          )
        k <- k + 1
      }

      # 3b) Residualization check: residual(outcome ~ r2) ~ gaze
      resid_fit <- tryCatch(lm(stats::as.formula(paste0(outcome_name, " ~ r2_z")), data = ds), error = function(e) NULL)
      if (!is.null(resid_fit)) {
        ds_resid <- ds %>% mutate(outcome_resid = residuals(resid_fit))
        if (sum(is.finite(ds_resid$outcome_resid)) >= 10) {
          fit_resid <- run_lmer(outcome_resid ~ gaze_z + (1 | subjectID), ds_resid)
          t_resid <- safe_term(fit_resid, "gaze_z")
          if (nrow(t_resid) > 0) {
            detail_rows[[k]] <- t_resid %>%
              mutate(
                task = task,
                spec_id = s$spec_id,
                electrodes = s$electrodes,
                latency_ms = s$latency_ms,
                gaze_measure = s$gaze_measure,
                baseline_gaze = s$baseline_gaze,
                outcome = measure_label,
                check = "gaze_on_r2_residualized_outcome",
                threshold = NA_real_,
                n_rows = nrow(ds_resid),
                n_subjects = n_distinct(ds_resid$subjectID)
              )
            k <- k + 1
          }
        }
      }

      # 4) Interaction: aperiodic ~ gaze * r2
      fit_int <- run_lmer(stats::as.formula(paste0(outcome_name, " ~ gaze_z * r2_z + (1 | subjectID)")), ds)
      t_int_gaze <- safe_term(fit_int, "gaze_z")
      if (nrow(t_int_gaze) > 0) {
        detail_rows[[k]] <- t_int_gaze %>%
          mutate(
            task = task,
            spec_id = s$spec_id,
            electrodes = s$electrodes,
            latency_ms = s$latency_ms,
            gaze_measure = s$gaze_measure,
            baseline_gaze = s$baseline_gaze,
            outcome = measure_label,
            check = "gaze_in_gaze_x_r2",
            threshold = NA_real_,
            n_rows = nrow(ds),
            n_subjects = n_distinct(ds$subjectID)
          )
        k <- k + 1
      }
      t_int_term <- safe_term(fit_int, "gaze_z:r2_z")
      if (nrow(t_int_term) > 0) {
        detail_rows[[k]] <- t_int_term %>%
          mutate(
            task = task,
            spec_id = s$spec_id,
            electrodes = s$electrodes,
            latency_ms = s$latency_ms,
            gaze_measure = s$gaze_measure,
            baseline_gaze = s$baseline_gaze,
            outcome = measure_label,
            check = "gaze_x_r2_interaction",
            threshold = NA_real_,
            n_rows = nrow(ds),
            n_subjects = n_distinct(ds$subjectID)
          )
        k <- k + 1
      }

      # 5) Strict threshold sensitivity for baseline gaze slope.
      for (thr in strict_r2_thresholds) {
        ds_thr <- ds %>% filter(r_squared >= thr)
        if (nrow(ds_thr) < 10) next
        ds_thr <- ds_thr %>%
          mutate(
            gaze_z = robust_z(gaze_value),
            outcome_z = if (outcome_name == "exponent_z") robust_z(aperiodic_exponent) else robust_z(aperiodic_offset)
          ) %>%
          filter(is.finite(gaze_z), is.finite(outcome_z))
        if (nrow(ds_thr) < 10) next

        fit_thr <- run_lmer(outcome_z ~ gaze_z + (1 | subjectID), ds_thr)
        t_thr <- safe_term(fit_thr, "gaze_z")
        if (nrow(t_thr) > 0) {
          detail_rows[[k]] <- t_thr %>%
            mutate(
              task = task,
              spec_id = s$spec_id,
              electrodes = s$electrodes,
              latency_ms = s$latency_ms,
              gaze_measure = s$gaze_measure,
              baseline_gaze = s$baseline_gaze,
              outcome = measure_label,
              check = "baseline_gaze_strict_r2",
              threshold = thr,
              n_rows = nrow(ds_thr),
              n_subjects = n_distinct(ds_thr$subjectID)
            )
          k <- k + 1
        }
      }
    }
  }

  detail <- bind_rows(detail_rows)
  if (nrow(detail) == 0) {
    message("No successful models for task: ", task)
    next
  }

  detail <- detail %>%
    mutate(
      significance = sig_label(estimate, p.value),
      threshold = ifelse(is.na(threshold), NA, as.numeric(threshold))
    ) %>%
    select(
      task, outcome, check, term, threshold, estimate, std.error, conf.low, conf.high, statistic, p.value,
      significance, spec_id, electrodes, latency_ms, gaze_measure, baseline_gaze, n_rows, n_subjects
    )

  summary_tbl <- detail %>%
    group_by(task, outcome, check, term, threshold) %>%
    summarise(
      n_models = n(),
      median_estimate = median(estimate, na.rm = TRUE),
      mean_estimate = mean(estimate, na.rm = TRUE),
      q25_estimate = quantile(estimate, 0.25, na.rm = TRUE),
      q75_estimate = quantile(estimate, 0.75, na.rm = TRUE),
      pct_positive_p05 = 100 * mean(p.value < 0.05 & estimate > 0, na.rm = TRUE),
      pct_negative_p05 = 100 * mean(p.value < 0.05 & estimate < 0, na.rm = TRUE),
      pct_nonsig = 100 * mean(!(p.value < 0.05), na.rm = TRUE),
      median_n_rows = median(n_rows, na.rm = TRUE),
      median_n_subjects = median(n_subjects, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(outcome, check, threshold, term)

  detail_path <- file.path(data_dir, paste0("aperiodic_fit_quality_", task, "_detail.csv"))
  summary_path <- file.path(data_dir, paste0("aperiodic_fit_quality_", task, "_summary.csv"))
  write.csv(detail, detail_path, row.names = FALSE)
  write.csv(summary_tbl, summary_path, row.names = FALSE)

  # ========== VISUALIZATION ==========
  check_labels <- c(
    "baseline_gaze" = "Aperiodic ~ gaze",
    "gaze_plus_r2" = "Aperiodic ~ gaze + R2",
    "gaze_on_r2_residualized_outcome" = "Residualized (Aperiodic~R2) ~ gaze",
    "gaze_in_gaze_x_r2" = "Aperiodic ~ gaze * R2 (gaze term)",
    "gaze_x_r2_interaction" = "Aperiodic ~ gaze * R2 (gaze:R2)",
    "r2_from_gaze" = "R2 ~ gaze",
    "baseline_gaze_strict_r2" = "Aperiodic ~ gaze (strict R2)"
  )

  plot_df <- summary_tbl %>%
    mutate(
      check_label = recode(check, !!!check_labels),
      outcome = factor(outcome, levels = c("Exponent", "Offset", "r_squared"))
    )

  # Figure 1: model check comparison (non-threshold checks)
  p1_df <- plot_df %>%
    filter(check != "baseline_gaze_strict_r2") %>%
    mutate(check_label = factor(check_label, levels = unique(check_label)))

  if (nrow(p1_df) > 0) {
    p_checks <- ggplot(p1_df, aes(x = median_estimate, y = check_label)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.3) +
      geom_segment(aes(x = q25_estimate, xend = q75_estimate, yend = check_label),
                   linewidth = 1.0, color = "#4C78A8", alpha = 0.9) +
      geom_point(aes(size = pct_positive_p05, color = pct_nonsig), alpha = 0.9) +
      scale_color_gradient(low = "#33CC66", high = "#bfbfbf", name = "% non-significant") +
      scale_size_continuous(range = c(2.5, 6), name = "% positive p<.05") +
      facet_wrap(~ outcome, scales = "free_y") +
      labs(
        title = paste0("Aperiodic/Gaze fit-quality checks - ", toupper(task)),
        subtitle = "Point: median beta, line: IQR, size: % positive significant",
        x = "Median standardized beta",
        y = "Model check"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom"
      )

    fig_checks <- file.path(fig_dir, paste0("AOC_aperiodic_fit_quality_", task, "_check_summary.png"))
    ggsave(fig_checks, p_checks, width = 12, height = 7, dpi = 300, bg = "white")
    message("Saved: ", basename(fig_checks))
  }

  # Figure 2: strict R2 threshold sensitivity
  p2_df <- plot_df %>%
    filter(check == "baseline_gaze_strict_r2", outcome %in% c("Exponent", "Offset"), is.finite(threshold))

  if (nrow(p2_df) > 0) {
    p_strict <- ggplot(p2_df, aes(x = threshold, y = median_estimate, color = outcome, group = outcome)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.3) +
      geom_line(linewidth = 1.0) +
      geom_point(aes(size = pct_positive_p05), alpha = 0.9) +
      geom_errorbar(aes(ymin = q25_estimate, ymax = q75_estimate), width = 0.01, alpha = 0.8) +
      scale_color_manual(values = c("Exponent" = "#4C78A8", "Offset" = "#E45756")) +
      scale_size_continuous(range = c(2.5, 6), name = "% positive p<.05") +
      labs(
        title = paste0("Strict R2 sensitivity - ", toupper(task)),
        subtitle = "Aperiodic ~ gaze after hard R2 thresholds",
        x = "R2 threshold",
        y = "Median standardized beta",
        color = "Outcome"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "bottom"
      )

    fig_strict <- file.path(fig_dir, paste0("AOC_aperiodic_fit_quality_", task, "_strict_r2.png"))
    ggsave(fig_strict, p_strict, width = 8, height = 6, dpi = 300, bg = "white")
    message("Saved: ", basename(fig_strict))
  }

  message(sprintf("Saved: %s (%d rows)", basename(detail_path), nrow(detail)))
  message(sprintf("Saved: %s (%d rows)", basename(summary_path), nrow(summary_tbl)))
}

message("=== Aperiodic fit-quality control DONE ===")
