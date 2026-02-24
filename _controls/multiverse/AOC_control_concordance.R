#!/usr/bin/env Rscript
# AOC Multiverse Control — Trial-Level vs Subject-Level Concordance
#
# Compares LMM effect estimates between trial-level and subject-level
# multiverse analyses. Matched universes (same branch values) should
# produce broadly consistent results; systematic divergences may indicate
# problems with aggregation or FOOOF averaging.
#
# Reads:  multiverse_{task}_results.csv           (trial-level)
#         multiverse_{task}_subject_results.csv    (subject-level)
#         multiverse_{task}_condition_results.csv
#         multiverse_{task}_subject_condition_results.csv
# Writes: figures  → .../figures/controls/multiverse/concordance/
#         data     → .../data/controls/multiverse/concordance/

library(tidyverse)

# ========== PATHS ==========
base      <- "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC"
csv_dir   <- file.path(base, "data", "multiverse")
fig_dir   <- file.path(base, "figures", "controls", "multiverse", "concordance")
data_dir  <- file.path(base, "data", "controls", "multiverse", "concordance")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

tasks <- c("sternberg", "nback")

# Branch columns used for matching universes across levels
branch_cols_main <- c("electrodes", "fooof", "latency_ms", "alpha_type",
                      "gaze_measure", "baseline_eeg", "baseline_gaze")
branch_cols_eeg  <- c("electrodes", "fooof", "latency_ms", "alpha_type", "baseline_eeg")

compare_levels <- function(trial_path, subject_path, branch_cols, term_filter, label, task) {
  if (!file.exists(trial_path) || !file.exists(subject_path)) {
    message(sprintf("  Skip %s: missing CSV(s).", label))
    return(NULL)
  }

  df_trial   <- read.csv(trial_path,   stringsAsFactors = FALSE)
  df_subject <- read.csv(subject_path, stringsAsFactors = FALSE)

  if (!is.null(term_filter)) {
    df_trial   <- df_trial   %>% filter(grepl(term_filter, term))
    df_subject <- df_subject %>% filter(grepl(term_filter, term))
  }

  # Keep the specific term with largest absolute mean estimate per universe
  if (n_distinct(df_trial$term) > 1) {
    best_term <- df_trial %>%
      group_by(term) %>% summarise(m = mean(abs(estimate), na.rm = TRUE), .groups = "drop") %>%
      slice_max(m, n = 1) %>% pull(term)
    df_trial   <- df_trial   %>% filter(term == best_term)
    df_subject <- df_subject %>% filter(term == best_term)
  }

  merged <- inner_join(
    df_trial   %>% select(all_of(branch_cols), estimate, std.error, p.value) %>%
      rename(est_trial = estimate, se_trial = std.error, p_trial = p.value),
    df_subject %>% select(all_of(branch_cols), estimate, std.error, p.value) %>%
      rename(est_subject = estimate, se_subject = std.error, p_subject = p.value),
    by = branch_cols
  )

  if (nrow(merged) == 0) {
    message(sprintf("  Skip %s: no matching universes.", label))
    return(NULL)
  }

  merged_complete <- merged %>%
    filter(is.finite(est_trial), is.finite(est_subject))

  if (nrow(merged_complete) < 2) {
    message(sprintf("  Skip %s: too few complete rows for concordance metrics.", label))
    return(NULL)
  }

  r_val <- cor(merged_complete$est_trial, merged_complete$est_subject, use = "complete.obs")
  fit <- lm(est_subject ~ est_trial, data = merged_complete)
  slope <- unname(coef(fit)[["est_trial"]])
  intercept <- unname(coef(fit)[["(Intercept)"]])
  mae <- mean(abs(merged_complete$est_subject - merged_complete$est_trial), na.rm = TRUE)
  sign_agree <- 100 * mean(sign(merged_complete$est_subject) == sign(merged_complete$est_trial), na.rm = TRUE)
  mean_x <- mean(merged_complete$est_trial, na.rm = TRUE)
  mean_y <- mean(merged_complete$est_subject, na.rm = TRUE)
  var_x <- var(merged_complete$est_trial, na.rm = TRUE)
  var_y <- var(merged_complete$est_subject, na.rm = TRUE)
  cov_xy <- cov(merged_complete$est_trial, merged_complete$est_subject, use = "complete.obs")
  ccc <- (2 * cov_xy) / (var_x + var_y + (mean_x - mean_y)^2)

  merged$label <- label
  merged$task  <- task
  merged$r     <- r_val

  list(
    data = merged, r = r_val, label = label, n = nrow(merged),
    slope = slope, intercept = intercept, mae = mae,
    sign_agree = sign_agree, ccc = ccc
  )
}

for (task in tasks) {
  Task <- paste0(toupper(substr(task, 1, 1)), substr(task, 2, nchar(task)))
  message(sprintf("=== %s ===", task))

  results <- list()

  # 1. Main: alpha ~ gaze_value (highest condition slope)
  results$main <- compare_levels(
    file.path(csv_dir, paste0("multiverse_", task, "_results.csv")),
    file.path(csv_dir, paste0("multiverse_", task, "_subject_results.csv")),
    branch_cols_main, "gaze_value", "alpha ~ gaze", task
  )

  # 2. Condition: alpha ~ condition
  results$condition <- compare_levels(
    file.path(csv_dir, paste0("multiverse_", task, "_condition_results.csv")),
    file.path(csv_dir, paste0("multiverse_", task, "_subject_condition_results.csv")),
    branch_cols_eeg, "Condition", "alpha ~ condition", task
  )

  # 3. Interaction
  results$interaction <- compare_levels(
    file.path(csv_dir, paste0("multiverse_", task, "_interaction_results.csv")),
    file.path(csv_dir, paste0("multiverse_", task, "_subject_interaction_results.csv")),
    branch_cols_main, NULL, "gaze x condition", task
  )

  valid <- results[!sapply(results, is.null)]
  if (length(valid) == 0) {
    message(sprintf("  No valid comparisons for %s.", task))
    next
  }
  label_order <- c("alpha ~ condition", "alpha ~ gaze", "gaze x condition")
  valid <- valid[order(match(sapply(valid, `[[`, "label"), label_order))]

  all_merged <- bind_rows(lapply(valid, `[[`, "data"))
  all_merged$label <- factor(all_merged$label, levels = label_order)
  write.csv(all_merged,
            file.path(data_dir, paste0("concordance_", task, ".csv")),
            row.names = FALSE)

  # --- Scatterplot: trial-level vs subject-level estimates ---
  subtitle_expr <- paste(
    sapply(valid, function(v)
      sprintf("'%s: r=%.3f ('*n[universes]*'=%d)'", v$label, v$r, v$n)),
    collapse = "*'  |  '*"
  )
  metrics_df <- bind_rows(lapply(valid, function(v) {
    d <- v$data %>% filter(is.finite(est_trial), is.finite(est_subject))
    xr <- range(d$est_trial, na.rm = TRUE)
    yr <- range(d$est_subject, na.rm = TRUE)
    dx <- max(diff(xr), 1e-9)
    dy <- max(diff(yr), 1e-9)
    tibble(
      label = v$label,
      x_pos = xr[1] + 0.02 * dx,
      y_pos = yr[1] - 0.06 * dy,
      metrics_label = sprintf(
        "slope=%.2f  intercept=%.2f\nMAE=%.3f  sign agreement=%.1f%%\nCCC=%.3f",
        v$slope, v$intercept, v$mae, v$sign_agree, v$ccc
      )
    )
  }))
  metrics_df$label <- factor(metrics_df$label, levels = label_order)

  p <- ggplot(all_merged, aes(x = est_trial, y = est_subject)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_point(alpha = 0.45, size = 1.5, color = "#336699") +
    geom_smooth(method = "lm", se = TRUE, color = "firebrick", linewidth = 0.8) +
    geom_text(
      data = metrics_df,
      aes(x = x_pos, y = y_pos, label = metrics_label),
      inherit.aes = FALSE, hjust = 0, vjust = 0,
      size = 3.8, color = "grey15", lineheight = 1.10
    ) +
    facet_wrap(~label, scales = "free") +
    coord_cartesian(clip = "off") +
    labs(
      title = paste0("Trial-level vs subject-level estimates — ", Task),
      subtitle = as.expression(parse(text = subtitle_expr)),
      x = expression(bold("Trial-level " * beta)),
      y = expression(bold("Subject-level " * beta))
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          axis.title.x = element_text(margin = margin(t = 30)),
          plot.margin = margin(10, 10, 80, 10))

  ggsave(file.path(fig_dir, paste0("AOC_concordance_", task, ".png")),
        plot = p, width = 14, height = 7, dpi = 300, bg = "white")

  # --- Summary report ---
  summary_df <- tibble(
    task     = task,
    analysis = sapply(valid, `[[`, "label"),
    r        = sapply(valid, `[[`, "r"),
    n_universes = sapply(valid, `[[`, "n"),
    slope = sapply(valid, `[[`, "slope"),
    intercept = sapply(valid, `[[`, "intercept"),
    mae = sapply(valid, `[[`, "mae"),
    sign_agreement_pct = sapply(valid, `[[`, "sign_agree"),
    ccc = sapply(valid, `[[`, "ccc")
  )
  write.csv(summary_df,
            file.path(data_dir, paste0("concordance_summary_", task, ".csv")),
            row.names = FALSE)

  for (v in valid) {
    message(sprintf(
      "  %s: r = %.3f, CCC = %.3f, slope = %.2f, intercept = %.2f, MAE = %.3f, sign agreement = %.1f%% (%d universes)",
      v$label, v$r, v$ccc, v$slope, v$intercept, v$mae, v$sign_agree, v$n
    ))
  }
}

message("=== Concordance Control DONE ===")
