# AOC Control — Multiverse with FOOOF R² Trial Filter
# Re-runs the full multiverse analysis and visualization for both Sternberg
# and N-back, excluding any trial whose FOOOF fit R² < 0.6.
#
# The filter is applied only to FOOOFed universes (nonFOOOFed universes and
# gaze-only models are unaffected since FOOOF was not used).
#
# Input:
#   multiverse_{task}.csv          — trial-level multiverse data
#   fooof_r2_{task}.csv            — per-trial FOOOF R² values
#
# Output figures (same as main multiverse, in FOOOF control directory):
#   AOC_multiverse_{task}_estimate.png
#   AOC_multiverse_{task}_grouped.png
#   AOC_multiverse_{task}_condition_alpha.png
#   AOC_multiverse_{task}_interaction.png
#   AOC_multiverse_{task}_condition_gaze.png
#   AOC_multiverse_{task}_aperiodic_exponent_spec.png
#   AOC_multiverse_{task}_aperiodic_offset_spec.png

library(lme4)
library(lmerTest)
library(broom.mixed)
library(tidyverse)
library(multiverse)
library(ggplot2)
library(patchwork)
library(cowplot)

# ========== PATHS ==========
csv_dir      <- Sys.getenv("AOC_MULTIVERSE_DIR",
                            unset = "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/multiverse")
r2_dir       <- Sys.getenv("AOC_R2_DIR",
                            unset = "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/controls/multiverse")
storage_plot <- "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/controls/multiverse/FOOOF"
if (!dir.exists(storage_plot)) dir.create(storage_plot, recursive = TRUE)

parse_threshold_grid <- function(x) {
  vals <- suppressWarnings(as.numeric(trimws(unlist(strsplit(x, ",")))))
  vals <- vals[is.finite(vals) & vals >= 0 & vals <= 1]
  vals <- sort(unique(vals))
  if (length(vals) == 0) vals <- c(0.5, 0.6, 0.7, 0.8, 0.9)
  vals
}

R2_THRESHOLD <- suppressWarnings(as.numeric(Sys.getenv("AOC_R2_THRESHOLD", unset = "0.6")))
if (!is.finite(R2_THRESHOLD)) R2_THRESHOLD <- 0.6
R2_THRESHOLD_GRID <- parse_threshold_grid(Sys.getenv("AOC_R2_THRESHOLD_GRID", unset = "0.5,0.6,0.7,0.8,0.9"))

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
  "BCEA", "Scan Path Length", "Gaze Velocity", "Gaze Deviation",
  "IAF", "Canonical",
  "No SpecParam", "SpecParam",
  "Occipital", "Posterior",
  "1000\u20132000 ms", "0\u20132000 ms", "0\u20131000 ms", "0\u2013500 ms"
)

v_p2_group_order <- c("Latency", "Electrodes", "Spectral\nParameterization",
                       "EEG Baseline", "Alpha", "Gaze Measure", "Gaze Baseline")

elec_order <- c("posterior", "occipital")
lat_order  <- c("0_500ms", "0_1000ms", "0_2000ms", "1000_2000ms")
gaze_order <- c("gaze_deviation", "gaze_velocity", "scan_path_length", "BCEA")

# ========== HELPER FUNCTIONS ==========
robust_z <- function(x, winsor = 3) {
  med <- median(x, na.rm = TRUE)
  ma  <- mad(x, na.rm = TRUE)
  if (is.na(ma) || ma == 0) return(rep(NaN, length(x)))
  z <- (x - med) / ma
  z[z >  winsor] <-  winsor
  z[z < -winsor] <- -winsor
  z
}

add_sig <- function(df) {
  df %>% mutate(
    condition = factor(case_when(
      p.value < 0.05 & estimate > 0 ~ "Positive",
      p.value < 0.05 & estimate < 0 ~ "Negative",
      TRUE ~ "Non-significant"
    ), levels = c("Positive", "Negative", "Non-significant"))
  )
}

safe_extract <- function(results, var_name) {
  r <- results[[var_name]]
  if (is.null(r) || !is.data.frame(r) || nrow(r) == 0) return(NULL)
  r
}

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

aperiodic_group_labels <- c(
  "latency_ms" = "Latency", "electrodes" = "Electrodes",
  "gaze_measure" = "Gaze Measure", "baseline_gaze" = "Gaze Baseline"
)
aperiodic_group_order <- c("Latency", "Electrodes", "Gaze Measure", "Gaze Baseline")

make_aperiodic_panel_long <- function(df) {
  df %>%
    pivot_longer(cols = c(latency_ms, electrodes, gaze_measure, baseline_gaze),
                 names_to = "Variable", values_to = "value") %>%
    mutate(
      Variable = recode(Variable, !!!aperiodic_group_labels),
      Variable = factor(Variable, levels = aperiodic_group_order),
      value = rename_opts(value),
      value = factor(value, levels = value_levels)
    )
}

branch_cols     <- c("electrodes", "fooof", "latency_ms", "alpha_type",
                     "gaze_measure", "baseline_eeg", "baseline_gaze")
eeg_branch_cols <- c("electrodes", "fooof", "latency_ms", "alpha_type", "baseline_eeg")
gaze_branch_cols <- c("latency_ms", "gaze_measure", "baseline_gaze")
ap_gaze_branch_cols <- c("electrodes", "latency_ms", "gaze_measure", "baseline_gaze")
ap_cond_branch_cols <- c("electrodes", "latency_ms")

# ===========================================================================
#                           TASK LOOP
# ===========================================================================
for (task in c("sternberg", "nback")) {

  message(sprintf("\n========== %s (FOOOF R\u00b2 \u2265 %.2f filter) ==========\n", toupper(task), R2_THRESHOLD))

  # ========== LOAD DATA ==========
  csv_path <- file.path(csv_dir, paste0("multiverse_", task, ".csv"))
  r2_path  <- file.path(r2_dir, paste0("fooof_r2_", task, ".csv"))
  if (!file.exists(csv_path)) { message("CSV not found: ", csv_path, " — skipping."); next }
  if (!file.exists(r2_path))  { message("R\u00b2 CSV not found: ", r2_path, " — skipping."); next }

  dat <- read.csv(csv_path, stringsAsFactors = FALSE, na.strings = c("NA", "NaN", ""))
  r2  <- read.csv(r2_path,  stringsAsFactors = FALSE, na.strings = c("NA", "NaN", ""))

  dat$subjectID <- as.factor(dat$subjectID)
  dat$Condition <- as.factor(dat$Condition)

  message(sprintf("Loaded multiverse: %d rows, %d universes.", nrow(dat), n_distinct(dat$universe_id)))
  message(sprintf("Loaded R\u00b2 data: %d FOOOF fits.", nrow(r2)))

  # ========== APPLY FOOOF R² FILTER ==========
  join_cols <- c("subjectID", "Condition", "Trial", "latency_ms", "electrodes")
  r2_slim   <- r2 %>% select(all_of(join_cols), r_squared)

  # Ensure matching types
  dat$subjectID_num <- as.numeric(as.character(dat$subjectID))
  r2_slim$subjectID_num <- as.numeric(r2_slim$subjectID)
  r2_slim$Condition <- as.factor(r2_slim$Condition)

  dat_qc <- dat %>%
    left_join(r2_slim %>% select(subjectID_num, Condition, Trial, latency_ms, electrodes, r_squared),
              by = c("subjectID_num", "Condition", "Trial", "latency_ms", "electrodes"))

  # Threshold sensitivity summary (retention profile across R² cutoffs).
  threshold_summary <- bind_rows(lapply(R2_THRESHOLD_GRID, function(thr) {
    removed_mask <- dat_qc$fooof == "FOOOFed" & !is.na(dat_qc$r_squared) & dat_qc$r_squared < thr
    kept <- dat_qc[!removed_mask, ]
    trial_key <- paste(kept$subjectID, kept$Condition, kept$Trial, sep = "_")
    tibble(
      task = task,
      r2_threshold = thr,
      rows_total = nrow(dat_qc),
      rows_retained = nrow(kept),
      rows_removed = nrow(dat_qc) - nrow(kept),
      rows_removed_pct = 100 * (nrow(dat_qc) - nrow(kept)) / max(nrow(dat_qc), 1),
      fooof_rows_total = sum(dat_qc$fooof == "FOOOFed", na.rm = TRUE),
      fooof_rows_removed = sum(removed_mask, na.rm = TRUE),
      fooof_rows_removed_pct = 100 * sum(removed_mask, na.rm = TRUE) / max(sum(dat_qc$fooof == "FOOOFed", na.rm = TRUE), 1),
      subjects_retained = n_distinct(kept$subjectID),
      trials_retained = n_distinct(trial_key),
      universes_retained = n_distinct(kept$universe_id)
    )
  }))

  summary_path <- file.path(csv_dir, paste0("multiverse_", task, "_fooof_r2_threshold_sensitivity.csv"))
  write.csv(threshold_summary, summary_path, row.names = FALSE)
  message("Saved threshold sensitivity summary: ", summary_path)

  n_before <- nrow(dat_qc)
  n_fooof_total <- sum(dat_qc$fooof == "FOOOFed", na.rm = TRUE)
  n_fooof_low_r2 <- sum(dat_qc$fooof == "FOOOFed" & !is.na(dat_qc$r_squared) & dat_qc$r_squared < R2_THRESHOLD, na.rm = TRUE)
  n_fooof_no_r2 <- sum(dat_qc$fooof == "FOOOFed" & is.na(dat_qc$r_squared), na.rm = TRUE)

  dat <- dat_qc %>%
    filter(!(fooof == "FOOOFed" & !is.na(r_squared) & r_squared < R2_THRESHOLD))

  n_after <- nrow(dat)
  message(sprintf("R\u00b2 filter (threshold = %.2f):", R2_THRESHOLD))
  message(sprintf("  FOOOFed rows: %d total, %d with R\u00b2 < %.2f removed (%.1f%%)",
                  n_fooof_total, n_fooof_low_r2, R2_THRESHOLD,
                  100 * n_fooof_low_r2 / max(n_fooof_total, 1)))
  if (n_fooof_no_r2 > 0) {
    message(sprintf("  %d FOOOFed rows had no matching R\u00b2 (kept).", n_fooof_no_r2))
  }
  message(sprintf("  Total rows: %d -> %d (%.1f%% removed)",
                  n_before, n_after, 100 * (n_before - n_after) / max(n_before, 1)))

  dat <- dat %>% select(-r_squared, -subjectID_num)

  # ========== STANDARD FILTERS (same as original) ==========
  if ("gaze_measure" %in% names(dat)) {
    dat <- dat %>% filter(gaze_measure != "gaze_density" | is.na(gaze_measure))
  }
  dat <- dat %>% filter(!electrodes %in% c("all", "parietal"))
  dat <- dat %>% filter(baseline_eeg != "dB", baseline_gaze != "dB")
  message(sprintf("After standard filters: %d rows, %d universes.", nrow(dat), n_distinct(dat$universe_id)))

  # ========== CONDITION SETUP ==========
  cond_levels  <- sort(unique(dat$Condition))
  if (task == "sternberg") {
    cond_labels <- setNames(paste0("Set size ", cond_levels), as.character(cond_levels))
  } else {
    cond_labels <- setNames(paste0(cond_levels, "-back"), as.character(cond_levels))
  }
  highest_cond <- max(as.numeric(as.character(cond_levels)))
  highest_alpha_term <- paste0("Condition", highest_cond)
  highest_int_term   <- paste0("gaze_value:Condition", highest_cond)
  highest_gaze_term  <- paste0("Condition", highest_cond)

  # ======================================================================
  #                         ANALYSIS
  # ======================================================================

  # ========== MAIN MULTIVERSE (7 dimensions) ==========
  message("Setting up main multiverse (7 dimensions)...")

  M <- multiverse()

  inside(M, {
    .elec   <- branch(electrodes,    "posterior", "occipital")
    .fooof  <- branch(fooof,         "FOOOFed", "nonFOOOFed")
    .lat    <- branch(latency_ms,    "0_500ms", "0_1000ms", "0_2000ms", "1000_2000ms")
    .alpha  <- branch(alpha_type,    "canonical", "IAF")
    .gaze   <- branch(gaze_measure,  "scan_path_length", "gaze_velocity", "microsaccades", "BCEA", "gaze_deviation")
    .bleeg  <- branch(baseline_eeg,  "raw", "pct_change")
    .blgaze <- branch(baseline_gaze, "raw", "pct_change")

    df <- dat %>%
      filter(electrodes == .elec, fooof == .fooof, latency_ms == .lat,
             alpha_type == .alpha, gaze_measure == .gaze,
             baseline_eeg == .bleeg, baseline_gaze == .blgaze) %>%
      filter(complete.cases(alpha, gaze_value, Condition, subjectID))

    df$gaze_value <- robust_z(df$gaze_value)
    df$alpha      <- robust_z(df$alpha)
    valid <- nrow(df) >= 10 && !any(is.nan(df$gaze_value)) && !any(is.nan(df$alpha))

    tid_int <- if (valid) {
      fit <- tryCatch(
        lmer(alpha ~ gaze_value * Condition + (1 | subjectID), data = df,
             control = lmerControl(optimizer = "bobyqa")),
        error = function(e) NULL)
      if (!is.null(fit)) broom.mixed::tidy(fit, conf.int = TRUE) else tibble()
    } else tibble()

    tid_cond <- if (valid) {
      bind_rows(lapply(as.character(cond_levels), function(cl) {
        dc <- df[df$Condition == cl, ]
        if (nrow(dc) < 5) return(tibble())
        fit_c <- tryCatch(
          lmer(alpha ~ gaze_value + (1 | subjectID), data = dc,
               control = lmerControl(optimizer = "bobyqa")),
          error = function(e) NULL)
        if (is.null(fit_c)) return(tibble())
        tid_c <- broom.mixed::tidy(fit_c, conf.int = TRUE) %>% filter(term == "gaze_value")
        tid_c$cond_label <- cond_labels[cl]
        tid_c
      }))
    } else tibble()
  })

  message("Executing main multiverse...")
  execute_multiverse(M)

  message("Extracting results from main multiverse...")
  M_expanded <- expand(M)

  M_int <- M_expanded %>%
    mutate(tid = map(.results, safe_extract, "tid_int")) %>%
    filter(!map_lgl(tid, is.null)) %>%
    select(.universe, all_of(branch_cols), tid) %>%
    unnest(tid) %>%
    filter(grepl("gaze_value", term))

  M_cond <- M_expanded %>%
    mutate(tid = map(.results, safe_extract, "tid_cond")) %>%
    filter(!map_lgl(tid, is.null)) %>%
    select(.universe, all_of(branch_cols), tid) %>%
    unnest(tid)

  if (nrow(M_cond) == 0L) { message("No successful LMM fits — skipping task."); next }

  se_thresh <- quantile(M_cond$std.error, 0.95, na.rm = TRUE)
  bad_ids   <- M_cond %>% filter(std.error > se_thresh) %>% pull(.universe) %>% unique()
  M_cond    <- M_cond %>% filter(!.universe %in% bad_ids)
  M_int     <- M_int  %>% filter(!.universe %in% bad_ids)
  message(sprintf("Dropped %d unstable universes (SE > %.4f). %d remain.",
                  length(bad_ids), se_thresh, n_distinct(M_cond$.universe)))

  M_cond <- add_sig(M_cond)
  M_int  <- add_sig(M_int)

  M_interaction <- M_int %>% filter(term == highest_int_term)

  # ========== EEG-ONLY MULTIVERSE (5 dimensions) ==========
  message("Setting up EEG-only multiverse (5 dimensions)...")

  dat_eeg <- dat %>%
    filter(!electrodes %in% c("all", "parietal"), baseline_eeg != "dB") %>%
    select(subjectID, Condition, alpha, electrodes, fooof, latency_ms, alpha_type, baseline_eeg) %>%
    distinct()

  M_eeg <- multiverse()

  inside(M_eeg, {
    .elec  <- branch(electrodes,   "posterior", "occipital")
    .fooof <- branch(fooof,        "FOOOFed", "nonFOOOFed")
    .lat   <- branch(latency_ms,   "0_500ms", "0_1000ms", "0_2000ms", "1000_2000ms")
    .alpha <- branch(alpha_type,   "canonical", "IAF")
    .bleeg <- branch(baseline_eeg, "raw", "pct_change")

    de <- dat_eeg %>%
      filter(electrodes == .elec, fooof == .fooof, latency_ms == .lat,
             alpha_type == .alpha, baseline_eeg == .bleeg) %>%
      filter(complete.cases(alpha, Condition, subjectID))

    de$alpha <- robust_z(de$alpha)
    valid <- nrow(de) >= 10 && !any(is.nan(de$alpha))

    tid_ca <- if (valid) {
      fit <- tryCatch(
        lmer(alpha ~ Condition + (1 | subjectID), data = de,
             control = lmerControl(optimizer = "bobyqa")),
        error = function(e) NULL)
      if (!is.null(fit)) {
        broom.mixed::tidy(fit, conf.int = TRUE) %>% filter(term == highest_alpha_term)
      } else tibble()
    } else tibble()
  })

  message("Executing EEG-only multiverse...")
  execute_multiverse(M_eeg)

  M_eeg_expanded <- expand(M_eeg)
  M_ca <- M_eeg_expanded %>%
    mutate(tid = map(.results, safe_extract, "tid_ca")) %>%
    filter(!map_lgl(tid, is.null)) %>%
    select(.universe, all_of(eeg_branch_cols), tid) %>%
    unnest(tid)

  if (nrow(M_ca) > 0L) M_ca <- add_sig(M_ca)

  # ========== GAZE-ONLY MULTIVERSE (3 dimensions) ==========
  message("Setting up gaze-only multiverse (3 dimensions)...")

  dat_gaze <- dat %>%
    select(subjectID, Condition, gaze_value, latency_ms, gaze_measure, baseline_gaze) %>%
    distinct()

  M_gaze <- multiverse()

  inside(M_gaze, {
    .lat    <- branch(latency_ms,    "0_500ms", "0_1000ms", "0_2000ms", "1000_2000ms")
    .gaze   <- branch(gaze_measure,  "scan_path_length", "gaze_velocity", "microsaccades", "BCEA", "gaze_deviation")
    .blgaze <- branch(baseline_gaze, "raw", "pct_change")

    dg <- dat_gaze %>%
      filter(latency_ms == .lat, gaze_measure == .gaze, baseline_gaze == .blgaze) %>%
      filter(complete.cases(gaze_value, Condition, subjectID))

    dg$gaze_value <- robust_z(dg$gaze_value)
    valid <- nrow(dg) >= 10 && !any(is.nan(dg$gaze_value))

    tid_cg <- if (valid) {
      fit <- tryCatch(
        lmer(gaze_value ~ Condition + (1 | subjectID), data = dg,
             control = lmerControl(optimizer = "bobyqa")),
        error = function(e) NULL)
      if (!is.null(fit)) {
        broom.mixed::tidy(fit, conf.int = TRUE) %>% filter(term == highest_gaze_term)
      } else tibble()
    } else tibble()
  })

  message("Executing gaze-only multiverse...")
  execute_multiverse(M_gaze)

  M_gaze_expanded <- expand(M_gaze)
  M_cg <- M_gaze_expanded %>%
    mutate(tid = map(.results, safe_extract, "tid_cg")) %>%
    filter(!map_lgl(tid, is.null)) %>%
    select(.universe, all_of(gaze_branch_cols), tid) %>%
    unnest(tid)

  if (nrow(M_cg) > 0L) M_cg <- add_sig(M_cg)

  # ========== APERIODIC MULTIVERSE ==========
  has_aperiodic <- "aperiodic_offset" %in% names(dat) && "aperiodic_exponent" %in% names(dat)
  M_ap_gaze_results <- tibble()
  M_ap_cond_results <- tibble()

  if (has_aperiodic) {
    message("Setting up aperiodic multiverse...")

    dat_ap <- dat %>%
      filter(fooof == "FOOOFed") %>%
      filter(complete.cases(aperiodic_offset, aperiodic_exponent))

    # --- Aperiodic ~ gaze (4D) ---
    dat_ap_gaze <- dat_ap %>%
      select(subjectID, Condition, aperiodic_offset, aperiodic_exponent,
             gaze_value, electrodes, latency_ms, gaze_measure, baseline_gaze) %>%
      distinct()

    M_ap_gaze <- multiverse()

    inside(M_ap_gaze, {
      .elec   <- branch(electrodes,    "posterior", "occipital")
      .lat    <- branch(latency_ms,    "0_500ms", "0_1000ms", "0_2000ms", "1000_2000ms")
      .gaze   <- branch(gaze_measure,  "scan_path_length", "gaze_velocity", "microsaccades", "BCEA", "gaze_deviation")
      .blgaze <- branch(baseline_gaze, "raw", "pct_change")

      dap <- dat_ap_gaze %>%
        filter(electrodes == .elec, latency_ms == .lat,
               gaze_measure == .gaze, baseline_gaze == .blgaze) %>%
        filter(complete.cases(aperiodic_exponent, aperiodic_offset, gaze_value, Condition, subjectID))

      dap$gaze_value <- robust_z(dap$gaze_value)
      dap$aperiodic_exponent <- robust_z(dap$aperiodic_exponent)
      dap$aperiodic_offset <- robust_z(dap$aperiodic_offset)
      valid <- nrow(dap) >= 10 && !any(is.nan(dap$gaze_value)) &&
               !any(is.nan(dap$aperiodic_exponent)) && !any(is.nan(dap$aperiodic_offset))

      tid_exp_gaze <- if (valid) {
        fit <- tryCatch(
          lmer(aperiodic_exponent ~ gaze_value + (1 | subjectID), data = dap,
               control = lmerControl(optimizer = "bobyqa")),
          error = function(e) NULL)
        if (!is.null(fit)) {
          broom.mixed::tidy(fit, conf.int = TRUE) %>% filter(term == "gaze_value") %>%
            mutate(aperiodic_measure = "Exponent")
        } else tibble()
      } else tibble()

      tid_off_gaze <- if (valid) {
        fit <- tryCatch(
          lmer(aperiodic_offset ~ gaze_value + (1 | subjectID), data = dap,
               control = lmerControl(optimizer = "bobyqa")),
          error = function(e) NULL)
        if (!is.null(fit)) {
          broom.mixed::tidy(fit, conf.int = TRUE) %>% filter(term == "gaze_value") %>%
            mutate(aperiodic_measure = "Offset")
        } else tibble()
      } else tibble()
    })

    message("Executing aperiodic ~ gaze multiverse...")
    execute_multiverse(M_ap_gaze)

    M_ap_gaze_exp <- expand(M_ap_gaze)
    M_ap_gaze_results <- bind_rows(
      M_ap_gaze_exp %>%
        mutate(tid = map(.results, safe_extract, "tid_exp_gaze")) %>%
        filter(!map_lgl(tid, is.null)) %>%
        select(.universe, all_of(ap_gaze_branch_cols), tid) %>%
        unnest(tid),
      M_ap_gaze_exp %>%
        mutate(tid = map(.results, safe_extract, "tid_off_gaze")) %>%
        filter(!map_lgl(tid, is.null)) %>%
        select(.universe, all_of(ap_gaze_branch_cols), tid) %>%
        unnest(tid)
    )
    if (nrow(M_ap_gaze_results) > 0) M_ap_gaze_results <- add_sig(M_ap_gaze_results)

    # --- Aperiodic ~ condition (2D) ---
    dat_ap_eeg <- dat_ap %>%
      select(subjectID, Condition, aperiodic_offset, aperiodic_exponent, electrodes, latency_ms) %>%
      distinct()

    M_ap_cond <- multiverse()

    inside(M_ap_cond, {
      .elec <- branch(electrodes, "posterior", "occipital")
      .lat  <- branch(latency_ms, "0_500ms", "0_1000ms", "0_2000ms", "1000_2000ms")

      dap <- dat_ap_eeg %>%
        filter(electrodes == .elec, latency_ms == .lat) %>%
        filter(complete.cases(aperiodic_exponent, aperiodic_offset, Condition, subjectID))

      dap$aperiodic_exponent <- robust_z(dap$aperiodic_exponent)
      dap$aperiodic_offset <- robust_z(dap$aperiodic_offset)
      valid <- nrow(dap) >= 10 && !any(is.nan(dap$aperiodic_exponent)) && !any(is.nan(dap$aperiodic_offset))

      tid_exp_cond <- if (valid) {
        fit <- tryCatch(
          lmer(aperiodic_exponent ~ Condition + (1 | subjectID), data = dap,
               control = lmerControl(optimizer = "bobyqa")),
          error = function(e) NULL)
        if (!is.null(fit)) {
          broom.mixed::tidy(fit, conf.int = TRUE) %>%
            filter(term == highest_alpha_term) %>%
            mutate(aperiodic_measure = "Exponent")
        } else tibble()
      } else tibble()

      tid_off_cond <- if (valid) {
        fit <- tryCatch(
          lmer(aperiodic_offset ~ Condition + (1 | subjectID), data = dap,
               control = lmerControl(optimizer = "bobyqa")),
          error = function(e) NULL)
        if (!is.null(fit)) {
          broom.mixed::tidy(fit, conf.int = TRUE) %>%
            filter(term == highest_alpha_term) %>%
            mutate(aperiodic_measure = "Offset")
        } else tibble()
      } else tibble()
    })

    message("Executing aperiodic ~ condition multiverse...")
    execute_multiverse(M_ap_cond)

    M_ap_cond_exp <- expand(M_ap_cond)
    M_ap_cond_results <- bind_rows(
      M_ap_cond_exp %>%
        mutate(tid = map(.results, safe_extract, "tid_exp_cond")) %>%
        filter(!map_lgl(tid, is.null)) %>%
        select(.universe, all_of(ap_cond_branch_cols), tid) %>%
        unnest(tid),
      M_ap_cond_exp %>%
        mutate(tid = map(.results, safe_extract, "tid_off_cond")) %>%
        filter(!map_lgl(tid, is.null)) %>%
        select(.universe, all_of(ap_cond_branch_cols), tid) %>%
        unnest(tid)
    )
    if (nrow(M_ap_cond_results) > 0) M_ap_cond_results <- add_sig(M_ap_cond_results)
  } else {
    message("Skipping aperiodic multiverse: columns not found in CSV.")
  }

  message(sprintf("=== %s analysis complete ===", toupper(task)))

  # ======================================================================
  #                         VISUALIZATION
  # ======================================================================

  M_cond_viz <- M_cond %>% filter(gaze_measure %in% gaze_order)
  M_cond_viz$condition <- factor(M_cond_viz$condition, levels = sig_levels)

  cond_labels_in_data <- unique(M_cond_viz$cond_label)
  cond_nums <- as.numeric(gsub("[^0-9]", "", cond_labels_in_data))
  highest_label <- cond_labels_in_data[which.max(cond_nums)]
  message(sprintf("Highest condition label: %s", highest_label))

  M_high <- M_cond_viz %>% filter(cond_label == highest_label)

  # ---------- FIGURE 1: alpha ~ gaze (sorted by estimate) ----------
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
  ggsave(file.path(storage_plot, paste0("AOC_multiverse_", task, "_estimate.png")),
         plot = p_combined, width = 14, height = 12, dpi = 600, bg = "white")
  message(sprintf("Saved: AOC_multiverse_%s_estimate.png", task))

  # ---------- FIGURE 2: alpha ~ gaze (grouped by processing hierarchy) ----------
  df_grouped <- M_high %>%
    mutate(
      .elec_ord = match(electrodes, elec_order),
      .lat_ord  = match(latency_ms, lat_order),
      .gaze_ord = match(gaze_measure, gaze_order)
    ) %>%
    arrange(fooof, .lat_ord, .elec_ord, baseline_eeg, alpha_type, .gaze_ord, baseline_gaze) %>%
    select(-.elec_ord, -.lat_ord, -.gaze_ord)
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
  ggsave(file.path(storage_plot, paste0("AOC_multiverse_", task, "_grouped.png")),
         plot = p_grp_combined, width = 14, height = 12, dpi = 600, bg = "white")
  message(sprintf("Saved: AOC_multiverse_%s_grouped.png", task))

  # ---------- FIGURE 3: alpha ~ condition (EEG-only) ----------
  if (nrow(M_ca) > 0L) {
    M_ca$condition <- factor(M_ca$condition, levels = sig_levels)
    M_ca_viz <- M_ca %>%
      mutate(.lat_ord = match(latency_ms, lat_order),
             .elec_ord = match(electrodes, elec_order)) %>%
      arrange(.lat_ord, .elec_ord, fooof, baseline_eeg, alpha_type) %>%
      select(-.lat_ord, -.elec_ord) %>%
      mutate(ordered_universe = row_number())

    ymax_ca <- max(abs(c(M_ca_viz$conf.low, M_ca_viz$conf.high)), na.rm = TRUE) * 1.05
    ylim_ca <- c(-ymax_ca, ymax_ca)

    p_ca_curve <- ggplot(M_ca_viz, aes(x = ordered_universe, y = estimate, color = condition)) +
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

    df_ca_specs <- M_ca_viz %>%
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
    ggsave(file.path(storage_plot, paste0("AOC_multiverse_", task, "_condition_alpha.png")),
           plot = p_ca_combined, width = 14, height = 10, dpi = 600, bg = "white")
    message(sprintf("Saved: AOC_multiverse_%s_condition_alpha.png", task))
  } else {
    message("Skipping Figure 3: no condition -> alpha results.")
  }

  # ---------- FIGURE 4: alpha ~ gaze x condition (interaction) ----------
  if (nrow(M_interaction) > 0) {
    M_interaction$condition <- factor(M_interaction$condition, levels = sig_levels)
    M_int_viz <- M_interaction %>% filter(gaze_measure %in% gaze_order)

    M_int_viz <- M_int_viz %>%
      mutate(.lat_ord = match(latency_ms, lat_order),
             .elec_ord = match(electrodes, elec_order),
             .gaze_ord = match(gaze_measure, gaze_order)) %>%
      arrange(.lat_ord, .elec_ord, fooof, baseline_eeg, alpha_type, .gaze_ord, baseline_gaze) %>%
      select(-.lat_ord, -.elec_ord, -.gaze_ord) %>%
      mutate(ordered_universe = row_number())

    ymax_int <- max(abs(c(M_int_viz$conf.low, M_int_viz$conf.high)), na.rm = TRUE) * 1.05
    ylim_int <- c(-ymax_int, ymax_int)

    p_int_curve <- ggplot(M_int_viz, aes(x = ordered_universe, y = estimate, color = condition)) +
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

    df_int_specs <- M_int_viz %>%
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
    ggsave(file.path(storage_plot, paste0("AOC_multiverse_", task, "_interaction.png")),
           plot = p_int_combined, width = 14, height = 12, dpi = 600, bg = "white")
    message(sprintf("Saved: AOC_multiverse_%s_interaction.png", task))
  } else {
    message("Skipping Figure 4: no interaction results.")
  }

  # ---------- FIGURE 5: gaze ~ condition (gaze-only) ----------
  if (nrow(M_cg) > 0L) {
    M_cg$condition <- factor(M_cg$condition, levels = sig_levels)
    M_cg_viz <- M_cg %>% filter(gaze_measure %in% gaze_order)

    M_cg_viz <- M_cg_viz %>%
      mutate(.lat_ord = match(latency_ms, lat_order),
             .gaze_ord = match(gaze_measure, gaze_order)) %>%
      arrange(.lat_ord, .gaze_ord, baseline_gaze) %>%
      select(-.lat_ord, -.gaze_ord) %>%
      mutate(ordered_universe = row_number())

    ymax_cg <- max(abs(c(M_cg_viz$conf.low, M_cg_viz$conf.high)), na.rm = TRUE) * 1.05
    ylim_cg <- c(-ymax_cg, ymax_cg)

    p_cg_curve <- ggplot(M_cg_viz, aes(x = ordered_universe, y = estimate, color = condition)) +
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

    df_cg_specs <- M_cg_viz %>%
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
    ggsave(file.path(storage_plot, paste0("AOC_multiverse_", task, "_condition_gaze.png")),
           plot = p_cg_combined, width = 14, height = 8, dpi = 600, bg = "white")
    message(sprintf("Saved: AOC_multiverse_%s_condition_gaze.png", task))
  } else {
    message("Skipping Figure 5: no gaze condition results.")
  }

  # ---------- FIGURE 6: Aperiodic ~ gaze (specification curves) ----------
  if (nrow(M_ap_gaze_results) > 0) {
    M_ap_spec <- M_ap_gaze_results
    M_ap_spec$condition <- factor(M_ap_spec$condition, levels = sig_levels)
    M_ap_spec <- M_ap_spec %>% filter(gaze_measure %in% gaze_order)

    for (ap_measure in c("Exponent", "Offset")) {
      M_ap <- M_ap_spec %>%
        filter(aperiodic_measure == ap_measure) %>%
        mutate(
          .lat_ord  = match(latency_ms, lat_order),
          .elec_ord = match(electrodes, elec_order),
          .gaze_ord = match(gaze_measure, gaze_order)
        ) %>%
        arrange(.lat_ord, .elec_ord, .gaze_ord, baseline_gaze) %>%
        select(-.lat_ord, -.elec_ord, -.gaze_ord) %>%
        mutate(ordered_universe = row_number())

      if (nrow(M_ap) == 0) next

      ymax_ap <- max(abs(c(M_ap$conf.low, M_ap$conf.high)), na.rm = TRUE) * 1.05
      ylim_ap <- c(-ymax_ap, ymax_ap)

      title_text <- paste0("Aperiodic ", tolower(ap_measure), " ~ gaze")
      subtitle_text <- paste0("aperiodic_", tolower(ap_measure),
                              " ~ gaze_value + (1|subjectID)")

      p_ap_curve <- ggplot(M_ap, aes(x = ordered_universe, y = estimate,
                                      color = condition)) +
        geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                      width = 0.3, alpha = 0.7) +
        geom_point(size = 1.5, alpha = 0.8) +
        geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
        scale_color_manual(values = sig_colors, name = "Significance") +
        guides(alpha = "none") +
        labs(title = bquote(bold(.(title_text))),
             subtitle = subtitle_text,
             x = "Universe", y = expression(bold("Standardized " * beta))) +
        theme_minimal() + theme(legend.position = "none") + v_common_theme +
        coord_cartesian(ylim = ylim_ap)
      legend_ap_sc <- get_legend(p_ap_curve + theme(legend.position = "bottom"))

      df_ap_specs <- M_ap %>%
        select(ordered_universe, latency_ms, electrodes, gaze_measure,
               baseline_gaze, condition)
      df_ap_long <- make_aperiodic_panel_long(df_ap_specs)

      p_ap_panel <- ggplot(df_ap_long, aes(x = ordered_universe, y = value,
                                            fill = condition)) +
        geom_tile() +
        scale_fill_manual(values = sig_colors, name = "Significance") +
        facet_grid(Variable ~ ., scales = "free_y", space = "free") +
        theme_minimal() +
        theme(strip.text.y = element_text(angle = 0, size = 13, face = "bold",
                                          hjust = 0),
              legend.position = "bottom") +
        labs(x = "Universe", y = "Analysis Decision") + v_common_theme

      p_ap_sc_combined <- p_ap_curve / legend_ap_sc / p_ap_panel +
        plot_layout(heights = c(0.8, 0.1, 1.2))

      suffix <- tolower(ap_measure)
      fname <- paste0("AOC_multiverse_", task, "_aperiodic_", suffix, "_spec.png")
      ggsave(file.path(storage_plot, fname),
             plot = p_ap_sc_combined, width = 14, height = 10, dpi = 600, bg = "white")
      message(sprintf("Saved: %s", fname))
    }
  } else {
    message("Skipping Figure 6: no aperiodic gaze results.")
  }

  message(sprintf("=== %s VISUALIZATION complete ===\n", toupper(task)))
}

message("=== All tasks complete. Figures saved to: ", storage_plot, " ===")
