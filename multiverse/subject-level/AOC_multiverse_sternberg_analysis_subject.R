# AOC Multiverse — Sternberg: Subject-Level Analysis
# Loads pre-aggregated subject-level CSV (from AOC_multiverse_prep_subject.m),
# fits LMMs per universe. Directly comparable to existing raincloud analyses.
# Uses multiverse R package (Sarma et al., 2021).
#
# Models:
#   M (main, 7D):      alpha ~ gaze_value * Condition + (1|subjectID)
#                       alpha ~ gaze_value + (1|subjectID) [per condition]
#   M_eeg (EEG, 5D):   alpha ~ Condition + (1|subjectID)
#   M_gaze (gaze, 3D): gaze_value ~ Condition + (1|subjectID)
#
# Output CSVs (saved to csv_dir, with _subject_ prefix):
#   multiverse_sternberg_subject_results.csv
#   multiverse_sternberg_subject_conditions_results.csv
#   multiverse_sternberg_subject_condition_results.csv
#   multiverse_sternberg_subject_interaction_results.csv
#   multiverse_sternberg_subject_condition_gaze_results.csv

library(lme4)
library(lmerTest)
library(broom.mixed)
library(tidyverse)
library(multiverse)

# ========== PATHS ==========
csv_dir  <- Sys.getenv("AOC_MULTIVERSE_DIR",
                        unset = "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features")
csv_path <- file.path(csv_dir, "multiverse_sternberg_subject.csv")
if (!file.exists(csv_path)) stop("Subject-level CSV not found: ", csv_path,
                                  "\nRun AOC_multiverse_prep_subject.m first.")

# ========== LOAD & FILTER DATA ==========
dat <- read.csv(csv_path, stringsAsFactors = FALSE, na.strings = c("NA", "NaN", ""))
dat$subjectID <- as.factor(dat$subjectID)
dat$Condition <- as.factor(dat$Condition)

message(sprintf("Loaded: %d subject-level rows, %d subjects.",
                nrow(dat), n_distinct(dat$subjectID)))

if ("gaze_measure" %in% names(dat)) {
  dat <- dat %>% filter(gaze_measure != "gaze_density" | is.na(gaze_measure))
}
dat <- dat %>% filter(!electrodes %in% c("all", "parietal"))
dat <- dat %>% filter(baseline_eeg != "dB", baseline_gaze != "dB")
message(sprintf("After filtering: %d rows.", nrow(dat)))

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

branch_cols <- c("electrodes", "fooof", "latency_ms", "alpha_type",
                 "gaze_measure", "baseline_eeg", "baseline_gaze")

# Condition setup
cond_levels  <- sort(unique(dat$Condition))
cond_labels  <- setNames(paste0("Set size ", cond_levels), as.character(cond_levels))
highest_cond <- max(as.numeric(as.character(cond_levels)))

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
        lm(alpha ~ gaze_value, data = dc),
        error = function(e) NULL)
      if (is.null(fit_c)) return(tibble())
      tid_c <- broom::tidy(fit_c, conf.int = TRUE) %>% filter(term == "gaze_value")
      tid_c$cond_label <- cond_labels[cl]
      tid_c
    }))
  } else tibble()
})

message("Executing main multiverse...")
execute_multiverse(M)

# ========== EXTRACT RESULTS ==========
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

if (nrow(M_cond) == 0L) stop("No successful LMM fits.")

# Drop unstable universes (SE > 95th percentile)
se_thresh <- quantile(M_cond$std.error, 0.95, na.rm = TRUE)
bad_ids   <- M_cond %>% filter(std.error > se_thresh) %>% pull(.universe) %>% unique()
M_cond    <- M_cond %>% filter(!.universe %in% bad_ids)
M_int     <- M_int  %>% filter(!.universe %in% bad_ids)
message(sprintf("Dropped %d unstable universes (SE > %.4f). %d remain.",
                length(bad_ids), se_thresh, n_distinct(M_cond$.universe)))

M_cond <- add_sig(M_cond)
M_int  <- add_sig(M_int)

write.csv(M_int,  file.path(csv_dir, "multiverse_sternberg_subject_results.csv"), row.names = FALSE)
write.csv(M_cond, file.path(csv_dir, "multiverse_sternberg_subject_conditions_results.csv"), row.names = FALSE)
message("Saved: multiverse_sternberg_subject_results.csv, multiverse_sternberg_subject_conditions_results.csv")

# Extract interaction term
highest_int_term <- paste0("gaze_value:Condition", highest_cond)
M_interaction    <- M_int %>% filter(term == highest_int_term)
if (nrow(M_interaction) > 0) {
  write.csv(M_interaction, file.path(csv_dir, "multiverse_sternberg_subject_interaction_results.csv"), row.names = FALSE)
  message(sprintf("Saved: multiverse_sternberg_subject_interaction_results.csv (%d universes)", nrow(M_interaction)))
} else {
  message("WARNING: No interaction terms found for ", highest_int_term)
}

# ========== EEG-ONLY MULTIVERSE (5 dimensions) ==========
message("Setting up EEG-only multiverse (5 dimensions)...")

dat_eeg <- dat %>%
  filter(!electrodes %in% c("all", "parietal"), baseline_eeg != "dB") %>%
  select(subjectID, Condition, alpha, electrodes, fooof, latency_ms, alpha_type, baseline_eeg) %>%
  distinct()

highest_alpha_term <- paste0("Condition", highest_cond)

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

eeg_branch_cols  <- c("electrodes", "fooof", "latency_ms", "alpha_type", "baseline_eeg")
M_eeg_expanded   <- expand(M_eeg)

M_ca <- M_eeg_expanded %>%
  mutate(tid = map(.results, safe_extract, "tid_ca")) %>%
  filter(!map_lgl(tid, is.null)) %>%
  select(.universe, all_of(eeg_branch_cols), tid) %>%
  unnest(tid)

if (nrow(M_ca) > 0L) {
  M_ca <- add_sig(M_ca)
  write.csv(M_ca, file.path(csv_dir, "multiverse_sternberg_subject_condition_results.csv"), row.names = FALSE)
  message(sprintf("Saved: multiverse_sternberg_subject_condition_results.csv (%d universes)", nrow(M_ca)))
} else {
  message("WARNING: No successful condition \u2192 alpha fits.")
}

# ========== GAZE-ONLY MULTIVERSE (3 dimensions) ==========
message("Setting up gaze-only multiverse (3 dimensions)...")

dat_gaze <- dat %>%
  select(subjectID, Condition, gaze_value, latency_ms, gaze_measure, baseline_gaze) %>%
  distinct()

highest_gaze_term <- paste0("Condition", highest_cond)

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

gaze_branch_cols <- c("latency_ms", "gaze_measure", "baseline_gaze")
M_gaze_expanded  <- expand(M_gaze)

M_cg <- M_gaze_expanded %>%
  mutate(tid = map(.results, safe_extract, "tid_cg")) %>%
  filter(!map_lgl(tid, is.null)) %>%
  select(.universe, all_of(gaze_branch_cols), tid) %>%
  unnest(tid)

if (nrow(M_cg) > 0L) {
  M_cg <- add_sig(M_cg)
  write.csv(M_cg, file.path(csv_dir, "multiverse_sternberg_subject_condition_gaze_results.csv"), row.names = FALSE)
  message(sprintf("Saved: multiverse_sternberg_subject_condition_gaze_results.csv (%d universes)", nrow(M_cg)))
} else {
  message("WARNING: No successful condition \u2192 gaze fits.")
}

# ========== APERIODIC MULTIVERSE ==========
if ("aperiodic_offset" %in% names(dat) && "aperiodic_exponent" %in% names(dat)) {
  message("Setting up aperiodic multiverse (subject-level)...")

  dat_ap <- dat %>%
    filter(fooof == "FOOOFed") %>%
    filter(complete.cases(aperiodic_offset, aperiodic_exponent))

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

  ap_gaze_branch_cols <- c("electrodes", "latency_ms", "gaze_measure", "baseline_gaze")
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

  if (nrow(M_ap_gaze_results) > 0) {
    M_ap_gaze_results <- add_sig(M_ap_gaze_results)
    write.csv(M_ap_gaze_results, file.path(csv_dir, "multiverse_sternberg_subject_aperiodic_gaze_results.csv"), row.names = FALSE)
    message(sprintf("Saved: aperiodic_gaze_results.csv (%d rows)", nrow(M_ap_gaze_results)))
  }

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

  ap_cond_branch_cols <- c("electrodes", "latency_ms")
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

  if (nrow(M_ap_cond_results) > 0) {
    M_ap_cond_results <- add_sig(M_ap_cond_results)
    write.csv(M_ap_cond_results, file.path(csv_dir, "multiverse_sternberg_subject_aperiodic_condition_results.csv"), row.names = FALSE)
    message(sprintf("Saved: aperiodic_condition_results.csv (%d rows)", nrow(M_ap_cond_results)))
  }

} else {
  message("Skipping aperiodic multiverse: columns not found in CSV.")
}

message("=== Sternberg SUBJECT-LEVEL multiverse ANALYSIS complete ===")
message("Result CSVs saved to: ", csv_dir)
