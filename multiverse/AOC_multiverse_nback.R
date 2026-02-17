# AOC Multiverse — N-Back: Specification Curve Analysis
# Uses multiverse R package (Sarma et al., 2021) for systematic decision branching.
# Reads multiverse_nback.csv, fits LMMs per universe, plots specification curves.
#
# Models:
#   Figure 1 (_estimate):        alpha ~ gaze_value + (1|subjectID) [per condition, sorted by estimate]
#   Figure 2 (_grouped):         alpha ~ gaze_value + (1|subjectID) [per condition, grouped by dimension]
#   Figure 3 (_condition_alpha): alpha ~ Condition + (1|subjectID)  [EEG-only universes]
#   Figure 4 (_interaction):     alpha ~ gaze_value * Condition + (1|subjectID) [interaction term]
#   Figure 5 (_condition_gaze):  gaze_value ~ Condition + (1|subjectID) [gaze-only universes]
#
# Figure output: '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/multiverse'

library(lme4)
library(lmerTest)
library(broom.mixed)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(cowplot)
library(multiverse)

# Paths
csv_dir <- Sys.getenv("AOC_MULTIVERSE_DIR", unset = "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features")
csv_path <- file.path(csv_dir, "multiverse_nback.csv")
if (!file.exists(csv_path)) stop("CSV not found: ", csv_path)
storage_plot <- Sys.getenv("AOC_MULTIVERSE_FIGURES", unset = "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/multiverse")
if (!dir.exists(storage_plot)) dir.create(storage_plot, recursive = TRUE)

# Theme
v_common_theme <- theme(
  plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
  axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
  axis.text.y = element_text(size = 14),
  axis.title.x = element_text(size = 14, face = "bold"),
  axis.title.y = element_text(size = 14, face = "bold"),
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 12, face = "bold"),
  plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
)

# ========== LABEL MAPPINGS ==========
group_labels <- c(
  "electrodes" = "Electrodes", "fooof" = "FOOOF", "latency_ms" = "Latency",
  "alpha_type" = "Alpha", "gaze_measure" = "Gaze Measure",
  "baseline_eeg" = "EEG Baseline", "baseline_gaze" = "Gaze Baseline"
)

rename_opts <- function(x) {
  recode(x,
    "posterior" = "Posterior", "occipital" = "Occipital",
    "FOOOFed" = "FOOOF", "nonFOOOFed" = "No FOOOF",
    "0_500ms" = "0\u2013500 ms", "0_1000ms" = "0\u20131000 ms",
    "0_2000ms" = "0\u20132000 ms", "1000_2000ms" = "1000\u20132000 ms",
    "canonical" = "Canonical", "IAF" = "IAF",
    "scan_path_length" = "Scan Path Length", "gaze_velocity" = "Gaze Velocity",
    "microsaccades" = "Microsaccades", "BCEA" = "BCEA",
    "raw" = "Raw", "dB" = "dB", "pct_change" = "% Change"
  )
}

value_levels <- c(
  "% Change", "Raw",
  "BCEA", "Microsaccades", "Gaze Velocity", "Scan Path Length",
  "IAF", "Canonical",
  "No FOOOF", "FOOOF",
  "Occipital", "Posterior",
  "1000\u20132000 ms", "0\u20132000 ms", "0\u20131000 ms", "0\u2013500 ms"
)

v_p2_group_order <- c("Latency", "Electrodes", "FOOOF",
                       "EEG Baseline", "Alpha", "Gaze Measure", "Gaze Baseline")

elec_order <- c("posterior", "occipital")
lat_order  <- c("0_500ms", "0_1000ms", "0_2000ms", "1000_2000ms")

# ========== LOAD & FILTER DATA ==========
dat <- read.csv(csv_path, stringsAsFactors = FALSE, na.strings = c("NA", "NaN", ""))
dat$subjectID <- as.factor(dat$subjectID)
dat$Condition <- as.factor(dat$Condition)

message(sprintf("Loaded: %d rows, %d universes.", nrow(dat), n_distinct(dat$universe_id)))

if ("gaze_measure" %in% names(dat)) {
  dat <- dat %>% filter(gaze_measure != "gaze_density" | is.na(gaze_measure))
}
dat <- dat %>% filter(!electrodes %in% c("all", "parietal"))
dat <- dat %>% filter(baseline_eeg != "dB", baseline_gaze != "dB")
message(sprintf("After filtering: %d rows, %d universes.", nrow(dat), n_distinct(dat$universe_id)))

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

sig_colors <- c("Positive" = "#33CC66", "Negative" = "#fe0000",
                "Non-significant" = "#d1d1d1")

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

safe_extract <- function(results, var_name) {
  r <- results[[var_name]]
  if (is.null(r) || !is.data.frame(r) || nrow(r) == 0) return(NULL)
  r
}

branch_cols <- c("electrodes", "fooof", "latency_ms", "alpha_type",
                 "gaze_measure", "baseline_eeg", "baseline_gaze")

# Condition setup
cond_levels <- sort(unique(dat$Condition))
cond_labels <- setNames(paste0(cond_levels, "-back"), as.character(cond_levels))
highest_cond <- max(as.numeric(as.character(cond_levels)))
highest_label <- cond_labels[as.character(highest_cond)]

# ========== MAIN MULTIVERSE (7 dimensions, 512 universes) ==========
message("Setting up main multiverse (7 dimensions)...")

M <- multiverse()

inside(M, {
  .elec   <- branch(electrodes,    "posterior", "occipital")
  .fooof  <- branch(fooof,         "FOOOFed", "nonFOOOFed")
  .lat    <- branch(latency_ms,    "0_500ms", "0_1000ms", "0_2000ms", "1000_2000ms")
  .alpha  <- branch(alpha_type,    "canonical", "IAF")
  .gaze   <- branch(gaze_measure,  "scan_path_length", "gaze_velocity", "microsaccades", "BCEA")
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

  # Interaction model: alpha ~ gaze_value * Condition + (1|subjectID)
  tid_int <- if (valid) {
    fit <- tryCatch(
      lmer(alpha ~ gaze_value * Condition + (1 | subjectID), data = df,
           control = lmerControl(optimizer = "bobyqa")),
      error = function(e) NULL)
    if (!is.null(fit)) broom.mixed::tidy(fit, conf.int = TRUE) else tibble()
  } else tibble()

  # Per-condition models: alpha ~ gaze_value + (1|subjectID)
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

message("Executing main multiverse (512 universes)...")
execute_multiverse(M)

# ========== EXTRACT RESULTS ==========
message("Extracting results from multiverse...")
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
bad_ids <- M_cond %>% filter(std.error > se_thresh) %>% pull(.universe) %>% unique()
M_cond <- M_cond %>% filter(!.universe %in% bad_ids)
M_int  <- M_int  %>% filter(!.universe %in% bad_ids)
message(sprintf("Dropped %d unstable universes (SE > %.4f). %d remain.",
                length(bad_ids), se_thresh, n_distinct(M_cond$.universe)))

M_cond <- add_sig(M_cond)
M_int  <- add_sig(M_int)

# Save CSVs
write.csv(M_int, file.path(csv_dir, "multiverse_nback_results.csv"), row.names = FALSE)
write.csv(M_cond, file.path(csv_dir, "multiverse_nback_conditions_results.csv"), row.names = FALSE)

# ========== HIGHEST-CONDITION DATA ==========
M_high <- M_cond %>% filter(cond_label == highest_label)
ord_high <- M_high %>% arrange(fooof, estimate) %>% pull(.universe)
ord_df <- data.frame(.universe = ord_high, ordered_universe = seq_along(ord_high))
M_high <- M_high %>% left_join(ord_df, by = ".universe")

# ========== FIGURE 1: SPECIFICATION CURVE (sorted by estimate) ==========
ymax_est <- max(abs(c(M_high$conf.low, M_high$conf.high)), na.rm = TRUE) * 1.05
ylim_est <- c(-ymax_est, ymax_est)

p_curve <- ggplot(M_high, aes(x = ordered_universe, y = estimate, color = condition)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.3, alpha = 0.7) +
  geom_point(size = 1.2, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  scale_color_manual(values = sig_colors, name = "Significance") +
  guides(alpha = "none") +
  labs(title = expression(bold(alpha ~ "~" ~ gaze)), x = "Universe",
       y = expression(bold("Standardized " * beta))) +
  theme_minimal() + theme(legend.position = "none") + v_common_theme +
  coord_cartesian(ylim = ylim_est)
legend_spec <- get_legend(p_curve + theme(legend.position = "bottom") + guides(alpha = "none"))

df_specs <- M_high %>%
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
ggsave(file.path(storage_plot, "AOC_multiverse_nback_estimate.png"),
       plot = p_combined, width = 14, height = 12, dpi = 600)
message("Saved AOC_multiverse_nback_estimate.png to ", storage_plot)

# ========== FIGURE 2: DIMENSION-GROUPED SPECIFICATION CURVE ==========
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
ggsave(file.path(storage_plot, "AOC_multiverse_nback_grouped.png"),
       plot = p_grp_combined, width = 14, height = 12, dpi = 600)
message("Saved AOC_multiverse_nback_grouped.png to ", storage_plot)

# ========== FIGURE 3: CONDITION EFFECT ON ALPHA (EEG-only universes) ==========
message("Setting up EEG-only multiverse (5 dimensions)...")

eeg_group_labels <- c(
  "electrodes" = "Electrodes", "fooof" = "FOOOF", "latency_ms" = "Latency",
  "alpha_type" = "Alpha", "baseline_eeg" = "EEG Baseline"
)
eeg_group_order <- c("Latency", "Electrodes", "FOOOF", "EEG Baseline", "Alpha")

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

eeg_branch_cols <- c("electrodes", "fooof", "latency_ms", "alpha_type", "baseline_eeg")
M_eeg_expanded <- expand(M_eeg)

M_ca <- M_eeg_expanded %>%
  mutate(tid = map(.results, safe_extract, "tid_ca")) %>%
  filter(!map_lgl(tid, is.null)) %>%
  select(.universe, all_of(eeg_branch_cols), tid) %>%
  unnest(tid)

if (nrow(M_ca) == 0L) {
  message("No successful condition \u2192 alpha fits.")
} else {
  M_ca <- add_sig(M_ca)
  M_ca <- M_ca %>%
    mutate(.lat_ord = match(latency_ms, lat_order),
           .elec_ord = match(electrodes, elec_order)) %>%
    arrange(.lat_ord, .elec_ord, fooof, baseline_eeg, alpha_type) %>%
    select(-.lat_ord, -.elec_ord) %>%
    mutate(ordered_universe = row_number())
  message(sprintf("Condition \u2192 alpha: %d EEG universes fitted.", nrow(M_ca)))

  ymax_ca <- max(abs(c(M_ca$conf.low, M_ca$conf.high)), na.rm = TRUE) * 1.05
  ylim_ca <- c(-ymax_ca, ymax_ca)

  p_ca_curve <- ggplot(M_ca, aes(x = ordered_universe, y = estimate, color = condition)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.3, alpha = 0.7) +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    scale_color_manual(values = sig_colors, name = "Significance") +
    guides(alpha = "none") +
    labs(title = expression(bold(alpha ~ "~" ~ condition)),
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
  ggsave(file.path(storage_plot, "AOC_multiverse_nback_condition_alpha.png"),
         plot = p_ca_combined, width = 14, height = 10, dpi = 600)
  message("Saved AOC_multiverse_nback_condition_alpha.png to ", storage_plot)

  write.csv(M_ca, file.path(csv_dir, "multiverse_nback_condition_results.csv"), row.names = FALSE)
}

# ========== FIGURE 4: INTERACTION (gaze x condition -> alpha) ==========
message("Building interaction specification curve...")

highest_int_term <- paste0("gaze_value:Condition", highest_cond)
M_interaction <- M_int %>% filter(term == highest_int_term)

if (nrow(M_interaction) > 0) {
  M_interaction <- add_sig(M_interaction)
  ord_int <- M_interaction %>% arrange(fooof, estimate) %>% pull(.universe)
  ord_int_df <- data.frame(.universe = ord_int, ordered_universe = seq_along(ord_int))
  M_interaction <- M_interaction %>% left_join(ord_int_df, by = ".universe")

  ymax_int <- max(abs(c(M_interaction$conf.low, M_interaction$conf.high)), na.rm = TRUE) * 1.05
  ylim_int <- c(-ymax_int, ymax_int)

  p_int_curve <- ggplot(M_interaction, aes(x = ordered_universe, y = estimate, color = condition)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.3, alpha = 0.7) +
    geom_point(size = 1.2, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    scale_color_manual(values = sig_colors, name = "Significance") +
    guides(alpha = "none") +
    labs(title = expression(bold(alpha ~ "~" ~ gaze ~ "\u00D7" ~ condition)),
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
  ggsave(file.path(storage_plot, "AOC_multiverse_nback_interaction.png"),
         plot = p_int_combined, width = 14, height = 12, dpi = 600)
  message("Saved AOC_multiverse_nback_interaction.png to ", storage_plot)

  write.csv(M_interaction, file.path(csv_dir, "multiverse_nback_interaction_results.csv"), row.names = FALSE)
} else {
  message("No interaction terms found for ", highest_int_term)
}

# ========== FIGURE 5: CONDITION EFFECT ON GAZE (gaze-only universes) ==========
message("Setting up gaze-only multiverse (3 dimensions)...")

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

dat_gaze <- dat %>%
  select(subjectID, Condition, gaze_value, latency_ms, gaze_measure, baseline_gaze) %>%
  distinct()

highest_gaze_term <- paste0("Condition", highest_cond)

M_gaze <- multiverse()

inside(M_gaze, {
  .lat    <- branch(latency_ms,    "0_500ms", "0_1000ms", "0_2000ms", "1000_2000ms")
  .gaze   <- branch(gaze_measure,  "scan_path_length", "gaze_velocity", "microsaccades", "BCEA")
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
M_gaze_expanded <- expand(M_gaze)

M_cg <- M_gaze_expanded %>%
  mutate(tid = map(.results, safe_extract, "tid_cg")) %>%
  filter(!map_lgl(tid, is.null)) %>%
  select(.universe, all_of(gaze_branch_cols), tid) %>%
  unnest(tid)

if (nrow(M_cg) == 0L) {
  message("No successful condition \u2192 gaze fits.")
} else {
  M_cg <- add_sig(M_cg)
  M_cg <- M_cg %>%
    mutate(.lat_ord = match(latency_ms, lat_order)) %>%
    arrange(.lat_ord, gaze_measure, baseline_gaze) %>%
    select(-.lat_ord) %>%
    mutate(ordered_universe = row_number())
  message(sprintf("Condition \u2192 gaze: %d gaze universes fitted.", nrow(M_cg)))

  ymax_cg <- max(abs(c(M_cg$conf.low, M_cg$conf.high)), na.rm = TRUE) * 1.05
  ylim_cg <- c(-ymax_cg, ymax_cg)

  p_cg_curve <- ggplot(M_cg, aes(x = ordered_universe, y = estimate, color = condition)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.3, alpha = 0.7) +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    scale_color_manual(values = sig_colors, name = "Significance") +
    guides(alpha = "none") +
    labs(title = expression(bold(gaze ~ "~" ~ condition)),
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
  ggsave(file.path(storage_plot, "AOC_multiverse_nback_condition_gaze.png"),
         plot = p_cg_combined, width = 14, height = 8, dpi = 600)
  message("Saved AOC_multiverse_nback_condition_gaze.png to ", storage_plot)

  write.csv(M_cg, file.path(csv_dir, "multiverse_nback_condition_gaze_results.csv"), row.names = FALSE)
}

message("=== N-back multiverse analysis complete ===")
