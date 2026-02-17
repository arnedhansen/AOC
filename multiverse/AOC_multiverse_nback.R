# AOC Multiverse — N-Back: alpha ~ gaze * condition + (1|subjectID)
# Reads multiverse_nback.csv, fits LMM per universe, plots specification curve + panel.
# Default figures show the slope at the HIGHEST WM load (3-back).
# Figure output: '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/multiverse'

library(lme4)
library(lmerTest)
library(broom.mixed)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(cowplot)

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

# Factor levels for panel y-axis (bottom → top, so first level = bottom)
# Ordered by processing stage: Latency → Electrodes → FOOOF →
#   EEG Baseline → Alpha → Gaze Measure → Gaze Baseline
value_levels <- c(
  "% Change", "Raw",
  "BCEA", "Microsaccades", "Gaze Velocity", "Scan Path Length",
  "IAF", "Canonical",
  "No FOOOF", "FOOOF",
  "Occipital", "Posterior",
  "1000\u20132000 ms", "0\u20132000 ms", "0\u20131000 ms", "0\u2013500 ms"
)

# Dimension ordering for facet strips (top → bottom = processing order)
v_p2_group_order <- c("Latency", "Electrodes", "FOOOF",
                       "EEG Baseline", "Alpha", "Gaze Measure", "Gaze Baseline")

elec_order <- c("posterior", "occipital")
lat_order  <- c("0_500ms", "0_1000ms", "0_2000ms", "1000_2000ms")

# ========== LOAD & FILTER DATA ==========
dat <- read.csv(csv_path, stringsAsFactors = FALSE, na.strings = c("NA", "NaN", ""))
dat$subjectID <- as.factor(dat$subjectID)
dat$Condition <- as.factor(dat$Condition)

message(sprintf("Loaded: %d rows, %d universes.", nrow(dat), n_distinct(dat$universe_id)))

# Drop gaze_density if present (legacy CSVs) and exclude "all" electrodes
if ("gaze_measure" %in% names(dat)) {
  dat <- dat %>% filter(gaze_measure != "gaze_density" | is.na(gaze_measure))
}
dat <- dat %>% filter(!electrodes %in% c("all", "parietal"))
# Only show % Change baseline on plots (drop dB); raw + pct_change remain
dat <- dat %>% filter(baseline_eeg != "dB", baseline_gaze != "dB")
message(sprintf("After filtering: %d rows, %d universes.", nrow(dat), n_distinct(dat$universe_id)))

robust_z <- function(x, winsor = 3) {
  med <- median(x, na.rm = TRUE)
  ma  <- mad(x, na.rm = TRUE)
  if (is.na(ma) || ma == 0) return(rep(NaN, length(x)))
  z <- (x - med) / ma
  z[z >  winsor] <-  winsor
  z[z < -winsor] <- -winsor
  z
}

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

# ========== FIT MODELS ==========
cond_levels <- sort(unique(dat$Condition))
cond_labels <- setNames(paste0(cond_levels, "-back"), cond_levels)
highest_cond <- max(as.numeric(as.character(cond_levels)))
highest_label <- cond_labels[as.character(highest_cond)]

universes <- unique(dat$universe_id)
M_int_list  <- list()
M_cond_list <- list()

for (i in seq_along(universes)) {
  u <- universes[i]
  du <- dat[dat$universe_id == u, ]
  du <- du[complete.cases(du[, c("alpha", "gaze_value", "Condition", "subjectID")]), ]
  if (nrow(du) < 10) next
  du$gaze_value <- robust_z(du$gaze_value)
  du$alpha      <- robust_z(du$alpha)
  if (any(is.nan(du$gaze_value)) || any(is.nan(du$alpha))) next
  r1 <- du[1, ]

  fit_int <- tryCatch(
    lmer(alpha ~ gaze_value * Condition + (1 | subjectID), data = du,
         control = lmerControl(optimizer = "bobyqa")),
    error = function(e) NULL
  )
  if (!is.null(fit_int)) {
    tid <- broom.mixed::tidy(fit_int, conf.int = TRUE)
    tid$universe_id <- u
    tid$electrodes <- r1$electrodes; tid$fooof <- r1$fooof
    tid$latency_ms <- r1$latency_ms; tid$alpha_type <- r1$alpha_type
    tid$gaze_measure <- r1$gaze_measure; tid$baseline_eeg <- r1$baseline_eeg
    tid$baseline_gaze <- r1$baseline_gaze
    M_int_list[[length(M_int_list) + 1L]] <- tid
  }

  for (cl in cond_levels) {
    dc <- du[du$Condition == cl, ]
    if (nrow(dc) < 5) next
    fit_c <- tryCatch(
      lmer(alpha ~ gaze_value + (1 | subjectID), data = dc,
           control = lmerControl(optimizer = "bobyqa")),
      error = function(e) NULL
    )
    if (is.null(fit_c)) next
    tid_c <- broom.mixed::tidy(fit_c, conf.int = TRUE) %>% filter(term == "gaze_value")
    tid_c$universe_id <- u
    tid_c$cond_label <- cond_labels[as.character(cl)]
    tid_c$electrodes <- r1$electrodes; tid_c$fooof <- r1$fooof
    tid_c$latency_ms <- r1$latency_ms; tid_c$alpha_type <- r1$alpha_type
    tid_c$gaze_measure <- r1$gaze_measure; tid_c$baseline_eeg <- r1$baseline_eeg
    tid_c$baseline_gaze <- r1$baseline_gaze
    M_cond_list[[length(M_cond_list) + 1L]] <- tid_c
  }
}

M_int <- bind_rows(M_int_list)
terms_keep <- c("gaze_value", grep("gaze_value:Condition", M_int$term, value = TRUE))
M_int <- M_int %>% filter(term %in% unique(terms_keep))

M_cond <- bind_rows(M_cond_list)
if (nrow(M_cond) == 0L) stop("No successful LMM fits.")

se_thresh <- quantile(M_cond$std.error, 0.95, na.rm = TRUE)
bad_ids <- M_cond %>% filter(std.error > se_thresh) %>% pull(universe_id) %>% unique()
M_cond <- M_cond %>% filter(!universe_id %in% bad_ids)
M_int  <- M_int  %>% filter(!universe_id %in% bad_ids)
message(sprintf("Dropped %d unstable universes (SE > %.4f). %d remain.",
                length(bad_ids), se_thresh, n_distinct(M_cond$universe_id)))

add_sig <- function(df) {
  df %>% mutate(
    condition = factor(case_when(
      p.value < 0.05 & estimate > 0 ~ "Positive",
      p.value < 0.05 & estimate < 0 ~ "Negative",
      TRUE ~ "Non-significant"
    ), levels = c("Positive", "Negative", "Non-significant"))
  )
}
M_cond <- add_sig(M_cond)
M_int  <- add_sig(M_int)

sig_colors <- c("Positive" = "#33CC66", "Negative" = "#fe0000",
                "Non-significant" = "#d1d1d1")

# ========== HIGHEST-CONDITION DATA ==========
M_high <- M_cond %>% filter(cond_label == highest_label)
ord_high <- M_high %>% arrange(fooof, estimate) %>% pull(universe_id)
ord_df <- data.frame(universe_id = ord_high, ordered_universe = seq_along(ord_high))
M_high <- M_high %>% left_join(ord_df, by = "universe_id")

# ========== FIGURE 1: SPECIFICATION CURVE (highest condition, sorted by estimate) ==========
ymax_est <- max(abs(c(M_high$conf.low, M_high$conf.high)), na.rm = TRUE) * 1.05
ylim_est <- c(-ymax_est, ymax_est)

p_curve <- ggplot(M_high, aes(x = ordered_universe, y = estimate, color = condition)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.3, alpha = 0.7) +
  geom_point(size = 1.2, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  scale_color_manual(values = sig_colors, name = "Significance") +
  guides(alpha = "none") +
  labs(title = expression(bold(alpha ~ "~" ~ gaze)), x = "Universe", y = "Estimate") +
  theme_minimal() + theme(legend.position = "none") + v_common_theme +
  coord_cartesian(ylim = ylim_est)
legend_spec <- get_legend(p_curve + theme(legend.position = "bottom") + guides(alpha = "none"))

df_specs <- M_high %>%
  select(ordered_universe, universe_id, electrodes, fooof, latency_ms, alpha_type,
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
ggsave(file.path(storage_plot, "AOC_multiverse_nback_estimate.png"), plot = p_combined, width = 14, height = 12, dpi = 600)
message("Saved AOC_multiverse_nback_estimate.png to ", storage_plot)

write.csv(M_int, file.path(csv_dir, "multiverse_nback_results.csv"), row.names = FALSE)

write.csv(M_cond, file.path(csv_dir, "multiverse_nback_conditions_results.csv"), row.names = FALSE)

# ========== FIGURE 2: DIMENSION-GROUPED SPECIFICATION CURVE (highest condition) ==========
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
       x = "Universe", y = "Estimate") +
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

# ========== FIGURE 4: CONDITION EFFECT ON ALPHA (EEG-only universes) ==========
# Model: alpha ~ Condition_numeric + (1|subjectID)
# Only EEG dimensions matter (gaze dimensions don't affect alpha), so we deduplicate.
message("Fitting condition → alpha models for EEG-only universes...")

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
dat_eeg$eeg_uid <- interaction(dat_eeg$electrodes, dat_eeg$fooof, dat_eeg$latency_ms,
                                dat_eeg$alpha_type, dat_eeg$baseline_eeg, drop = TRUE)
dat_eeg$Condition_num <- as.numeric(as.character(dat_eeg$Condition))

eeg_universes <- unique(dat_eeg$eeg_uid)
M_cond_alpha_list <- list()

for (i in seq_along(eeg_universes)) {
  uid <- eeg_universes[i]
  de <- dat_eeg[dat_eeg$eeg_uid == uid, ]
  de <- de[complete.cases(de[, c("alpha", "Condition_num", "subjectID")]), ]
  if (nrow(de) < 10) next
  de$alpha <- robust_z(de$alpha)
  if (any(is.nan(de$alpha))) next
  r1 <- de[1, ]

  fit_ca <- tryCatch(
    lmer(alpha ~ Condition_num + (1 | subjectID), data = de,
         control = lmerControl(optimizer = "bobyqa")),
    error = function(e) NULL
  )
  if (is.null(fit_ca)) next
  tid_ca <- broom.mixed::tidy(fit_ca, conf.int = TRUE) %>% filter(term == "Condition_num")
  tid_ca$eeg_uid <- uid
  tid_ca$electrodes <- r1$electrodes; tid_ca$fooof <- r1$fooof
  tid_ca$latency_ms <- r1$latency_ms; tid_ca$alpha_type <- r1$alpha_type
  tid_ca$baseline_eeg <- r1$baseline_eeg
  M_cond_alpha_list[[length(M_cond_alpha_list) + 1L]] <- tid_ca
}

M_ca <- bind_rows(M_cond_alpha_list)
if (nrow(M_ca) == 0L) {
  message("No successful condition → alpha fits.")
} else {
  M_ca <- add_sig(M_ca)
  M_ca <- M_ca %>%
    mutate(.lat_ord = match(latency_ms, lat_order),
           .elec_ord = match(electrodes, elec_order)) %>%
    arrange(.lat_ord, .elec_ord, fooof, baseline_eeg, alpha_type) %>%
    select(-.lat_ord, -.elec_ord) %>%
    mutate(ordered_universe = row_number())
  message(sprintf("Condition → alpha: %d EEG universes fitted.", nrow(M_ca)))

  ymax_ca <- max(abs(c(M_ca$conf.low, M_ca$conf.high)), na.rm = TRUE) * 1.05
  ylim_ca <- c(-ymax_ca, ymax_ca)

  p_ca_curve <- ggplot(M_ca, aes(x = ordered_universe, y = estimate, color = condition)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.3, alpha = 0.7) +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    scale_color_manual(values = sig_colors, name = "Significance") +
    guides(alpha = "none") +
    labs(title = expression(bold(alpha ~ "~" ~ condition)),
         x = "Universe", y = "Estimate") +
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
  ggsave(file.path(storage_plot, "AOC_multiverse_nback_condition.png"),
         plot = p_ca_combined, width = 14, height = 10, dpi = 600)
  message("Saved AOC_multiverse_nback_condition.png to ", storage_plot)

  write.csv(M_ca, file.path(csv_dir, "multiverse_nback_condition_results.csv"), row.names = FALSE)
}
