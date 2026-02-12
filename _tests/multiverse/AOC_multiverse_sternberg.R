# AOC Multiverse — Sternberg: alpha ~ gaze * condition + (1|subjectID)
# Reads multiverse_sternberg.csv, fits LMM per universe, plots specification curve + panel.
# Figure output: '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/tests/multiverse'

library(lme4)
library(lmerTest)
library(broom.mixed)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(cowplot)   # for get_legend()

# Paths: CSV from MATLAB; figures path configurable via AOC_MULTIVERSE_FIGURES
script_dir <- Sys.getenv("AOC_MULTIVERSE_DIR", unset = "")
if (nchar(script_dir) == 0L) script_dir <- getwd()
csv_path <- file.path(script_dir, "multiverse_sternberg.csv")
if (!file.exists(csv_path)) stop("CSV not found: ", csv_path, ". Run AOC_multiverse_prep.m or set AOC_MULTIVERSE_DIR.")
storage_plot <- Sys.getenv("AOC_MULTIVERSE_FIGURES", unset = "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/tests/multiverse")
if (!dir.exists(storage_plot)) dir.create(storage_plot, recursive = TRUE)

# Theme (DALAS-style)
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

# Load long data (treat MATLAB NaN as missing so complete.cases and LMM work correctly)
dat <- read.csv(csv_path, stringsAsFactors = FALSE, na.strings = c("NA", "NaN", ""))
dat$subjectID <- as.factor(dat$subjectID)
dat$Condition <- as.factor(dat$Condition)

# Fit LMM per universe; z-score gaze_value and alpha within universe so estimates
# are comparable across gaze measures (fixations ~0-20, SPL ~0-5000, etc.)
universes <- unique(dat$universe_id)
M_list <- list()
for (i in seq_along(universes)) {
  u <- universes[i]
  du <- dat[dat$universe_id == u, ]
  du <- du[complete.cases(du[, c("alpha", "gaze_value", "Condition", "subjectID")]), ]
  if (nrow(du) < 10) next
  # Z-score so estimates are in SD units (comparable across gaze measures)
  du$gaze_value <- as.numeric(scale(du$gaze_value))
  du$alpha      <- as.numeric(scale(du$alpha))
  # Skip universes with zero variance (all identical values)
  if (any(is.nan(du$gaze_value)) || any(is.nan(du$alpha))) next
  fit <- tryCatch(
    lmer(alpha ~ gaze_value * Condition + (1 | subjectID), data = du,
         control = lmerControl(optimizer = "bobyqa")),
    error = function(e) NULL
  )
  if (is.null(fit)) next
  tid <- broom.mixed::tidy(fit, conf.int = TRUE)
  tid$universe_id <- u
  r1 <- du[1, ]
  tid$electrodes <- r1$electrodes
  tid$fooof <- r1$fooof
  tid$latency_ms <- r1$latency_ms
  tid$alpha_type <- r1$alpha_type
  tid$gaze_measure <- r1$gaze_measure
  M_list[[length(M_list) + 1L]] <- tid
}
M_results <- bind_rows(M_list)
if (nrow(M_results) == 0L) stop("No successful LMM fits. Check data and convergence.")

# Keep only gaze_value (main effect) and gaze_value:Condition (interaction) for the curve
terms_keep <- c("gaze_value", grep("gaze_value:Condition", M_results$term, value = TRUE))
terms_keep <- unique(terms_keep)
M_results <- M_results %>% filter(term %in% terms_keep)

# Order universes by estimate of gaze_value (primary term)
term_primary <- "gaze_value"
if (term_primary %in% M_results$term) {
  ord <- M_results %>% filter(term == term_primary) %>% arrange(estimate) %>% pull(universe_id)
} else {
  ord <- unique(M_results$universe_id)
}
ord_df <- data.frame(universe_id = ord, ordered_universe = seq_along(ord))
M_results <- M_results %>% left_join(ord_df, by = "universe_id")

# Significance and color
M_results <- M_results %>%
  mutate(
    condition = case_when(
      p.value < 0.05 & estimate > 0 ~ "Significant Positive",
      p.value < 0.05 & estimate < 0 ~ "Significant Negative",
      TRUE ~ "Non-significant"
    ),
    color = case_when(
      condition == "Significant Positive" ~ "#33CC66",
      condition == "Significant Negative" ~ "#fe0000",
      TRUE ~ "#d1d1d1"
    ),
    alpha_pt = if_else(condition == "Non-significant", 0.5, 0.8)
  )
M_results$condition <- factor(M_results$condition,
  levels = c("Significant Positive", "Significant Negative", "Non-significant"))

# Specification curve: one plot per term (gaze_value, Condition4, Condition6, interactions)
terms_plot <- unique(M_results$term)
max_u <- max(M_results$ordered_universe, na.rm = TRUE)
ylim_est <- c(min(M_results$estimate, na.rm = TRUE) - 0.1, max(M_results$estimate, na.rm = TRUE) + 0.1)

plist <- lapply(terms_plot, function(tr) {
  df <- M_results %>% filter(term == tr)
  ggplot(df, aes(x = ordered_universe, y = estimate, color = condition)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.3, alpha = 0.7) +
    geom_point(aes(alpha = alpha_pt), size = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    scale_color_manual(values = c("#33CC66", "#fe0000", "#d1d1d1"), name = "Significance") +
    guides(alpha = "none") +
    labs(title = tr, x = "Universe", y = "Estimate") +
    theme_minimal() + theme(legend.position = "none") + v_common_theme +
    coord_cartesian(ylim = ylim_est)
})
legend_spec <- get_legend(
  plist[[1]] + theme(legend.position = "bottom") + guides(alpha = "none")
)
p_curve <- wrap_plots(plist, ncol = 2) + plot_layout(guides = "collect")

# Specification panel: tiles (electrodes, fooof, latency_ms, alpha_type, gaze_measure) + N per universe
df_specs <- M_results %>% filter(term == term_primary) %>%
  select(ordered_universe, universe_id, electrodes, fooof, latency_ms, alpha_type, gaze_measure, condition, color)
# N (number of trials) per universe from original data
N_per_u <- dat %>% group_by(universe_id) %>% summarise(N = n(), .groups = "drop")
df_specs <- df_specs %>% left_join(N_per_u, by = "universe_id")
df_long <- df_specs %>%
  pivot_longer(cols = c(electrodes, fooof, latency_ms, alpha_type, gaze_measure),
               names_to = "Variable", values_to = "value")
v_p2_group_order <- c("electrodes", "fooof", "latency_ms", "alpha_type", "gaze_measure")
df_long$Variable <- factor(df_long$Variable, levels = v_p2_group_order)

p_panel <- ggplot(df_long, aes(x = ordered_universe, y = value, fill = condition)) +
  geom_tile() +
  scale_fill_manual(values = c("Significant Positive" = "#33CC66",
                               "Significant Negative" = "#fe0000",
                               "Non-significant" = "#d1d1d1"), name = "Significance") +
  facet_grid(Variable ~ ., scales = "free_y", space = "free") +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0, size = 10), legend.position = "bottom") +
  labs(x = "Universe", y = "Option") + v_common_theme

# N panel: show sample size per universe (DALAS-style: report N in panel)
df_N <- df_specs %>% distinct(ordered_universe, N)
p_N <- ggplot(df_N, aes(x = ordered_universe, y = "N (trials)")) +
  geom_tile(aes(fill = N), width = 1, height = 0.5) +
  scale_fill_gradient(low = "#f0f0f0", high = "#2166ac", name = "N") +
  theme_minimal() + theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 10)) +
  labs(x = "Universe")

# Combined figure (curve, legend, panel, N)
p_combined <- p_curve / legend_spec / p_panel / p_N + plot_layout(heights = c(1.2, 0.15, 1.5, 0.3))
ggsave(file.path(storage_plot, "AOC_multiverse_sternberg.svg"), plot = p_combined, width = 14, height = 12, dpi = 300)
ggsave(file.path(storage_plot, "AOC_multiverse_sternberg.png"), plot = p_combined, width = 14, height = 12, dpi = 300)
message("Saved AOC_multiverse_sternberg.svg and .png to ", storage_plot)

# Optional: save model results CSV for inspection
write.csv(M_results, file.path(script_dir, "multiverse_sternberg_results.csv"), row.names = FALSE)
