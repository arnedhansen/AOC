# AOC Multiverse — N-Back (OLD 5-dimension grid, 192 universes)
# Reads multiverse_nback.csv from the original AOC_multiverse_prep.m (no baseline/freq_method columns).
# Fits LMM per universe, plots specification curve + panel (DALAS-style).
# Figure output: '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/tests/multiverse'

library(lme4)
library(lmerTest)
library(broom.mixed)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(cowplot)   # for get_legend()

# Paths: CSVs live on the server features directory; figures path configurable
csv_dir <- Sys.getenv("AOC_MULTIVERSE_DIR", unset = "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features")
csv_path <- file.path(csv_dir, "multiverse_nback_OLD.csv")
if (!file.exists(csv_path)) stop("CSV not found: ", csv_path, ". Run AOC_multiverse_prep.m or set AOC_MULTIVERSE_DIR.")
storage_plot <- Sys.getenv("AOC_MULTIVERSE_FIGURES", unset = "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/tests/multiverse")
if (!dir.exists(storage_plot)) dir.create(storage_plot, recursive = TRUE)

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

# Treat MATLAB NaN as missing so complete.cases and LMM work correctly
dat <- read.csv(csv_path, stringsAsFactors = FALSE, na.strings = c("NA", "NaN", ""))
dat$subjectID <- as.factor(dat$subjectID)
dat$Condition <- as.factor(dat$Condition)

# Z-score gaze_value and alpha within each universe so estimates are in SD units
universes <- unique(dat$universe_id)
M_list <- list()
for (i in seq_along(universes)) {
  u <- universes[i]
  du <- dat[dat$universe_id == u, ]
  du <- du[complete.cases(du[, c("alpha", "gaze_value", "Condition", "subjectID")]), ]
  if (nrow(du) < 10) next
  du$gaze_value <- as.numeric(scale(du$gaze_value))
  du$alpha      <- as.numeric(scale(du$alpha))
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

term_primary <- "gaze_value"
if (term_primary %in% M_results$term) {
  ord <- M_results %>% filter(term == term_primary) %>% arrange(estimate) %>% pull(universe_id)
} else {
  ord <- unique(M_results$universe_id)
}
ord_df <- data.frame(universe_id = ord, ordered_universe = seq_along(ord))
M_results <- M_results %>% left_join(ord_df, by = "universe_id")

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

# Specification curve: single plot for gaze_value main effect only
# (interaction terms are still available in the results CSV for inspection)
df_curve <- M_results %>% filter(term == term_primary)
ylim_est <- c(min(df_curve$estimate, na.rm = TRUE) - 0.1, max(df_curve$estimate, na.rm = TRUE) + 0.1)

p_curve <- ggplot(df_curve, aes(x = ordered_universe, y = estimate, color = condition)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.3, alpha = 0.7) +
  geom_point(aes(alpha = alpha_pt), size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  scale_color_manual(values = c("#33CC66", "#fe0000", "#d1d1d1"), name = "Significance") +
  guides(alpha = "none") +
  labs(title = "gaze_value (main effect)", x = "Universe", y = "Estimate") +
  theme_minimal() + theme(legend.position = "none") + v_common_theme +
  coord_cartesian(ylim = ylim_est)
legend_spec <- get_legend(
  p_curve + theme(legend.position = "bottom") + guides(alpha = "none")
)

df_specs <- M_results %>% filter(term == term_primary) %>%
  select(ordered_universe, universe_id, electrodes, fooof, latency_ms, alpha_type, gaze_measure, condition, color)
N_per_u <- dat %>% group_by(universe_id) %>% summarise(N = n(), .groups = "drop")
df_specs <- df_specs %>% left_join(N_per_u, by = "universe_id")
df_long <- df_specs %>%
  pivot_longer(cols = c(electrodes, fooof, latency_ms, alpha_type, gaze_measure),
               names_to = "Variable", values_to = "value")
df_long$Variable <- factor(df_long$Variable,
  levels = c("electrodes", "fooof", "latency_ms", "alpha_type", "gaze_measure"))

p_panel <- ggplot(df_long, aes(x = ordered_universe, y = value, fill = condition)) +
  geom_tile() +
  scale_fill_manual(values = c("Significant Positive" = "#33CC66",
                               "Significant Negative" = "#fe0000",
                               "Non-significant" = "#d1d1d1"), name = "Significance") +
  facet_grid(Variable ~ ., scales = "free_y", space = "free") +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0, size = 10), legend.position = "bottom") +
  labs(x = "Universe", y = "Option") + v_common_theme

df_N <- df_specs %>% distinct(ordered_universe, N)
p_N <- ggplot(df_N, aes(x = ordered_universe, y = "N (trials)")) +
  geom_tile(aes(fill = N), width = 1, height = 0.5) +
  scale_fill_gradient(low = "#f0f0f0", high = "#2166ac", name = "N") +
  theme_minimal() + theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 10)) +
  labs(x = "Universe")

p_combined <- p_curve / legend_spec / p_panel / p_N + plot_layout(heights = c(0.8, 0.1, 1.5, 0.3))
ggsave(file.path(storage_plot, "AOC_multiverse_nback_OLD.svg"), plot = p_combined, width = 14, height = 12, dpi = 300)
ggsave(file.path(storage_plot, "AOC_multiverse_nback_OLD.png"), plot = p_combined, width = 14, height = 12, dpi = 300)
message("Saved AOC_multiverse_nback_OLD.svg and .png to ", storage_plot)

write.csv(M_results, file.path(csv_dir, "multiverse_nback_OLD_results.csv"), row.names = FALSE)
