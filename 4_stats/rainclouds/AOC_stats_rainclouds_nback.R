# Raincloud plots for AOC N-back data (GLMM-driven asterisks)

# ----------------------------
# Load necessary libraries
# ----------------------------
library(ggplot2)
library(ggdist)
library(dplyr)
library(colorspace)   # darken(), lighten(), desaturate()
library(ggpubr)       # stat_pvalue_manual()
library(lme4)
library(lmerTest)
library(emmeans)
library(stringr)

# ----------------------------
# Colour palette (Perfect AOC pastels)
# ----------------------------
pal <- c("#93B8C4", "#82AD82", "#D998A2")

# ----------------------------
# Helper: map p to significance symbols
# ----------------------------
p_to_signif <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "n.s."
  )
}

# ----------------------------
# Helper: fit GLMM per variable and return df for stat_pvalue_manual
# - Model: response ~ Condition + (1|ID)
# - Pairwise emmeans contrasts with chosen p-adjust
# - y.position is computed relative to data range plus padding
# ----------------------------
glmm_contrasts_df <- function(dat, var, comparisons, y_min, y_max, delta,
                              pad_steps = 0.10, p_adjust = "bonferroni") {
  dvar <- dat %>% dplyr::filter(!is.na(.data[[var]]))
  if (length(unique(dvar$Condition)) < 2L || nrow(dvar) == 0L) {
    return(data.frame())
  }
  
  form <- as.formula(paste(var, "~ Condition + (1|ID)"))
  fit  <- lmer(form, data = dvar, REML = TRUE, na.action = na.omit)
  
  em   <- emmeans::emmeans(fit, ~ Condition)
  pw   <- emmeans::contrast(em, method = "pairwise", adjust = p_adjust)
  pwdf <- as.data.frame(pw)
  
  find_row <- function(df, g1, g2) {
    hits <- with(df, grepl(g1, contrast, fixed = TRUE) & grepl(g2, contrast, fixed = TRUE))
    df[hits, , drop = FALSE] |> dplyr::slice(1)
  }
  
  out <- dplyr::bind_rows(lapply(seq_along(comparisons), function(i) {
    cpair <- comparisons[[i]]
    row   <- find_row(pwdf, cpair[1], cpair[2])
    if (nrow(row) == 0) return(NULL)
    y_pos <- y_max + (i * pad_steps * (y_max - y_min) + delta)
    data.frame(
      group1     = cpair[1],
      group2     = cpair[2],
      p.adj      = row$p.value,
      p.signif   = p_to_signif(row$p.value),
      y.position = y_pos,
      stringsAsFactors = FALSE
    )
  }))
  out
}

# ----------------------------
# Read in the data
# ----------------------------
dat <- read.csv("/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/merged_data_nback.csv")

# Transform and derive columns
#dat$ReactionTime <- dat$ReactionTime * 1000
dat$GazeStd      <- (dat$GazeStdX + dat$GazeStdY) / 2
dat$Condition    <- factor(dat$Condition,
                           levels = c(1,2,3),
                           labels = c("1-back","2-back","3-back"))

# Ensure ID is a factor (needed for random intercepts)
if (!is.factor(dat$ID)) dat$ID <- factor(dat$ID)

# Variables, labels and save names
variables  <- c("Accuracy","ReactionTime","GazeDeviation","GazeStd",
                "MSRate","Fixations","Saccades","AlphaPower","IAF")
titles     <- c("Accuracy","Reaction Time","Gaze Deviation","Gaze Std",
                "Microsaccade Rate","Fixations","Saccades","Alpha Power","IAF")
y_labels   <- c("Accuracy [%]","Reaction Time [s]","Gaze Deviation [px]","Gaze Std [px]",
                "Microsaccade Rate [MS/s]","Fixations","Saccades",
                "Alpha Power [\u03BCV²/Hz]","IAF [Hz]")
save_names <- c("acc","rt","gazedev","gazestd","ms","fix","sacc","pow","iaf")

# Remove outliers by 1.5×IQR within each condition
dat <- dat %>%
  group_by(Condition) %>%
  mutate(across(all_of(variables), ~{
    lo <- quantile(.x, .25, na.rm=TRUE) - 1.5*IQR(.x, na.rm=TRUE)
    hi <- quantile(.x, .75, na.rm=TRUE) + 1.5*IQR(.x, na.rm=TRUE)
    replace(.x, .x<lo|.x>hi, NA)
  })) %>%
  ungroup()

# Create a tiny jittered copy of Accuracy for half-eye density only
set.seed(123)
dat2 <- dat %>%
  filter(!is.na(Accuracy)) %>%
  mutate(Accuracy_jit = Accuracy + runif(n(), -1.15, 1.15))

# Output directory
output_dir <- "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/stats/rainclouds"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Pairwise comparisons for stats
comparisons <- list(
  c("1-back","2-back"),
  c("1-back","3-back"),
  c("2-back","3-back")
)

# ----------------------------
# Loop through variables
# ----------------------------
for (i in seq_along(variables)) {
  var       <- variables[i]
  y_lab     <- y_labels[i]
  save_name <- save_names[i]
  # work on non-missing rows for this variable
  dvar <- dat %>% dplyr::filter(!is.na(.data[[var]]))
  if (nrow(dvar) == 0L) next  # nothing to plot
  
  # data-driven limits and padding from dvar
  y_min   <- min(dvar[[var]], na.rm = TRUE)
  y_max   <- max(dvar[[var]], na.rm = TRUE)
  rng     <- y_max - y_min
  delta   <- 0.05 * rng
  y_extra <- 0.25 * rng
  
  # ----------------------------
  # Compute data-driven limits and padding
  # ----------------------------
  y_min  <- min(dat[[var]], na.rm = TRUE)
  y_max  <- max(dat[[var]], na.rm = TRUE)
  rng    <- y_max - y_min
  delta  <- 0.025 * rng
  y_extra <- 0.25 * rng   # vertical headroom so brackets don't get cropped
  
  # Variable-specific lower bounds (to keep your previous look)
  lower_bound <-
    if (var == "Accuracy")       65 else
        if (var == "GazeDeviation")   5 else y_min
  
  # Optional variable-specific nominal uppers; we still ensure headroom
  nominal_upper <-
      if (var == "GazeDeviation")  65 else
        if (var == "Accuracy")      102 else NA_real_
  
  upper_bound <- max(c(y_max + y_extra, nominal_upper), na.rm = TRUE)
  
  # ----------------------------
  # BASE PLOT (no stats)
  # ----------------------------
  p_base <- dvar %>%
    ggplot(aes(x = Condition, y = .data[[var]])) +
    
    # half-eye densities (use jittered Accuracy for nicer density)
    {
      if (var == "Accuracy") {
        stat_halfeye(
          data = dat2,
          aes(y     = Accuracy_jit,
              colour = Condition,
              fill   = after_scale(lighten(colour, 0.5))),
          adjust        = 0.5,
          width         = 0.6,
          scale         = 1,
          justification = 1.5,
          side          = "left",
          alpha         = 0.5
        )
      } else {
        stat_halfeye(
          aes(colour = Condition,
              fill   = after_scale(lighten(colour, 0.5))),
          adjust        = 0.5,
          width         = 0.5,
          justification = 1.55,
          side          = "left",
          alpha         = 0.5
        )
      }
    } +
    
    # boxplot on the real data
    geom_boxplot(
      aes(color = Condition, fill = after_scale(lighten(color, 0.5))),
      width = 0.35, outlier.shape = NA
    ) +
    
    # points on the real data
    geom_point(
      aes(color = Condition, fill = after_scale(lighten(color, 0.5))),
      shape = 21, stroke = 0.4, size = 3,
      position = position_jitter(seed = 1, width = 0.125),
      alpha = 0.65
    ) +
    geom_point(
      aes(fill = Condition),
      color = "transparent", shape = 21, stroke = 0.4, size = 0.25,
      alpha = 0.4, position = position_jitter(seed = 1, width = 0.125)
    ) +
    
    # colours & guides
    scale_color_manual(values = pal, guide = "none") +
    scale_fill_manual(values  = pal, guide = "none") +
    
    # labels
    labs(
      x        = "Condition",
      y        = y_lab,
      title    = titles[i],
      subtitle = paste("N-back", titles[i], "by Condition")
    ) +
    
    # theme
    theme_minimal(base_family = "Zilla Slab", base_size = 15) +
    theme(
      plot.background    = element_rect(fill = "white", colour = NA),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.ticks         = element_blank(),
      axis.text.x        = element_text(family = "Roboto Mono"),
      axis.text.y        = element_text(family = "Roboto Mono", size = 15),
      axis.title.y       = element_text(margin = margin(r = 10), size = 16),
      plot.subtitle      = element_text(colour = "grey40", hjust = 0,
                                        margin = margin(0,0,20,0)),
      plot.title.position= "plot",
      plot.margin        = margin(15,15,10,15)
    ) +
    
    # y-scale and zoom (base plot: no extra headroom needed)
    scale_y_continuous(expand = expansion(add = 0)) +
    coord_cartesian(
      ylim = c(lower_bound, ifelse(is.finite(nominal_upper), nominal_upper, y_max)),
      clip = "off"
    )
  
  # save the base raincloud
  ggsave(
    filename = file.path(output_dir,
                         paste0("AOC_stats_rainclouds_",
                                save_name, "_nback.png")),
    plot   = p_base, width = 8, height = 6, dpi = 300
  )
  
  # ----------------------------
  # GLMM + emmeans contrasts for this variable
  # ----------------------------
  glmm_df <- glmm_contrasts_df(
    dat         = dvar,
    var         = var,
    comparisons = comparisons,
    y_min       = y_min,
    y_max       = y_max,
    delta       = delta,
    pad_steps   = 0.10,          # similar to step.increase
    p_adjust    = "bonferroni"
  )
  # ---- ensure the panel is tall enough for brackets ----
  max_bracket <- if (!is.null(glmm_df) && nrow(glmm_df) > 0) {
    max(glmm_df$y.position, na.rm = TRUE)
  } else {
    y_max + y_extra
  }
  upper_needed <- max_bracket + 0.02 * (y_max - y_min)
  upper_bound  <- max(c(upper_needed, nominal_upper, y_max + y_extra), na.rm = TRUE)
  
  # ----------------------------
  # STATS PLOT (driven by GLMM p-values)
  # - Give headroom so brackets are never cropped
  # ----------------------------
  p_stats <- p_base +
    labs(title = "", subtitle = "") +
    coord_cartesian(ylim = c(lower_bound, upper_bound), clip = "off") +
    theme(plot.margin = margin(20 + delta*10, 15, 10, 15))
  
  if (!is.null(glmm_df) && nrow(glmm_df) > 0) {
    p_stats <- p_stats +
      ggpubr::stat_pvalue_manual(
        data        = glmm_df,
        label       = "p.signif",
        xmin        = "group1",
        xmax        = "group2",
        y.position  = "y.position",
        tip.length  = 0.01,
        bracket.size= 0.6,
        size        = 5,
        family      = "Roboto Mono",
        colour      = "grey30"
      )
  }
  
  # save the stats raincloud
  ggsave(
    filename = file.path(output_dir,
                         paste0("AOC_stats_rainclouds_",
                                save_name, "_nback_stats.png")),
    plot   = p_stats, width = 8, height = 6, dpi = 300
  )
}
