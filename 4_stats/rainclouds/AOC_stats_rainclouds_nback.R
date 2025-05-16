# Raincloud plots for AOC N-back data

# Load necessary libraries
library(ggplot2)
library(ggdist)
library(dplyr)
library(colorspace)   # provides darken(), lighten(), desaturate()
library(rstatix)      # for rm-ANOVA and pairwise tests
library(ggpubr)       # for stat_compare_means()

# Define colour palette (Perfect AOC pastels)
pal <- c("#93B8C4", "#82AD82", "#D998A2")

# Helper to add sample size as text (unused in this version)
add_sample <- function(x) {
  offset <- 0.025 * diff(range(x))
  return(c(y = max(x) + offset, label = length(x)))
}

# Read in the data
dat <- read.csv("/Volumes/methlab/Students/Arne/AOC/data/features/merged_data_nback.csv")

# Transform and derive columns
dat$ReactionTime <- dat$ReactionTime * 1000
dat$GazeStd      <- (dat$GazeStdX + dat$GazeStdY) / 2
dat$Condition    <- factor(dat$Condition,
                           levels = c(1,2,3),
                           labels = c("1-back","2-back","3-back"))

# Variables, labels and save names
variables  <- c("Accuracy","ReactionTime","GazeDeviation","GazeStd",
                "MSRate","Fixations","Saccades","AlphaPower","IAF")
titles     <- c("Accuracy","Reaction Time","Gaze Deviation","Gaze Std",
                "Microsaccade Rate","Fixations","Saccades","Alpha Power","IAF")
y_labels   <- c("Accuracy [%]","Reaction Time [ms]","Gaze Deviation [px]","Gaze Std [px]",
                "Microsaccade Rate [ms/s]","Fixations","Saccades",
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

# Output directory
output_dir <- "/Volumes/methlab/Students/Arne/AOC/figures/stats/rainclouds"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Pairwise comparisons for stats
comparisons <- list(
  c("1-back","2-back"),
  c("1-back","3-back"),
  c("2-back","3-back")
)

# Loop through variables
for (i in seq_along(variables)) {
  
  var       <- variables[i]
  y_lab     <- y_labels[i]
  save_name <- save_names[i]
  
  ##### BASE PLOT #####
  p_base <- dat %>%
    ggplot(aes(x = Condition, y = .data[[var]])) +
    stat_halfeye(
      aes(color = Condition, fill = after_scale(lighten(color,0.5))),
      adjust        = 0.5, width = 0.5, height = 0.6,
      .width        = 0,
      justification = 1.55, side = "left"
    ) +
    geom_boxplot(
      aes(color = Condition, fill = after_scale(lighten(color,0.5))),
      width = 0.35, outlier.shape = NA
    ) +
    geom_point(
      aes(color = Condition, fill = after_scale(lighten(color,0.5))),
      shape = 21, stroke = 0.4, size = 3,
      position = position_jitter(seed=1, width=0.125),
      alpha = 0.65
    ) +
    geom_point(
      aes(fill = Condition),
      color = "transparent", shape = 21, stroke = 0.4, size = 0.25,
      alpha = 0.4, position = position_jitter(seed=1, width=0.125)
    ) +
    scale_color_manual(values = pal, guide = "none") +
    scale_fill_manual(values = pal, guide = "none") +
    labs(
      x        = "Condition",
      y        = y_lab,
      title    = titles[i],
      subtitle = paste("N-back", titles[i], "by Condition")
    ) +
    theme_minimal(base_family="Zilla Slab", base_size=15) +
    theme(
      plot.background    = element_rect(fill="white", colour=NA),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.ticks         = element_blank(),
      axis.text.x        = element_text(family="Roboto Mono"),
      axis.text.y        = element_text(family="Roboto Mono", size=15),
      axis.title.y       = element_text(margin=margin(r=10), size=16),
      plot.subtitle      = element_text(colour="grey40", hjust=0,
                                        margin=margin(0,0,20,0)),
      plot.title.position= "plot",
      plot.margin        = margin(15,15,10,15)
    )
  
  # Custom y-limits for base if desired
  if (var=="Accuracy") {
    p_base <- p_base +
      scale_y_continuous(
        limits = c(70,100),
        breaks = seq(75,100,5),
        expand = c(0.001,0.001)
      )
  }
  # (repeat for GazeDeviation, MSRate, Fixations, IAF as in your script…)
  
  # Save base plot
  ggsave(
    filename = file.path(output_dir,
                         paste0("AOC_stats_rainclouds_", save_name, "_nback.png")),
    plot   = p_base, width = 8, height = 6, dpi = 300
  )
  
  ##### STATS PLOT #####
  # Run repeated-measures ANOVA and pairwise tests
  anova_res <- dat %>%
    anova_test(dv = .data[[var]], wid = ID, within = Condition)
  print(glue::glue("ANOVA for {var}:")); print(anova_res)
  
  pwc <- dat %>%
    pairwise_t_test(
      formula         = as.formula(paste(var, "~ Condition")),
      paired          = TRUE,
      p.adjust.method = "bonferroni"
    )
  print(glue::glue("Pairwise tests for {var}:")); print(pwc)
  
  # Compute data range and buffer for brackets
  y_min <- min(dat[[var]], na.rm=TRUE)
  y_max <- max(dat[[var]], na.rm=TRUE)
  delta <- 0.05 * (y_max - y_min)
  
  # Build stats plot from p_base
  p_stats <- p_base +
    labs(title="", subtitle="") +
    coord_cartesian(ylim = c(y_min, y_max), clip = "off") +
    stat_compare_means(
      comparisons     = comparisons,
      method          = "t.test",
      paired          = TRUE,
      p.adjust.method = "bonferroni",
      label           = "p.signif",
      symnum.args     = list(
        cutpoints = c(0,0.001,0.01,0.05,1),
        symbols   = c("***","**","*","n.s.")
      ),
      label.size      = 5,
      family          = "Roboto Mono",
      colour          = "grey30",
      tip.length      = 0.01,
      bracket.size    = 0.6,
      step.increase   = 0.1,
      hide.ns         = FALSE
    ) +
    # give 10% breathing room above max for brackets
    scale_y_continuous(expand = expansion(mult = c(0, .10))) +
    # adjust top margin for the annotation strip
    theme(plot.margin = margin(20 + delta*10, 15, 10, 15))
  
  # Special margin and y-limits for RT if desired
  if (var=="ReactionTime") {
    p_stats <- p_stats +
      theme(plot.margin = margin(20,30,10,15)) +
      scale_y_continuous(
        limits = c(550,1300),
        breaks = seq(600,1200,100),
        expand = c(0.001,0.001)
      )
  }
  
  # Save stats plot
  ggsave(
    filename = file.path(output_dir,
                         paste0("AOC_stats_rainclouds_", save_name, "_nback_stats.png")),
    plot   = p_stats, width = 8, height = 6, dpi = 300
  )
}
