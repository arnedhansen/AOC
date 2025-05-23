# Raincloud plots for AOC Sternberg data

# Load necessary libraries
library(ggplot2)
library(ggdist)
library(dplyr)
library(colorspace)   # provides darken(), lighten(), desaturate()
library(rstatix)      # for rm-ANOVA and pairwise tests
library(ggpubr)       # for stat_compare_means()

# Define colour palette
pal <- c("#93B8C4", "#82AD82", "#D998A2") # Perfect AOC pastel colours

# Function to add sample size as text (offset adjusted by a fraction of the range)
add_sample <- function(x) {
  offset <- 0.025 * diff(range(x))
  return(c(y = max(x) + offset, label = length(x)))
}

# Read in the data (adjust the file path if needed)
dat <- read.csv("/Volumes/methlab/Students/Arne/AOC/data/features/merged_data_sternberg.csv")

# Change condition from 2, 4, 6 to 1, 2, 4
dat$Condition <- dat$Condition / 2

# Transform reaction times to milliseconds
dat$ReactionTime <- dat$ReactionTime * 1000

# Create the GazeStd column as the average of GazeStdX and GazeStdY
dat$GazeStd <- (dat$GazeStdX + dat$GazeStdY) / 2

# Convert Condition to a factor with labels 
dat$Condition <- factor(dat$Condition, levels = c(1, 2, 3), labels = c("WM load 2", "WM load 4", "WM load 6"))

# Define the list of variables, y-axis labels and save names
variables  <- c("Accuracy", "ReactionTime", "GazeDeviation", "GazeStd", "MSRate", "Fixations", "Saccades", "AlphaPower", "IAF")
titles     <- c("Accuracy", "Reaction Time", "Gaze Deviation", "Gaze Standard Deviation",
                "Microsaccade Rate", "Fixations", "Saccades", "Alpha Power", "Individual Alpha Frequency")
y_labels   <- c("Accuracy [%]", "Reaction Time [ms]", "Gaze Deviation [px]", "Gaze Std [px]",
                "Microsaccade Rate [MS/s]", "Fixations", "Saccades", "Alpha Power [\u03BCV²/Hz]", "IAF [Hz]")
save_names <- c("acc", "rt", "gazedev", "gazestd", "ms", "fix", "sacc", "pow", "iaf")

# Remove outliers using the IQR method (1.5 * IQR rule)
dat <- dat %>%
  group_by(Condition) %>%
  mutate(across(all_of(variables), ~{
    lower <- quantile(.x, 0.25, na.rm = TRUE) - 1.5 * IQR(.x, na.rm = TRUE)
    upper <- quantile(.x, 0.75, na.rm = TRUE) + 1.5 * IQR(.x, na.rm = TRUE)
    replace(.x, .x < lower | .x > upper, NA)
  })) %>%
  ungroup()

# Define the output directory and create it if it doesn't exist
output_dir <- "/Volumes/methlab/Students/Arne/AOC/figures/stats/rainclouds"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Specify all pairwise comparisons
comparisons <- list(
  c("WM load 2", "WM load 4"),
  c("WM load 2", "WM load 6"),
  c("WM load 4", "WM load 6")
)

# Loop over each variable to create and save raincloud plots
for(i in seq_along(variables)) {
  
  var       <- variables[i]
  y_lab     <- y_labels[i]
  save_name <- save_names[i]
  
  ##### BASE PLOT #####
  # Build the base raincloud plot (no stats, no significance annotations)
  p_base <- dat %>% 
    ggplot(aes(x = Condition, y = .data[[var]])) +
    
    # Raincloud (half-eye)
    stat_halfeye(
      aes(color = Condition,
          fill  = after_scale(lighten(color, 0.5))),
      adjust        = 0.5,
      width         = 0.5,
      height        = 0.6,
      .width        = 0,
      justification = 1.55,
      side          = "left"
    ) +
    
    # Boxplot
    geom_boxplot(
      aes(color = Condition,
          fill  = after_scale(lighten(color, 0.5))),
      width         = 0.35, 
      outlier.shape = NA
    ) +
    
    # Jittered points
    geom_point(
      aes(color = Condition,
          fill  = after_scale(lighten(color, 0.5))),
      shape         = 21,
      stroke        = 0.4,
      size          = 3,
      position      = position_jitter(seed = 1, width = 0.125),
      alpha         = 0.65
    ) +
    geom_point(
      aes(fill = Condition),
      color    = "transparent",
      shape    = 21,
      stroke   = 0.4,
      size     = 0.25,
      alpha    = 0.4,
      position = position_jitter(seed = 1, width = 0.125)
    ) +
    
    # Labels, title and theme
    scale_color_manual(values = pal, guide = "none") + 
    scale_fill_manual(values = pal, guide = "none") +
    labs(
      x        = "Condition",
      y        = y_lab,
      title    = titles[i],
      subtitle = paste("Sternberg", titles[i], "by Condition")
    ) +
    theme_minimal(base_family = "Zilla Slab", base_size = 15) +
    theme(
      plot.background      = element_rect(fill = "white", colour = NA),
      panel.grid.minor     = element_blank(),
      panel.grid.major.x   = element_blank(),
      axis.ticks           = element_blank(),
      axis.text.x          = element_text(family = "Arial"),
      axis.text.y          = element_text(family = "Arial", size = 15),
      axis.title.y         = element_text(margin = margin(r = 10), size = 16),
      plot.subtitle        = element_text(colour = "grey40", hjust = 0,
                                          margin = margin(0, 0, 20, 0)),
      plot.title.position  = "plot",
      plot.margin          = margin(15, 15, 10, 15)
    )
  
  # Save the base plot
  ggsave(
    filename = file.path(output_dir, paste0("AOC_stats_rainclouds_", save_name, "_sternberg.png")),
    plot     = p_base,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
  ##### STATS #####
  # Repeated‐measures ANOVA (printed to console)
  anova_res <- dat %>%
    anova_test(
      dv      = .data[[var]],
      wid     = ID,
      within  = Condition
    )
  print(glue::glue("ANOVA for {var}:"))
  print(anova_res)
  
  # Pairwise paired t-tests with Bonferroni correction (printed to console)
  pwc <- dat %>%
    pairwise_t_test(
      formula         = as.formula(paste(var, "~ Condition")),
      paired          = TRUE,
      p.adjust.method = "bonferroni"
    )
  print(glue::glue("Pairwise tests for {var}:"))
  print(pwc)
  
  # Determine natural y-limits and a small delta for the annotation strip
  y_min <- min(dat[[var]], na.rm = TRUE)
  y_max <- max(dat[[var]], na.rm = TRUE)
  delta <- 0.05 * (y_max - y_min)
  
  ##### STATS PLOT #####
  # Build the stats plot starting from p_base
  p_stats <- p_base +
    labs(title    = "",
         subtitle = "") +
    # fix y to data min/max and allow drawing outside
    coord_cartesian(ylim = c(y_min, y_max), clip = "off") +
    # Significance annotations
    stat_compare_means(
      comparisons      = comparisons,
      method           = "t.test",
      paired           = TRUE,
      p.adjust.method  = "bonferroni",
      label            = "p.signif",
      symnum.args      = list(
        cutpoints = c(0, 0.001, 0.01, 0.05, 1),
        symbols   = c("***", "**", "*", "n.s.")
      ),
      label.size       = 5,
      family           = "Roboto Mono",
      colour           = "grey30",
      tip.length       = 0.01,
      bracket.size     = 0.6,
      step.increase    = 0.1,
      hide.ns          = FALSE
    )
  
  # Add plot margins
  if (var == "Accuracy") {
    p_stats <- p_stats +
      theme(plot.margin = margin(30, 50, 10, 15)) +
      scale_y_continuous(
        limits = c(65, 140),
        breaks = seq(70, 100, by = 5),
        expand = c(0.001, 0.001)
      )
  } else if (var == "ReactionTime") {
    p_stats <- p_stats +
      theme(plot.margin = margin(50, 75, 50, 75)) +
      coord_cartesian(ylim = c(300, 1300), clip = "off") +
      scale_y_continuous(
        breaks = seq(400, 1200, by = 100),
        expand = c(0.001, 0.001)
      )
   } else if (var == "GazeDeviation") {
      p_stats <- p_stats +
        theme(plot.margin = margin(50, 75, 10, 15)) +
        scale_y_continuous(
          limits = c(5, 175),
          breaks = seq(25, 125, by = 25),
          expand = c(0.001, 0.001)
        )
  } else {
    p_stats <- p_stats +
      theme(plot.margin = margin(20 + delta * 10, 15, 10, 15)) +
      scale_y_continuous(expand = expansion(mult = c(0, .10)))
  }
  
  # Save the stats plot
  ggsave(
    filename = file.path(output_dir, paste0("AOC_stats_rainclouds_", save_name, "_sternberg_stats.png")),
    plot     = p_stats,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}
