# Load necessary libraries
library(ggplot2)
library(ggdist)
library(dplyr)
library(colorspace)  # provides darken(), lighten(), desaturate()

# Define colour palette
pal <- c("#FF8C00", "#A034F0", "#159090")

# Function to add sample size as text (offset adjusted by a fraction of the range)
add_sample <- function(x) {
  offset <- 0.025 * diff(range(x))
  return(c(y = max(x) + offset, label = length(x)))
}

# Read in the data (adjust the file path if needed)
dat <- read.csv("/Volumes/methlab/Students/Arne/AOC/data/features/merged_data_nback.csv")

# Transform reaction times to milliseconds
dat$ReactionTime <- dat$ReactionTime * 1000

# Create the GazeStd column as the average of GazeStdX and GazeStdY
dat$GazeStd <- (dat$GazeStdX + dat$GazeStdY) / 2

# Convert Condition to a factor with labels "1-back", "2-back", "3-back"
dat$Condition <- factor(dat$Condition, levels = c(1, 2, 3), labels = c("1-back", "2-back", "3-back"))

# Define the list of variables, x-axis labels and save names (to be used in file naming)
variables <- c("Accuracy", "ReactionTime", "GazeDeviation", "GazeStd", "MSRate", "Fixations", "Saccades", "AlphaPower", "IAF")
x_labels  <- c("Accuracy [%]", "Reaction Time [ms]", "Gaze Deviation [px]", "Gaze Std [px]",
               "Microsaccade Rate [ms/s]", "Fixations", "Saccades", "Alpha Power [dB]", "IAF [Hz]")
save_names <- c("acc", "rt", "gazedev", "ms", "blink", "fix", "sacc", "pow", "iaf")

# Define the output directory and create it if it doesn't exist
output_dir <- "/Volumes/methlab/Students/Arne/AOC/figures/stats/rainclouds"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop over each variable to create and save raincloud plots
for(i in seq_along(variables)) {
  
  var <- variables[i]
  x_lab <- x_labels[i]
  save_name <- save_names[i]
  
  p <- dat %>% 
    group_by(Condition) %>% 
    ggplot(aes(x = .data[[var]], y = Condition)) +
    # Raincloud (half-eye) plot
    ggdist::stat_halfeye(
      aes(color = Condition,
          fill = after_scale(colorspace::lighten(color, 0.5))),
      adjust = 0.5,
      width = 0.5,
      height = 0.6,
      .width = 0,  # removes the middle line along the data
      justification = -0.4,
      point_colour = NA
    ) +
    # Boxplot layer
    geom_boxplot(
      aes(color = after_scale(colorspace::darken(color, 0.1, space = "HLS")),
          fill = after_scale(colorspace::desaturate(colorspace::lighten(color, 0.8), 0.4))),
      width = 0.35, 
      outlier.shape = NA
    ) +
    # Display data points with jitter on the y-axis
    geom_point(
      aes(color = after_scale(colorspace::darken(color, 0.1, space = "HLS"))),
      fill = "white",
      shape = 21,
      stroke = 0.4,
      size = 0.5,
      position = position_jitter(seed = 1, height = 0.125),
      alpha = 0.4
    ) +
    geom_point(
      aes(fill = Condition),
      color = "transparent",
      shape = 21,
      stroke = 0.4,
      size = 0.3,
      alpha = 0.1,
      position = position_jitter(seed = 1, height = 0.125)
    ) +
    # Add median value as text
    stat_summary(
      geom = "text",
      fun = "median",
      aes(label = round(after_stat(x), 2),
          color = after_scale(colorspace::darken(color, 0.1, space = "HLS"))),
      family = "Roboto Mono",
      fontface = "bold",
      size = 4.5,
      vjust = -3.5
    ) +
    # Add sample size as text
    stat_summary(
      geom = "text",
      fun.data = add_sample,
      aes(label = paste("n =", after_stat(label)),
          color = after_scale(colorspace::darken(color, 0.1, space = "HLS"))),
      family = "Roboto Mono",
      size = 4,
      hjust = -0.5
    ) +
    # Set manual colours for consistency
    scale_color_manual(values = pal, guide = "none") + 
    scale_fill_manual(values = pal, guide = "none") +
    # Define labels and titles
    labs(
      x = x_lab,
      y = NULL,
      title = var,
      subtitle = paste("N-back", var, "by Condition")
    ) +
    theme_minimal(base_family = "Zilla Slab", base_size = 15) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(family = "Roboto Mono"),
      axis.text.y = element_text(
        colour = colorspace::darken(pal, 0.1, space = "HLS"), 
        size = 15
      ),
      axis.title.x = element_text(margin = margin(t = 10), size = 16),
      plot.subtitle = element_text(
        colour = "grey40", hjust = 0,
        margin = margin(0, 0, 20, 0)
      ),
      plot.title.position = "plot",
      plot.margin = margin(15, 15, 10, 15)
    )
  
  # Optionally, adjust the x-axis for ReactionTime specifically (if required)
  if(var == "ReactionTime") {
    p <- p + scale_x_continuous(
      limits = c(200, 2250),
      breaks = seq(250, 2000, by = 250),
      expand = c(0.001, 0.001)
    )
  }
  
  # Save the plot as a PNG file
  ggsave(filename = file.path(output_dir, paste0("AOC_stats_rainclouds_", save_name, "_nback.png")),
         plot = p, width = 8, height = 6, dpi = 300)
}
