# AOC Sternberg EEG Alpha Power Lateralization

# Load necessary libraries
library(ggplot2)
library(ggdist)
library(dplyr)
library(colorspace)  # provides darken(), lighten(), desaturate()

# Define colour palette
pal <- c("#93B8C4", "#82AD82", "#D998A2") 

# Read in the data
dat <- read.csv("/Volumes/methlab/Students/Arne/AOC/data/features/merged_data_sternberg.csv")

# Remove outliers by setting any points over 0.55 or -0.55 in dat$Lateralization to NA
dat$Lateralization[dat$Lateralization > 0.55 | dat$Lateralization < -0.55] <- NA

# Convert Condition to a factor with labels 
dat$Condition <- dat$Condition / 2
dat$Condition <- factor(dat$Condition, levels = c(1, 2, 3), labels = c("WM load 2", "WM load 4", "WM load 6"))

# Set variables for lateralisation visualisation
var <- "Lateralization"
title <- "Lateralization"
subtit <- "Alpha Power Lateralization"
x_lab <- "Lateralization [%]"  # continuous variable on x-axis
y_lab <- "Condition"
save_name <- "lat"

# Define the output directory and create it if it doesn't exist
output_dir <- "/Volumes/methlab/Students/Arne/AOC/figures/eeg/lateralization"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Create horizontal raincloud plot (lateralisation on x-axis)
dat$Condition_num <- as.numeric(dat$Condition)
p <- dat %>% 
  group_by(Condition) %>% 
  ggplot(aes(x = .data[[var]], y = Condition_num - 0.125)) +
  
  # Half-eye density layer adjusted for horizontal orientation
  ggdist::stat_halfeye(
    aes(color = Condition,
        fill = after_scale(lighten(color, 0.5))),
    adjust = 0.25,
    width = 0.2,
    height = .55,
    .width = 0,               # removes the central line
    justification = -.5,     # shifts the half-eye relative to the boxplot
    side = "right"            # density drawn to the right of the boxplot
  ) +
  
  # Horizontal boxplot layer
  geom_boxplot(
    aes(color = Condition,
        fill = after_scale(lighten(color, 0.5))),
    width = 0.35, 
    outlier.shape = NA,
    #position = position_nudge(y = -0.125), 
  ) +
  
  # Jittered data points; jitter applied in the y-axis (Condition) to avoid overlap
  geom_point(
    aes(color = Condition,
        fill = after_scale(lighten(color, 0.5))),
    shape = 21,
    stroke = 0.4,
    size = 3,
    position = position_jitter(seed = 1, height = .125),
    alpha = 0.65
  ) +
  geom_point(
    aes(fill = Condition),
    color = "transparent",
    shape = 20,
    stroke = 0.4,
    size = .75,
    alpha = 0.4,
    position = position_jitter(seed = 1, height = .125)
  ) +
  
  # Set manual colours for consistency
  scale_color_manual(values = pal, guide = "none") + 
  scale_fill_manual(values = pal, guide = "none") +
  
  # Adjust y-axis to show original condition labels
  scale_y_continuous(
    breaks = 1:length(levels(dat$Condition)),
    labels = levels(dat$Condition)
  ) +
  
  # Define labels and titles
  labs(
    x = x_lab,
    y = y_lab,
    title = title,
    subtitle = paste("Sternberg", subtit, "by Condition")
  ) +
  
  theme_minimal(base_family = "Zilla Slab", base_size = 15) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),  # remove horizontal grid lines
    axis.ticks = element_blank(),
    axis.text.x = element_text(family = "Roboto Mono", size = 15),
    axis.text.y = element_text(family = "Roboto Mono"),
    axis.title.x = element_text(margin = margin(b = 10), size = 16),
    plot.subtitle = element_text(
      colour = "grey40", hjust = 0,
      margin = margin(0, 0, 20, 0)
    ),
    plot.title.position = "plot",
    plot.margin = margin(15, 15, 10, 15)
  )

# Adjust x-axis limits and breaks
p <- p + scale_x_continuous(
  limits = c(-.55, .55),
  breaks = seq(-.5, .5, by = .25),
  expand = c(0.001, 0.001)
) 

# Save the plot as a PNG file
ggsave(filename = file.path(output_dir, paste0("AOC_stats_rainclouds_", save_name, "_sternberg.png")),
       plot = p, width = 8, height = 6, dpi = 300)
