# AOC Power Analysis — Gaze (SIMR)
# Simulates power for mixed model (alpha ~ gaze + (1|Subject)) with simr; effect size d=0.4 (Popov et al.). Power curve over N=5–125. Saves figure.
#
# Key outputs:
#   AOC_power_analysis_gaze.png; power curve (console)

# Tutorial: https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12504
library(lme4)
library(simr)
library(ggplot2)

script_dir <- dirname(normalizePath(
  sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
))
if (length(script_dir) == 0 || is.na(script_dir) || script_dir == "") {
  script_dir <- getwd()
}
source(file.path(script_dir, "AOC_power_curve_plot.R"))

fig_dir <- "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/power_analysis"

# Example dataset
set.seed(123)
Subject <- factor(rep(1:10))
gaze <- rnorm(20, mean=0, sd=1)
alpha_power <- rnorm(20, mean=0, sd=1)
data <- data.frame(Subject, gaze, alpha_power)

# Fit the mixed model
model <- lmer(alpha_power ~ gaze + (1 | Subject), data=data)
summary(model)

# Extend the model to allow for sample size calculation
extended_model <- extend(model, along = "Subject", n=100) # Extend to a larger sample size

# Set the fixed effect size for condition (as estimated from the model)
# Effect size of gaze on alpha: d = 0.4 (from Popov et al., 2021b)
fixef(extended_model)["gaze"] <- 0.4  

# Power analysis
sim <- powerSim(extended_model, nsim = 1000, progress = FALSE)

# Compute power curve for a range of sample sizes
power_curve <- powerCurve(extended_model, along = "Subject", nsim = 1000, breaks = seq(5, 125, by=5))

# Print, plot and save the power curve
print(power_curve)
summary(power_curve)
power_curve_summary <- summary(power_curve)
power_curve_df <- data.frame(
  Subjects = power_curve_summary$nlevels,
  Power = power_curve_summary$mean,
  Lower = power_curve_summary$lower,
  Upper = power_curve_summary$upper
)
write.csv(
  power_curve_df,
  file.path(fig_dir, "AOC_power_analysis_gaze.csv"),
  row.names = FALSE
)
save_power_curve_plot(
  power_curve_df,
  file.path(fig_dir, "AOC_power_analysis_gaze.png")
)