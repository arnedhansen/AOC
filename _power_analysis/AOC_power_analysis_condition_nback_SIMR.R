# AOC Power Analysis — N-Back Condition (SIMR)
# Simulates power for mixed model (alpha ~ condition + (1|Subject)) with simr; effect size d=1.38 (Gevins & Smith 2000). Power curve over N=5–125. Saves figure.
#
# Key outputs:
#   AOC_power_analysis_nback.png; power curve (console)

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
Subject <- factor(rep(1:30, each=10))
condition <- factor(rep(c("A", "B"), 150)) # Lowest and highest condition
alpha_power <- rnorm(300, mean=0, sd=1) + ifelse(condition == "A", 1, 0)
data <- data.frame(Subject, condition, alpha_power)

# Fit the mixed model
model <- lmer(alpha_power ~ condition + (1 | Subject), data=data)
summary(model)

# Extend the model to allow for sample size calculation
extended_model <- extend(model, along="Subject", n=100) # Extend to a larger sample size

# Set the fixed effect size for condition (as estimated from the model)
# Effect size d = 1.809 (0.45 n2p (partial eta squared) from Scharinger et al., 2017)
# Effect size d = 1.38 (0.32 n2p (partial eta squared) from Gevins & Smith, 2000)
fixef(extended_model)["conditionB"] <- 1.38  

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
  file.path(fig_dir, "AOC_power_analysis_nback.csv"),
  row.names = FALSE
)
save_power_curve_plot(
  power_curve_df,
  file.path(fig_dir, "AOC_power_analysis_nback.png")
)