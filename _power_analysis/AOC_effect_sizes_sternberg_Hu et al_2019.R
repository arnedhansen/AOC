## AOC Sternberg Effect Size (WM Load on Alpha) from Hu et al, 2019

# T-values from Hu et al., 2019 (UF data; digits)
# "During the retention period, alpha power increased as a function of memory load"
# Load 5 > Load 1, t(19) = −2.67, p = .015 -> d = 0.19
# Load 3 > Load 1, t(19) = −2.21, p = .039
# Load 5 > Load 3, t(19) = −2.15, p = .048
t_values <- c(-2.67, -2.21, -2.15)
n <- 20

# Function to calculate Cohen's d from paired t-values (adapted from Baguley, 2012)
cohens_d <- function(t, n) {
  d <- t / sqrt(n) * sqrt(2 / n)
  return(d)
}

# Calculate Cohen's d for each t-value
cohens_d_values <- sapply(t_values, cohens_d, n = n)

# Print the results
cohens_d_values
