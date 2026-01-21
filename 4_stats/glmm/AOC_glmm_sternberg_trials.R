# AOC GLMM — Sternberg (Trial-Level)
# Fits trial-level GLMMs on merged_data_sternberg_trials (Accuracy, RT, gaze, EEG) with Condition + (1|ID) or (1+Condition|ID). Anova, emmeans, model tables.
#
# Key outputs:
#   Model summaries; Anova tables; (optional) exported tables

#install.packages('lme4')
#install.packages('nlme')
#install.packages('emmeans')
#install.packages("ggplot2")
#install.packages("pbkrtest")
#install.packages("lmerTest")
#install.packages("sjPlot")
#install.packages("writexl")
#install.packages("drop1")
#install.packages("car")
#install.packages("glmmTMB")
library(lme4)
library(nlme)
library(emmeans)
library(ggplot2)
library(pbkrtest)
library(lmerTest)
library(sjPlot)
library("writexl")
library("lm.beta") 
#library(drop1)
library(car)
library(glmmTMB)
source("/Users/Arne/Documents/GitHub/functions/export_model_table.R")

# Set WD
setwd('/Users/Arne/Documents/GitHub/AOC/4_stats')

# no scientific notation, exact p-values
options(scipen=999) # default = 0
options(scipen=0)
# printing smaller p-values than e-16
# https://stackoverflow.com/questions/6970705/why-cant-i-get-a-p-value-smaller-than-2-2e-16

# load data 
dat0 <- read.csv('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/merged_data_sternberg_trials.csv')
needed <- c("ID", "Trial", "Condition") 
dat <- dat0[complete.cases(dat0[, needed]), ]
#cat(sprintf("Lost %d rows (%.2f%%) due to missing ID/Trial/Condition\n",
#            nrow(dat0) - nrow(dat),
#            100 * (nrow(dat0) - nrow(dat)) / nrow(dat0)))

# make sure, the vars are factors
dat$ID = as.factor(dat$ID)
dat$Trial = as.factor(dat$Trial)
dat$Condition = as.factor(dat$Condition)

z-score within participant

# When you include (1 | ID) or (1 + Condition | ID) in your GLMM, 
# you’re already accounting for the fact that multiple trials belong 
# to the same participant. The model therefore implicitly handles 
# trial-level dependencies by nesting observations within subjects. 
# You don’t need to include Trial as a random effect unless trial 
# numbers repeat across subjects and have some shared meaning (for 
# example, the same item or stimulus shown to all participants

######################
##### Behavioral #####
######################

# Model for Accuracy (glmer for binary data)
glmm_accuracy <- glmer(Accuracy ~ Condition + (1 + Condition | ID),
                       data = dat,
                       family = binomial(link = "logit"),
                       control = glmerControl(optimizer = "bobyqa",
                                              optCtrl = list(maxfun = 2e5)))
summary(glmm_accuracy)

# Estimated marginal means on the probability scale
emm_acc <- emmeans(glmm_accuracy, ~ Condition, type = "response")
summary(emm_acc)        # gives estimated probabilities per condition
pairs(emm_acc, adjust = "tukey")   # pairwise comparisons (odds ratios and prob diffs)

# Anova
Anova(glmm_accuracy, type = "II")

# Model for Reaction Time
glmm_rt <- lmer(log(ReactionTime) ~ Condition + (1 + Condition | ID), data = dat) # log-transform RT
summary(glmm_rt)

emm_rt <- emmeans(glmm_rt, ~ Condition)
summary(emm_rt, type = "response")   # back-transformed mean RTs (ms)
pairs(emm_rt, adjust = "tukey", type = "response")

# Anova
Anova(glmm_rt, type = "II")

# Extract coefficients
coeff_accuracy <- summary(glmm_accuracy)[["coefficients"]]
coeff_accuracy
coeff_rt <- summary(glmm_rt)[["coefficients"]]
coeff_rt

# Save model overview
export_model_table(glmm_accuracy, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/stats/glmm/AOC_sternberg_glmm_acc_trials.docx')
export_model_table(glmm_rt, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/stats/glmm/AOC_sternberg_glmm_rt_trials.docx')

## -------------------------
## GazeDeviation (Early, Late)
## -------------------------
glmm_gaze_early <- lmer(GazeDeviationEarly ~ Condition +
                          (1 + Condition | ID),
                        data = dat, REML = TRUE,
                        control = lmerControl(optimizer = "bobyqa"))
summary(glmm_gaze_early)
Anova(glmm_gaze_early, type = "II")
export_model_table(glmm_gaze_early, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/stats/glmm/AOC_sternberg_glmm_gazeDevEarly_trials.docx')

glmm_gaze_late <- lmer(GazeDeviationLate ~ Condition +
                         (1 + Condition | ID),
                       data = dat, REML = TRUE,
                       control = lmerControl(optimizer = "bobyqa"))
summary(glmm_gaze_late)
Anova(glmm_gaze_late, type = "II")
export_model_table(glmm_gaze_late, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/stats/glmm/AOC_sternberg_glmm_gazeDevLate_trials.docx')

## -------------------------
## Microsaccade Rate (Early, Late) — not BL
## -------------------------

glmm_ms_early_gamma <- lmer(MSRateEarly ~ Condition + (1 + Condition | ID), data = dat)

export_model_table(glmm_ms_early_gamma, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/stats/glmm/AOC_sternberg_glmm_MSRateEarly_trials.docx')

glmm_ms_late_gamma <- lmer(MSRateLate ~ Condition + (1 + Condition | ID), data = dat)

export_model_table(glmm_ms_late_gamma, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/stats/glmm/AOC_sternberg_glmm_MSRateLate_trials.docx')

## -------------------------
## ScanPathLength (Early, Late)
## -------------------------
glmm_scan_early <- lmer(log1p(ScanPathLengthEarly) ~ Condition +
                          (1 + Condition | ID),
                        data = dat, REML = TRUE,
                        control = lmerControl(optimizer = "bobyqa"))
summary(glmm_scan_early)
export_model_table(glmm_scan_early, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/stats/glmm/AOC_sternberg_glmm_scanEarly_trials.docx')

glmm_scan_late <- lmer(log1p(ScanPathLengthLate) ~ Condition +
                         (1 + Condition | ID),
                       data = dat, REML = TRUE,
                       control = lmerControl(optimizer = "bobyqa"))
summary(glmm_scan_late)
export_model_table(glmm_scan_late, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/stats/glmm/AOC_sternberg_glmm_scanLate_trials.docx')

## -------------------------
## ScanPathLength (Full window 0–2 s)
## -------------------------

# Construct Full measures (assumes Early [0–1] and Late [1–2] are non-overlapping)
dat$ScanPathLengthFull   <- dat$ScanPathLengthEarly + dat$ScanPathLengthLate

# Main effect of Condition (Full)
glmm_scan_full <- lmer(ScanPathLengthFull ~ Condition +
                         (1 + Condition | ID),
                       data = dat, REML = TRUE,
                       control = lmerControl(optimizer = "bobyqa"))
summary(glmm_scan_full)
Anova(glmm_scan_full, type = "II")
export_model_table(glmm_scan_full, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/stats/glmm/AOC_sternberg_glmm_scanFull_trials.docx')

## -------------------------
## Alpha power: main effect of Condition
## -------------------------
glmm_alpha_early_cond <- lmer(AlphaPowerEarly ~ Condition +
                                (1 + Condition | ID),
                              data = dat, REML = TRUE,
                              control = lmerControl(optimizer = "bobyqa"))
summary(glmm_alpha_early_cond)
Anova(glmm_alpha_early_cond, type = "II")
export_model_table(glmm_alpha_early_cond, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/stats/glmm/AOC_sternberg_glmm_alphaEarly_cond_trials.docx')

glmm_alpha_late_cond <- lmer(AlphaPowerLate ~ Condition +
                               (1 + Condition | ID),
                             data = dat, REML = TRUE,
                             control = lmerControl(optimizer = "bobyqa"))
summary(glmm_alpha_late_cond)
Anova(glmm_alpha_late_cond, type = "II")
export_model_table(glmm_alpha_late_cond, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/stats/glmm/AOC_sternberg_glmm_alphaLate_cond_trials.docx')
