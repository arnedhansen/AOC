#!/usr/bin/env Rscript
# AOC sternberg mixed models (trial-level)
# load, prep, models — separate from condition-level confirmatory script
# RE: (1|Subject) + (1|Trial); IQR: 3× within Subject×Load

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(car)
})

options(scipen = 999)

confirmatory_script_dir <- {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("--file=", "", args[grep("^--file=", args)])
  if (length(file_arg) > 0) {
    dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE))
  } else {
    getwd()
  }
}
source(file.path(confirmatory_script_dir, "..", "AOC_stats_glmm_helpers.R"))

######################
###### Load Data #####
######################

dat <- read.csv("SET_PATH_TO_AOC_MERGED_DATA_STERNBERG_TRIALS_CSV", stringsAsFactors = FALSE)


######################
###### Prep Data #####
######################

dat$ID <- as.factor(dat$ID)
dat$Subject <- dat$ID
dat$Trial <- as.factor(dat$Trial)
dat$Load <- factor(
  c(`2` = "WM load 2", `4` = "WM load 4", `6` = "WM load 6")[as.character(as.numeric(dat$Condition))],
  levels = c("WM load 2", "WM load 4", "WM load 6"),
  ordered = FALSE
)
dat <- dat[!is.na(dat$Load), , drop = FALSE]

dat$Accuracy <- as.numeric(dat$Accuracy)
dat$ReactionTime <- as.numeric(dat$ReactionTime)
dat$ERSD <- as.numeric(dat$ERSD_late)

# Fully standardized ERS/ERD for co-variation models (predictor + outcome).
ersd_z <- grand_mean_zscore(dat$ERSD)
dat$ERSD_z <- ersd_z$z
dat$GazeDeviationBL <- as.numeric(dat$GazeDeviationLateBL)
dat$MSRateBL <- as.numeric(dat$MSRateLateBL)

#########################################
###### Behavior and ERS/ERD ~ Load #######
#########################################

######################
###### Accuracy ######
######################

m_accuracy <- lmer(Accuracy ~ Load + (1 | Subject) + (1 | Trial), data = dat, REML = FALSE)
summary(m_accuracy)
Anova(m_accuracy, type = "II")

############################
###### Reaction Time #######
############################

m_rt <- lmer(ReactionTime ~ Load + (1 | Subject) + (1 | Trial), data = dat, REML = FALSE)
summary(m_rt)
Anova(m_rt, type = "II")

######################
######## ERSD ########
######################

m_ersd <- lmer(ERSD ~ Load + (1 | Subject) + (1 | Trial), data = dat, REML = FALSE)
summary(m_ersd)
Anova(m_ersd, type = "II")

#########################################################
###### Trial-level baselined gaze metrics only ##########
#########################################################

##########################
###### GazeDeviationBL ###
##########################

m_gaze_bl <- lmer(GazeDeviationBL ~ Load + (1 | Subject) + (1 | Trial), data = dat, REML = FALSE)
summary(m_gaze_bl)
Anova(m_gaze_bl, type = "II")

######################
###### MSRateBL ######
######################

m_ms_bl <- lmer(MSRateBL ~ Load + (1 | Subject) + (1 | Trial), data = dat, REML = FALSE)
summary(m_ms_bl)
Anova(m_ms_bl, type = "II")

######################################
###### ERSD x GazeDeviationBL ########
######################################

dat <- add_gaze_z_column(dat, "GazeDeviationBL")
m_ersd_gaze_bl_full <- lmer(ERSD_z ~ GazeDeviationBL_z * Load + (1 | Subject) + (1 | Trial), data = dat, REML = FALSE)
m_ersd_gaze_bl_add <- lmer(ERSD_z ~ GazeDeviationBL_z + Load + (1 | Subject) + (1 | Trial), data = dat, REML = FALSE)
lrt_ersd_gaze_bl <- anova(m_ersd_gaze_bl_add, m_ersd_gaze_bl_full)
if (is.finite(lrt_ersd_gaze_bl$`Pr(>Chisq)`[2]) && lrt_ersd_gaze_bl$`Pr(>Chisq)`[2] < 0.05) {
  m_ersd_gaze_bl_final <- m_ersd_gaze_bl_full
} else {
  m_ersd_gaze_bl_final <- m_ersd_gaze_bl_add
}
summary(m_ersd_gaze_bl_final)
Anova(m_ersd_gaze_bl_final, type = "III")

######################################
######## ERSD x MSRateBL #############
######################################

dat <- add_gaze_z_column(dat, "MSRateBL")
m_ersd_ms_bl_full <- lmer(ERSD_z ~ MSRateBL_z * Load + (1 | Subject) + (1 | Trial), data = dat, REML = FALSE)
m_ersd_ms_bl_add <- lmer(ERSD_z ~ MSRateBL_z + Load + (1 | Subject) + (1 | Trial), data = dat, REML = FALSE)
lrt_ersd_ms_bl <- anova(m_ersd_ms_bl_add, m_ersd_ms_bl_full)
if (is.finite(lrt_ersd_ms_bl$`Pr(>Chisq)`[2]) && lrt_ersd_ms_bl$`Pr(>Chisq)`[2] < 0.05) {
  m_ersd_ms_bl_final <- m_ersd_ms_bl_full
} else {
  m_ersd_ms_bl_final <- m_ersd_ms_bl_add
}
summary(m_ersd_ms_bl_final)
Anova(m_ersd_ms_bl_final, type = "III")
