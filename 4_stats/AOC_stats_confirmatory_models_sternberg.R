#!/usr/bin/env Rscript
# AOC sternberg mixed models
# load, prep, models

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
source(file.path(confirmatory_script_dir, "AOC_stats_glmm_helpers.R"))

######################
###### Load Data #####
######################

dat <- read.csv("SET_PATH_TO_AOC_MERGED_DATA_STERNBERG_CSV", stringsAsFactors = FALSE)

######################
###### Prep Data #####
######################

dat$ID <- as.factor(dat$ID)
dat$Subject <- dat$ID
dat$Load <- factor(
  c(`2` = "WM load 2", `4` = "WM load 4", `6` = "WM load 6")[as.character(as.numeric(dat$Condition))],
  levels = c("WM load 2", "WM load 4", "WM load 6"),
  ordered = FALSE
)
dat <- dat[!is.na(dat$Load), , drop = FALSE]

dat$Accuracy <- as.numeric(dat$Accuracy)
dat$ReactionTime <- as.numeric(dat$ReactionTime)
dat$GazeDeviation <- as.numeric(dat$GazeDeviation)
dat$MSRate <- as.numeric(dat$MSRate)
dat$ERSD <- as.numeric(dat$ERSD_late)
dat$GazeDeviationBL <- as.numeric(dat$GazeDeviationLateBL)
dat$MSRateBL <- as.numeric(dat$MSRateLateBL)

#########################################
###### Confirmatory Unbaselined Gaze ####
#########################################

######################
###### Accuracy ######
######################

m_accuracy <- lmer(Accuracy ~ Load + (1 | Subject), data = dat, REML = FALSE)
summary(m_accuracy)
Anova(m_accuracy, type = "II")

############################
###### Reaction Time #######
############################

m_rt <- lmer(ReactionTime ~ Load + (1 | Subject), data = dat, REML = FALSE)
summary(m_rt)
Anova(m_rt, type = "II")

######################
######## Gaze ########
######################

m_gaze <- lmer(GazeDeviation ~ Load + (1 | Subject), data = dat, REML = FALSE)
summary(m_gaze)
Anova(m_gaze, type = "II")

######################
###### MS Rate #######
######################

m_ms <- lmer(MSRate ~ Load + (1 | Subject), data = dat, REML = FALSE)
summary(m_ms)
Anova(m_ms, type = "II")

######################
######## ERSD ########
######################

m_ersd <- lmer(ERSD ~ Load + (1 | Subject), data = dat, REML = FALSE)
summary(m_ersd)
Anova(m_ersd, type = "II")

######################################
###### ERSD x GazeDeviation ##########
######################################

dat <- add_gaze_z_column(dat, "GazeDeviation")
m_ersd_gaze_full <- lmer(ERSD ~ GazeDeviation_z * Load + (1 | Subject), data = dat, REML = FALSE)
m_ersd_gaze_add <- lmer(ERSD ~ GazeDeviation_z + Load + (1 | Subject), data = dat, REML = FALSE)
lrt_ersd_gaze <- anova(m_ersd_gaze_add, m_ersd_gaze_full)
if (is.finite(lrt_ersd_gaze$`Pr(>Chisq)`[2]) && lrt_ersd_gaze$`Pr(>Chisq)`[2] < 0.05) {
  m_ersd_gaze_final <- m_ersd_gaze_full
} else {
  m_ersd_gaze_final <- m_ersd_gaze_add
}
summary(m_ersd_gaze_final)
Anova(m_ersd_gaze_final, type = "III")

######################################
######## ERSD x MSRate ###############
######################################

dat <- add_gaze_z_column(dat, "MSRate")
m_ersd_ms_full <- lmer(ERSD ~ MSRate_z * Load + (1 | Subject), data = dat, REML = FALSE)
m_ersd_ms_add <- lmer(ERSD ~ MSRate_z + Load + (1 | Subject), data = dat, REML = FALSE)
lrt_ersd_ms <- anova(m_ersd_ms_add, m_ersd_ms_full)
if (is.finite(lrt_ersd_ms$`Pr(>Chisq)`[2]) && lrt_ersd_ms$`Pr(>Chisq)`[2] < 0.05) {
  m_ersd_ms_final <- m_ersd_ms_full
} else {
  m_ersd_ms_final <- m_ersd_ms_add
}
summary(m_ersd_ms_final)
Anova(m_ersd_ms_final, type = "III")

############################################
###### Exploratory Baselined Gaze Data #####
############################################

##########################
###### GazeDeviationBL ###
##########################

m_gaze_bl <- lmer(GazeDeviationBL ~ Load + (1 | Subject), data = dat, REML = FALSE)
summary(m_gaze_bl)
Anova(m_gaze_bl, type = "II")

######################
###### MSRateBL ######
######################

m_ms_bl <- lmer(MSRateBL ~ Load + (1 | Subject), data = dat, REML = FALSE)
summary(m_ms_bl)
Anova(m_ms_bl, type = "II")

######################################
###### ERSD x GazeDeviationBL ########
######################################

dat <- add_gaze_z_column(dat, "GazeDeviationBL")
m_ersd_gaze_bl_full <- lmer(ERSD ~ GazeDeviationBL_z * Load + (1 | Subject), data = dat, REML = FALSE)
m_ersd_gaze_bl_add <- lmer(ERSD ~ GazeDeviationBL_z + Load + (1 | Subject), data = dat, REML = FALSE)
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
m_ersd_ms_bl_full <- lmer(ERSD ~ MSRateBL_z * Load + (1 | Subject), data = dat, REML = FALSE)
m_ersd_ms_bl_add <- lmer(ERSD ~ MSRateBL_z + Load + (1 | Subject), data = dat, REML = FALSE)
lrt_ersd_ms_bl <- anova(m_ersd_ms_bl_add, m_ersd_ms_bl_full)
if (is.finite(lrt_ersd_ms_bl$`Pr(>Chisq)`[2]) && lrt_ersd_ms_bl$`Pr(>Chisq)`[2] < 0.05) {
  m_ersd_ms_bl_final <- m_ersd_ms_bl_full
} else {
  m_ersd_ms_bl_final <- m_ersd_ms_bl_add
}
summary(m_ersd_ms_bl_final)
Anova(m_ersd_ms_bl_final, type = "III")
