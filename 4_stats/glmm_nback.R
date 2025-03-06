# GLMM Stats for AOC N-back

install.packages('lme4')
install.packages('nlme')
install.packages('emmeans')
install.packages("ggplot2")
install.packages("pbkrtest")
install.packages("lmerTest")
install.packages("sjPlot")
install.packages("writexl")
install.packages("drop1")
install.packages("car")
library(lme4)
library(nlme)
library(emmeans)
library(ggplot2)
library(pbkrtest)
library(lmerTest)
library(sjPlot)
library("writexl")
library("lm.beta") 
library(drop1)
library(car)

# Set WD
setwd('/Users/Arne/Documents/GitHub/AOC/4_stats')

# no scientific notation, exact p-values
options(scipen=999) # default = 0
options(scipen=0)
# printing smaller p-values than e-16
# https://stackoverflow.com/questions/6970705/why-cant-i-get-a-p-value-smaller-than-2-2e-16

# load data 
dat = read.csv('/Volumes/methlab/Students/Arne/AOC/data/features/merged_data_nback.csv')

# make sure, the vars are factors
dat$ID = as.factor(dat$ID)
dat$Condition = as.factor(dat$Condition)

######################
##### Behavioral #####
######################

# Model for Accuracy
glmm_accuracy <- lmer(dat$Accuracy ~ dat$Condition + (1|dat$ID), data = dat)
summary(glmm_accuracy)

# Anova
Anova(glmm_accuracy, type = "II")

# Model for Reaction Time
glmm_rt <- lmer(dat$ReactionTime ~ dat$Condition + (1|dat$ID), data = dat)
summary(glmm_rt)

# Anova
Anova(glmm_rt, type = "II")

# Extract coefficients
coeff_accuracy <- summary(glmm_accuracy)[["coefficients"]]
coeff_accuracy
coeff_rt <- summary(glmm_rt)[["coefficients"]]
coeff_rt

# Display models

######################
######## Gaze ########
######################

# Model
glmm_gaze <- lmer(dat$GazeDeviation ~ dat$Condition + (1|dat$ID), data = dat)
summary(glmm_gaze)

# Anova
Anova(glmm_gaze, type = "II")

# Extract coefficients
coeff_gaze <- summary(glmm_gaze)[["coefficients"]]
coeff_gaze

# Display models

######################
####### Alpha ########
######################

# Model
glmm_alpha <- lmer(dat$AlphaPower ~ dat$GazeDeviation * dat$Condition + (1|dat$ID), data = dat)
summary(glmm_alpha)

# Anova
Anova(glmm_alpha, type = "III")

# Dropping terms and finalizing model if needed
s_alpha <- drop1(glmm_alpha)
glmm_alpha_final <- get_model(s_alpha)
summary(glmm_alpha_final)
#Anova(glmm_alpha_final, type = "II")

# Extract coefficients
coeff_alpha <- summary(glmm_alpha_final)[["coefficients"]]

# Display models
tab_model(glmm_alpha_final)