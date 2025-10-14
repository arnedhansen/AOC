# GLMM Stats for AOC nback

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
source("/Users/Arne/Documents/GitHub/functions/export_model_table.R")


# Set WD
setwd('/Users/Arne/Documents/GitHub/AOC/4_stats')

# no scientific notation, exact p-values
options(scipen=999) # default = 0
options(scipen=0)
# printing smaller p-values than e-16
# https://stackoverflow.com/questions/6970705/why-cant-i-get-a-p-value-smaller-than-2-2e-16

# load data 
dat = read.csv('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/merged_data_nback.csv')
dat = na.omit(dat)

# make sure, the vars are factors
dat$ID = as.factor(dat$ID)
dat$Condition = as.factor(dat$Condition)

######################
##### Behavioral #####
######################

# Model for Accuracy
glmm_alphaMS <- lmer(dat$AlphaPower ~ dat$MSRate + (1|dat$ID), data = dat)
summary(glmm_alphaMS)

# Anova
Anova(glmm_alphaMS, type = "II")

# Extract coefficients
coeff_accuracy <- summary(glmm_alphaMS)[["coefficients"]]
coeff_accuracy

# Save model overview
export_model_table(glmm_alphaMS, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/stats/glmm/AOC_nback_glmm_TEST_alphaMS.docx')

######################
######## XXXX ########
######################

