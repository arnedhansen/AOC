



setwd('V:\hswanb\01-InfCry\02-IC_EyeTracking\03-scripts\01-analysis\02-behavioural')

library(lme4)
library(afex)
library(nlme)
library(psych)
library(lsmeans)
library(ggplot2)
library(stringi)
library(stringr)
library(plyr)
library(ISLR)
library(multcomp)
library(sjPlot)
library(purrr)
#install.packages("ggeffects")
library(ggeffects)
library(wesanderson)
library(ggsci)
library(svglite)
library(MMeM)
library(MCMCglmm)
library(ordinal)
library(bayesplot)
library(cv)
library(lmeresampler)

dataStr <- 'V:/hswanb/01-InfCry/02-IC_EyeTracking/02-analysis/01-ET/bstrap/ET_BootStrap_ready.csv'

#Turn data types into factors for the modelling.
data  <- read.csv(dataStr)
data$subjID = as.factor(data$subjID)
data$trials = as.factor(data$trials)
data$intensity = as.factor(data$intensity)
data$contour = as.factor(data$contour)
data$sourceID = as.factor(data$sourceID)
data$sex = as.factor(data$sex)
data$parental = as.factor(data$parental)

data$ordSalient = factor(data$salient,ordered = TRUE)
data$ordAversive = factor(data$aversive,ordered = TRUE)
data$ordinalPosition = factor(data$ordinalPosition,ordered = TRUE)



## LMM Without Demographic Information
aversiveTallyModel <- lmer(tallyWin ~   mps  * contours * aversion + (1|subjID), data,  REML = FALSE)
p_val_AvTally <- bootstrap_pvals(model = aversiveTallyModel, type = "parametric", B = 1000)


salienceTallyModel <- lmer(tallyWin ~   mps  * contours * salience + (1|subjID), data,  REML = FALSE)
p_val_SaTally <- bootstrap_pvals(model = salienceTallyModel, type = "parametric", B = 1000)

























tallyContourOnly <- lmer(tallyWin ~   contours * sourceLang + (1|subjID), data,  REML = FALSE)
p_val_MPS<- bootstrap_pvals(model = tallyContourOnly, type = "parametric", B = 1000)





tallyContourOnly <- lmer(tallyWin ~   contours * sourceLang + (1|subjID), data,  REML = FALSE)
p_val_MPS<- bootstrap_pvals(model = tallyContourOnly, type = "parametric", B = 1000)




aversiveModel <- lmer(aversion ~   mps  * contours + (1|subjID), data,  REML = FALSE)






p_values <- bootstrap_pvals(model = aversiveModel, type = "parametric", B = 1000)




salienceModel <- lmer(salience ~   mps  * contours + (1|subjID), data,  REML = FALSE)

salienceSexModel <- lmer(salience ~   mps  * contours * sex + (1|subjID), data,  REML = FALSE)
aversiveSexModel <- lmer(aversion ~   mps  * contours * sex + (1|subjID), data,  REML = FALSE)





p_valuesSal <- bootstrap_pvals(model = salienceModel, type = "parametric", B = 1000)

p_valuesAvSex <- bootstrap_pvals(model = aversiveSexModel, type = "parametric", B = 1000)
p_valuesSalSex <- bootstrap_pvals(model = salienceSexModel, type = "parametric", B = 1000)















cfData <-groupdata2::fold(data)
cross_validate(data = cfData,family = 'gaussian', formulas = 'aversive ~ stimMPS  * contour + (1|subjID)+(1|ordinalPosition)')





aversiveModel <- lmer(aversive ~   stimMPS  * contour + (1|subjID)+(1|ordinalPosition), data,  REML = FALSE)


cross_validate(data = data,family = 'gaussian', formulas = stimMPS  * contour + (1|subjID)+(1|ordinalPosition))

cross_validate(data = data, family = 'gaussian', 
               
               avSummary <- summary(aversiveModel)
               
               salienceModel <- lmer(salient ~   stimMPS  * contour + (1|subjID)+(1|ordinalPosition), data, REML = FALSE)
               salSummary <- summary(salienceModel)
               res<-resid(salienceModel,na.action=na.omit) 
               