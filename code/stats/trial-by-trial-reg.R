
#  This R function runs trial-by-trial logistic / mixed-model regressions for antisac response outcome (correct/error) and latency with 
#  preparatory power changes in the beta / theta band with subjects as random effect

#  Copyright (C) April 2020, last modified April 2022
#   J.Waldthaler, A. Sperlich, D. Pedrosa
#   University Hospital of Gießen and Marburg
#
#   This software may be used, copied, or redistributed as long as it is
#   not sold and this copyright notice is reproduced on each copy made.
#   This routine is provided as is without any express or implied
#   warranties whatsoever.

# load required packages
library(tidyverse) # data wrangling and visualization
library(sjPlot)    # to visualizing mixed-effects models
library(effects)   # to visualizing mixed-effects models
library(lme4)      # "golden standard" for mixed-effects modelling in R (no p-values)
library(lmerTest)  # p-values for MEMs based on the Satterthwaite approximation
library(report)    # mainly for an "report" function
library(emmeans)   # post-hoc analysis
library(knitr)     # beautifying tables
library(sjstats)   # ICC - intraclass-correlation coefficient
library(caret)     # ML, model comparison & utility functions
library(ggeffects)

# select folder and files 
setwd("~/Documents")
theta_filename <- paste("theta_log_reg.csv")
beta_filename <- paste("beta_log_reg.csv")
theta_lat_filename <- paste("theta_log_reg_lat.csv")
beta_lat_filename <- paste("beta_log_reg_lat.csv")
theta <- read.csv(theta_filename)
beta <- read.csv(beta_filename)
theta_lat <-read.csv(theta_lat_filename)
beta_lat <-read.csv(beta_lat_filename)

# set off DBS condition to be the reference
theta <- within(theta, condition <- relevel(condition, ref = "off"))
beta <- within(beta, condition <- relevel(condition, ref = "off"))

# ----------------first for THETA-------------------------
# ---logistic regression for trial outcome (correct / error)---

lmer.bin.max = glmer(response ~ theta * condition + (1 | subj), data=theta, family = binomial)
# description of the maximal model (used for reporting estimates and standard errors)
summary(lmer.bin.max)
tab_model(lmer.bin.max,   show.aic = T)

# plot
pr <- ggpredict(lmer.bin.max, c("theta", "condition"))
plot(pr)

# evaluate the contribution of the INTERACTION to the overall fit
# remove interaction term from model and compare to max model
lmer.bin.int = glmer(response ~ theta + condition + (1 | subj), data=theta, family = binomial)
anova(lmer.bin.max,lmer.bin.int)

# evaluate the contribution of the main effect of condition to the overall fit
lmer.bin.con = glmer(response ~ theta  + (1 | subj), data=theta, family = binomial)
anova(lmer.bin.int,lmer.bin.con)

# to compare pairwise, run original model again with new reference
# and correct p values (here: 0.25)
theta <- within(theta, condition <- relevel(condition, ref = "130"))
lmer.bin.max2 = glmer(response ~ theta * condition + (1 | subj), data=theta, family = binomial)
summary(lmer.bin.max2)
tab_model(lmer.bin.max2,   show.aic = T)

# evaluate the contribution of the main effect of cortical power to the overall fit
lmer.bin.eta = glmer(response ~ condition  + (1 | subj), data=theta, family = binomial)
anova(lmer.bin.int,lmer.bin.eta)

#---mixed linear model for latency---

lmer.bin.max = lmer(latency ~ theta  * condition + (1 | subj), data=theta_lat)
# description of the maximal model (used for reporting estimates and standard errors)
summary(lmer.bin.max)
tab_model(lmer.bin.max,   show.aic = T)

# plot
pr <- ggpredict(lmer.bin.max, c("theta", "condition"))
plot(pr, use.theme=FALSE, colors=c("black", "red", "blue"), show.title = FALSE, 
     add.data=TRUE, limit.range=FALSE)

# evaluate the contribution of the INTERACTION to the overall fit
# remove interaction term from model and compare to max model
lmer.bin.int = lmer(latency ~ theta + condition + (1 | subj), data=theta_lat)
anova(lmer.bin.max,lmer.bin.int)

# evaluate the contribution of the main effect of condition to the overall fit
lmer.bin.con = lmer(latency ~ theta  + (1 | subj), data=theta_lat)
anova(lmer.bin.int,lmer.bin.con)

# to compare pairwise, run original model again with new reference
# and correct p values (here: 0.25) 
theta_lat <- within(theta_lat, condition <- relevel(condition, ref = "130"))
lmer.bin.max2 = lmer(latency ~ theta * condition + (1 | subj), data=theta_lat)
summary(lmer.bin.max2)
tab_model(lmer.bin.max2,   show.aic = T)

# evaluate the contribution of the main effect of cortical power to the overall fit
theta_lat <- within(theta_lat, condition <- relevel(condition, ref = "off"))
lmer.bin.eta = lmer(latency ~ condition  + (1 | subj), data=theta_lat)
anova(lmer.bin.int,lmer.bin.eta)


# --------------------same for BETA-------------------------
# ---logistic regression for trial outcome (correct / error)---

lmer.bin.max = glmer(response ~ beta * condition + (1 | subj), data=beta, family = binomial)
# description of the maximal model (used for reporting estimates and standard errors)
summary(lmer.bin.max)
tab_model(lmer.bin.max,   show.aic = T)

# plot
pr <- ggpredict(lmer.bin.max, c("beta", "condition"))
plot(pr)

# evaluate the contribution of the INTERACTION to the overall fit
# remove interaction term from model and compare to max model
lmer.bin.int = glmer(response ~ beta + condition + (1 | subj), data=beta, family = binomial)
anova(lmer.bin.max,lmer.bin.int)

# evaluate the contribution of the main effect of condition to the overall fit
lmer.bin.con = glmer(response ~ beta  + (1 | subj), data=beta, family = binomial)
anova(lmer.bin.int,lmer.bin.con)

# to compare pairwise, run original model again with new reference
# and correct p values (here: 0.25)
beta <- within(beta, condition <- relevel(condition, ref = "130"))
lmer.bin.max2 = glmer(response ~ beta * condition + (1 | subj), data=beta, family = binomial)
summary(lmer.bin.max2)
tab_model(lmer.bin.max2,   show.aic = T)

# evaluate the contribution of the main effect of cortical power to the overall fit
lmer.bin.eta = glmer(response ~ condition  + (1 | subj), data=beta, family = binomial)
anova(lmer.bin.int,lmer.bin.eta)

#---mixed linear model for latency---

lmer.bin.max = lmer(latency ~ beta  * condition + (1 | subj), data=beta_lat)
# description of the maximal model (used for reporting estimates and standard errors)
summary(lmer.bin.max)
tab_model(lmer.bin.max,   show.aic = T)

# plot
pr <- ggpredict(lmer.bin.max, c("theta", "condition"))
plot(pr, use.theme=FALSE, colors=c("black", "red", "blue"), show.title = FALSE, 
     add.data=TRUE, limit.range=FALSE)

# evaluate the contribution of the INTERACTION to the overall fit
# remove interaction term from model and compare to max model
lmer.bin.int = lmer(latency ~ beta + condition + (1 | subj), data=beta_lat)
anova(lmer.bin.max,lmer.bin.int)

# evaluate the contribution of the main effect of condition to the overall fit
lmer.bin.con = lmer(latency ~ beta  + (1 | subj), data=beta_lat)
anova(lmer.bin.int,lmer.bin.con)

# to compare pairwise, run original model again with new reference
# and correct p values (here: 0.25) 
beta_lat <- within(beta_lat, condition <- relevel(condition, ref = "130"))
lmer.bin.max2 = lmer(latency ~ beta * condition + (1 | subj), data=beta_lat)
summary(lmer.bin.max2)
tab_model(lmer.bin.max2,   show.aic = T)

# evaluate the contribution of the main effect of cortical power to the overall fit
beta_lat <- within(theta_lat, condition <- relevel(condition, ref = "off"))
lmer.bin.eta = lmer(latency ~ condition  + (1 | subj), data=beta_lat)
anova(lmer.bin.int,lmer.bin.eta)
