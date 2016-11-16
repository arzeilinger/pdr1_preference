#### ANALYSIS OF BGSS PREFERENCE ON PDR1 RESISTANT AND SUSCEPTIBLE INFECTED PLANTS
#### USING CONSUMER MOVEMENT MODEL
# Code includes maximum likelihood estimation of 4 model variants, model selection,
# estimation of variance using normal approximation method
# model averaging, and calculation of equilibrial probabilities, Pj*

rm(list = ls())
# libraries
my.packages <- c("xlsx", "tidyr", "dplyr", "data.table", "ggplot2", "lattice", "optimx", "bbmle", "numDeriv")
lapply(my.packages, require, character.only = TRUE)

## Input functions for P1 and P2 equations, 
## negative log likelihood (NLL) functions, 
## gradient functions, and AIC function
source("R_functions/consumer_movement_model_source_functions.R")


###################################################################################################
#### Importing data 

prefdata <- read.xlsx("data/pdr1_preference_data.xlsx", sheetName = "data")
str(prefdata)

# Separate leaf and symptom data from preference count data
leafdata <- prefdata[,13:ncol(prefdata)]
prefdata <- prefdata[,1:12]
prefdata$cage <- paste(prefdata$trt, prefdata$rep, sep = "")
head(prefdata)

#### Calculate total counts among cages for each week, genotype, and time point
# n1 = source plant (Xylella infected)
# n2 = test plant
# n3 = neutral space

cmmData <- prefdata %>% group_by(week, trt, time_from_start_hr) %>% 
  summarise(n1 = sum(source_plant),
            n2 = sum(test_plant),
            n3 = sum(neutral_space)) %>%
  as.data.frame()
cmmData$t <- cmmData$time_from_start_hr
cmmData$N <- with(cmmData, n1 + n2 + n3)
cmmData$week.trt <- factor(paste(cmmData$week, cmmData$trt, sep = ""))

r12Data <- dplyr::filter(cmmData, week == 12, trt == "R") 
s3Data <- dplyr::filter(cmmData, week == 3, trt == "S")
s12Data <- dplyr::filter(cmmData, week == 12, trt == "S")
s8Data <- dplyr::filter(cmmData, week == 8, trt == "S")
r8Data <- dplyr::filter(cmmData, week == 8, trt == "R")


#### Fit models to each of the trt-week combinations
modelFits <- lapply(levels(cmmData$week.trt), function(x) optimizeCMM(dat = cmmData[cmmData$week.trt == x,], upperConstraint = 64, aiccN = 8))

#### Calculate variance-covariance and correlation matrices for each trt-week combination
matrices <- lapply(modelFits, getParCorrelations)

#### Extract and organize parameter estimates from all good models (dAICc <= 7)
paramResults <- lapply(modelFits, mleTable)

#### Average results for all good models
paramAverage <- lapply(paramResults, function(x) averageModels(x, dAIC.threshold = 7))

# Add week.trt combination to each element of the list and combine into one data.frame
for(i in 1:length(paramAverage)){
  paramAverage[[i]]$week.trt <- levels(cmmData$week.trt)[i]
} 
paramData <- rbindlist(paramAverage) %>% as.data.frame()
paramData$estimate <- factor2numeric(paramData$estimate)
paramData$variance <- factor2numeric(paramData$variance)

# Calculate SE and 95% CI
paramData$se <- sqrt(paramData$variance)
paramData$cil <- with(paramData, estimate - 1.96*se)
paramData$ciu <- with(paramData, estimate + 1.96*se)

#########################################################################
#### Calculating Confidence Intervals using Quadratic approximation
#########################################################################
# Method based on description in Bolker (2008) pgs. 196 - 201.
# Method can only be used if MLE is close to global maximum.
# If parameter estimates are on or near an inequality constraint (e.g., ~0.0001),
# then the jackkife method must be used to estimate variance.




#####################################################################################
#### Equilibrial probabilities for stink bug data
#####################################################################################

## Functions to calculate equilibrium probabilities, P1* and P2* from parameter estimates
# Equilibrium equations from Equation A12 in Appendix A of Zeilinger et al. (2014)
P1eq.func <- function(params){
  p1 <- params[1]
  p2 <- params[2]
  mu1 <- params[3]
  mu2 <- params[4]
  P1eq <- (p1*mu2)/(mu1*p2 + mu2*p1 + mu1*mu2)
  return(P1eq)
}
P2eq.func <- function(params){
  p1 <- params[1]
  p2 <- params[2]
  mu1 <- params[3]
  mu2 <- params[4]
  P2eq <- (p2*mu1)/(mu1*p2 + mu2*p1 + mu1*mu2)
  return(P2eq)
}

# S12
params12S <- paramData[paramData$week.trt == "12S", "estimate"]
P1eq.func(params12S) # P1 equilibrium
P2eq.func(params12S) # P2 equilibrium

# R12
params12R <- paramData[paramData$week.trt == "12R", "estimate"]
P1eq.func(params12R) # P1 equilibrium
P2eq.func(params12R) # P2 equilibrium



