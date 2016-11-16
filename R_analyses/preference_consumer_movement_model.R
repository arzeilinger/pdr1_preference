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





#########################################################################
#### Calculating Confidence Intervals using Quadratic approximation
#########################################################################
# Method based on description in Bolker (2008) pgs. 196 - 201.
# Method can only be used if MLE is close to global maximum.
# If parameter estimates are on or near an inequality constraint (e.g., ~0.0001),
# then the jackkife method must be used to estimate variance.

# Calculate SE and correlation matrices for each model
op.list <- list(op.fixed, op.p.choice, op.mu.choice, op.choice) # Place all optimx outputs in a list
se.list <- list(length(4)) # Create empty list for standard errors
cormat.list <- list(length(4)) # Create empty list for correlation matrices
for(i in 1:length(op.list)){
  mod.i <- op.list[[i]]
  atr.i <- attr(mod.i, "details") # Extract attributes from each model output
  hess.i <- atr.i[1,3][[1]] # Extract Hessian matrix from each output
  vcov.i <- solve(hess.i) # Calculate the Variance-Covariance Matrix
  cormat.list[[i]] <- cov2cor(vcov.i) # Calculate the the Corelation Matrix
  se.i <- sqrt(diag(vcov.i)) # Extract Variances and calculate standard errors
  se.list[[i]] <- se.i
}
# Extract variances for each model
se.fixed <- se.list[[1]]
se.p.choice <- se.list[[2]]
se.mu.choice <- se.list[[3]]
se.choice <- se.list[[4]]
# Create data frame of variances, with rows for each model and columns for each parameter
se.tab <- data.frame(rbind(c(se.fixed[1], se.fixed[1], se.fixed[2], se.fixed[2]),
                           c(se.p.choice[1], se.p.choice[2], se.p.choice[3], se.p.choice[3]),
                           c(se.mu.choice[1], se.mu.choice[1], se.mu.choice[2], se.mu.choice[3]),
                           c(se.choice)))
names(se.tab) <- c("se.p1", "se.p2", "se.mu1", "se.mu2") 
mle.tab <- cbind(mle.tab, se.tab) # Combine with parameter estimates and AIC results




################################################################################
#### Model-averaged parameter estimates and variances ##########################
################################################################################
# Note: Averaging only models with some information, i.e., those with dAICc < 7

# All models likelihoods and weights
amd <- mle.tab$dAIC
alv <- exp(-amd/2)
tml <- sum(alv) # Total marginal likelihood
awv <- alv/tml

# Good model likelihoods and weights
gmd <- mle.tab[mle.tab$dAIC < 7,] # select only good models
gmd$gml <- exp(-gmd$dAIC/2) # Relative or marginal likelihoods for models with some information (within)
gmd$wts <- gmd$gml/tml # Calculate relative weights for each model

# Calculate averaged parameter estimates and unconditional variances for each parameter
# p1
p1av <- sum(gmd$wts*gmd$p1)
p1var <- gmd$se.p1^2
p1var.av <- sum(gmd$wts*(p1var + (gmd$p1 - p1av)^2))
# p2
p2av <- sum(gmd$wts*gmd$p2)
p2var <- gmd$se.p2^2
p2var.av <- sum(gmd$wts*(p2var + (gmd$p2 - p2av)^2))
# mu1
mu1av <- sum(gmd$wts*gmd$mu1)
mu1var <- gmd$se.mu1^2
mu1var.av <- sum(gmd$wts*(mu1var + (gmd$mu1 - mu1av)^2))
# mu2
mu2av <- sum(gmd$wts*gmd$mu2)
mu2var <- gmd$se.mu2^2
mu2var.av <- sum(gmd$wts*(mu2var + (gmd$mu2 - mu2av)^2))

# Table of model-averaged parameter estimates and variances
resav <- data.frame(params = c("p1", "p2", "mu1", "mu2"),
                    estimates = c(p1av, p2av, mu1av, mu2av), 
                     vars = c(p1var.av, p2var.av, mu1var.av, mu2var.av))
resav$se <- sqrt(resav$vars)




#####################################################################################
#### Equilibrial probabilities for stink bug data
#####################################################################################

## Functions to calculate equilibrium probabilities, P1* and P2* from parameter estimates
# Equilibrium equations from Equation A12 in Appendix A of Zeilinger et al. (2014)
P1eq.func <- function(p1, p2, mu1, mu2){
  P1eq <- (p1*mu2)/(mu1*p2 + mu2*p1 + mu1*mu2)
  return(P1eq)
}
P2eq.func <- function(p1, p2, mu1, mu2){
  P2eq <- (p2*mu1)/(mu1*p2 + mu2*p1 + mu1*mu2)
  return(P2eq)
}

p1hat <- resav$estimates[1]; p2hat <- resav$estimates[2]
mu1hat <- resav$estimates[3]; mu2hat <- resav$estimates[4]

P1.star <- P1eq.func(p1hat, p2hat, mu1hat, mu2hat) # P1 equilibrium
P2.star <- P2eq.func(p1hat, p2hat, mu1hat, mu2hat) # P2 equilibrium



