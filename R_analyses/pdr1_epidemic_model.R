#### Pierce's Disease SECI epidemic Model
#### Developed by Adam Zeilinger, based on Zeilinger and Daugherty (2014)

rm(list = ls())
# Load packages
my.packages <- c("tidyr", "dplyr", "data.table", "deSolve")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/pdr1_epidemic_model_functions.R")
source("R_functions/simulateData.R")
source("R_functions/factor2numeric.R")

##########################################################################################
#### Compiling parameter estimates
##########################################################################################
#### Loading data sets
# Rate parameters
rateParams <- read.csv("results/pdr1_cmm_rate_parameter_estimates.csv", header = TRUE)
rateParams

# Transmission parameters
transdata <- readRDS("output/pdr1_transmission_preference_dataset.rds")
str(transdata)

# Construct source plant infection status column
transdata$source.cfu.per.g <- factor2numeric(transdata$source.cfu.per.g)
transdata$source.infection <- transdata$source.cfu.per.g
transdata$source.infection[transdata$source.infection > 0] <- 1

# Transmission summary table
transSummary <- transdata %>% group_by(., week, trt) %>% 
  summarise(propInfected = sum(test.plant.infection, na.rm = TRUE)/sum(!is.na(test.plant.infection)),
            sepropI = sqrt(propInfected*(1 - propInfected)/8),
            meanPD = mean(pd_index, na.rm = TRUE),
            sePD = sd(pd_index, na.rm = TRUE)/sqrt(sum(!is.na(pd_index))),
            meancfu = mean(log10(source.cfu.per.g+1)),
            secfu = sd(log10(source.cfu.per.g+1))/sqrt(sum(!is.na(source.cfu.per.g))),
            sourcepInfected = sum(source.infection, na.rm = TRUE)/sum(!is.na(source.infection)))

# Empty list of length 2 to hold the simulated parameter values
parList <- vector("list", 2)

for(i in 1:length(levels(rateParams$trt))){
  trt.i <- levels(rateParams$trt)[i]
  
  #### Rate parameters
  # phi_pc = 8S p1
  # phi_pi = 12S p1
  # phi_muc = 8S mu1
  # phi_mui = 12S mu1
  # Multiply by 24 to get rate in day-1
  
  phi_pcBar <- rateParams %>% dplyr::filter(., parameter == "p1" & week == 8 & trt == trt.i) %>% dplyr::select(., estimate, se)*24
  phi_pcVec <- simulateData(phi_pcBar$estimate, phi_pcBar$se)
  
  phi_piBar <- rateParams %>% dplyr::filter(., parameter == "p1" & week == 12 & trt == trt.i) %>% dplyr::select(., estimate, se)*24
  phi_piVec <- simulateData(phi_piBar$estimate, phi_piBar$se)
  
  phi_mucBar <- rateParams %>% dplyr::filter(., parameter == "mu1" & week == 8 & trt == trt.i) %>% dplyr::select(., estimate, se)*24
  phi_mucVec <- simulateData(phi_mucBar$estimate, phi_mucBar$se)
  
  phi_muiBar <- rateParams %>% dplyr::filter(., parameter == "mu1" & week == 12 & trt == trt.i) %>% dplyr::select(., estimate, se)*24
  phi_muiVec <- simulateData(phi_muiBar$estimate, phi_muiBar$se)
  
  #### Transmission parameters
  # beta_c ~ P(trans) 8S
  # beta_i ~ P(trans) 12S
  # Divide P(trans) by 8 to give the inoculation per vector
  # Assume that inoculation is 60% of acquisition, because I don't have acquisition data yet
  # alpha_c = beta_c/0.6
  # alpha_i = beta_i/0.6
  beta_cBar <- transSummary %>% dplyr::filter(., week == "8" & trt == trt.i) %>% dplyr::select(., propInfected, sepropI)
  beta_cVec <- simulateData(mean = beta_cBar$propInfected, SE = beta_cBar$sepropI) 
  
  beta_iBar <- transSummary %>% dplyr::filter(., week == "12" & trt == trt.i) %>% dplyr::select(., propInfected, sepropI)
  beta_iVec <- simulateData(mean = beta_iBar$propInfected, SE = beta_iBar$sepropI)
  
  alpha_cVec <- beta_cVec/0.6
  alpha_iVec <- beta_iVec/0.6
  
  #### Latent and incubation period
  # Latent period: take as uniform distribution between 4 and 21 days (or as a rate, 1/4 - 1/21),
  # 4 is based on Hill and Purcell 1997
  # 21 is based on significant transmission at 3 weeks
  deltaVec <- runif(n = 1000, min = 1/21, max = 1/4)
  
  # Incubation period: For susceptible line, assume incubation period is between 9 and 12 weeks,
  # Because I have no evidence of change in preference or in transmission at 8 weeks but I do have changes at 12 weeks
  gammaVec <- runif(n = 1000, min = 1/(12*7), max = 1/(9*7))
  
  #### Host and vector recovery
  # Take vector recovery from Fabien's experiment with WT plants
  lambdaVec <- simulateData(mean = 0.083, SE = 0.042)
  
  # Take host recovery as the change in proportion of source plants that were negative between 3 and 12 weeks
  # sourcepInfected for 12S - sourcepInfected for 3S
  prop.i <- as.numeric(transSummary[transSummary$week == 3 & transSummary$trt == trt.i, "sourcepInfected"])
  prop.f <- as.numeric(transSummary[transSummary$week == 12 & transSummary$trt == trt.i, "sourcepInfected"])
  nuVec <- ifelse(prop.i >= prop.f, prop.i - prop.f, 0)
  
  #### Feeding time 
  # Estimate of T is taken from Almeida and Backus 2004, Table 4, WDEI calculation for waveform C1 (ingestion)
  # Duration originaly given in minutes, convert to fraction of a day
  TVec <- simulateData(mean = 48.1/60, SE = 12.2/60)
  
  #### Combine all parameter values
  # Parameter vector order for simulation function
  # alpha_c = x[1], alpha_i = x[2], 
  # beta_c = x[3], beta_i = x[4],
  # phi_pc = x[5], phi_pi = x[6], phi_muc = x[7], phi_mui = x[8], 
  # delta = x[9], gamma = x[10], 
  # nu = x[11], lambda = x[12], 
  # T = x[13])
  par.i <- cbind(alpha_cVec, alpha_iVec,
                 beta_cVec, beta_iVec,
                 phi_pcVec, phi_piVec, phi_mucVec, phi_muiVec,
                 deltaVec, gammaVec,
                 nuVec, lambdaVec,
                 TVec)
  parList[[i]] <- par.i
}

###########################################################################################
#### Numerical Simulations
###########################################################################################

# Initial state variables
S0 <- 99; I0 <- 1; E0 <- 0; C0 <- 0
U0 <- 200; V0 <- 0

ti <- proc.time()
simList <- lapply(parList, function(x) apply(x, 1, SECIMSimulations) %>% rbindlist() %>% as.data.frame())
tf <- proc.time()
# Time the simulations took:
(tf-ti)/60


simList <- lapply(1:length(simList), function(x) simList[[x]][,-1])

summaryList <- lapply(simList, function(x) data.frame("mean" = apply(x, 2, median),
                                                      "sd" = apply(x, 2, sd),
                                                      "cil" = apply(x, 2, function(x) quantile(x, 0.025)),
                                                      "ciu" = apply(x, 2, function(x) quantile(x, 0.975)),
                                                      "max" = apply(x, 2, max)))
summaryList

maxi <- which(simS$V_c == max(simS$V_c))
maxRes <- simS[c(maxi-1, maxi, maxi+1),]
maxpar <- parS[c(maxi-1, maxi, maxi+1),]


meanParS <- colMeans(parS)

dynS <- SECIMDynamics(meanParS)
matplot(dynS[1:50,-1], type = "l",
        col = c("black", "pink", "green", "brown", "grey", "green", "brown"), 
        lty = c(1,1,1,1,3,3,3), lwd = 2)
legend("topright", c("S", "E", "C", "I", "U", "V_c", "V_i"), 
       col = c("black", "pink", "green", "brown", "grey", "green", "brown"), 
       lty = c(1,1,1,1,3,3,3), lwd = 2)
