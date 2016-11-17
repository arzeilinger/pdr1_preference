#### ANALYSIS OF BGSS PREFERENCE ON PDR1 RESISTANT AND SUSCEPTIBLE INFECTED PLANTS
#### USING CONSUMER MOVEMENT MODEL
# Code includes maximum likelihood estimation of 4 model variants, model selection,
# estimation of variance using normal approximation method
# model averaging, and calculation of equilibrial probabilities, Pj*

rm(list = ls())
# libraries
# loading dtplyr that replaces dplyr and data.table
my.packages <- c("xlsx", "tidyr", "dtplyr", "ggplot2", 
                 "lattice", "optimx", "bbmle", "numDeriv", "stringr")
lapply(my.packages, require, character.only = TRUE)

## Input functions for P1 and P2 equations, 
## negative log likelihood (NLL) functions, 
## gradient functions, and AIC function
source("R_functions/consumer_movement_model_source_functions.R")
# Plotting functions using lattice
source("R_functions/lattice_plotting_functions.R")

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


#### Fit models to each of the trt-week combinations
modelFits <- lapply(levels(cmmData$week.trt), function(x) optimizeCMM(dat = cmmData[cmmData$week.trt == x,], upperConstraint = 64, aiccN = 8))
names(modelFits) <- levels(cmmData$week.trt)
saveRDS(modelFits, file = "output/CMM_optimx_model_selection_output.rds")

modelFits <- readRDS("output/CMM_optimx_model_selection_output.rds")

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

dplyr::filter(paramData, week.trt == "12S" | week.trt == "12R")
testpardat <- dplyr::filter(paramData, week.trt == "12S")
testpardat


#####################################################################################
#### Plotting rate parameters
# Structure data.frame for plotting
plotPars <- dplyr::filter(paramData, week.trt != "12.2R" & week.trt != "12.2S") # Remove second 12-week trials
# Split week.trt into week and trt columns
plotPars$week <- plotPars$week.trt %>% str_extract(., "[0-9]+") %>% as.numeric() 
plotPars$trt <- plotPars$week.trt %>% str_extract(., "[aA-zZ]+")
# Split parameter into rate and choice columns, and replace values with more meaningful terms
plotPars$rate <- plotPars$parameter %>% str_extract(., "[aA-zZ]+") %>% 
  gsub("p", "attraction", ., fixed = TRUE) %>% 
  gsub("mu", "leaving", ., fixed = TRUE)
plotPars$choice <- plotPars$parameter %>% str_extract(., "[0-9]+") %>% as.numeric() %>%
  gsub(1, "infected", ., fixed = TRUE) %>%
  gsub(2, "Xf-free", ., fixed = TRUE)
plotPars$latticegroups <- with(plotPars, paste(trt, choice, sep=" ")) %>% factor()
plotPars

write.csv("plotPars", file = "results/pdr1_cmm_rate_parameter_estimates.csv", row.names = FALSE)

# Create dummy x variable to space out points
adj <- c(-0.5, -0.25, 0.25, 0.5)
plotPars$dummyx <- 0
for(i in 1:length(levels(plotPars$latticegroups))){
  group.i <- levels(plotPars$latticegroups)[i]
  data.i <- which(plotPars$latticegroups == group.i)
  adj.i <- adj[i]
  plotPars$dummyx[data.i] <- plotPars[plotPars$latticegroups == group.i,"week"] + adj.i
}


# Rate parameter plots together
tiff(filename = "results/figures/pdr1_cmm_rate_parameter_plot.tif")
  plot.new()
  with(plotPars,
       xyplot(estimate ~ dummyx|rate, groups = latticegroups,
              ly = cil, uy = ciu,
              scales = list(col = 1, alternating = 1, tck = c(1, 0), cex = 1.1, relation = "free"),
                            # x = list(limits = c(0.85, 5), at = c(2, 4), 
                            #          labels = list(c("",""), c("",""), c("Damaged", "Undamaged"), c("Damaged", "Undamaged"))),
                            # y = list(limits = list(c(0, 1.2), c(0, 0.6), 
                            #                        c(0, 0.3), c(0, 0.15)),
                            #          at = list(seq(0, 1.2, by = 0.4), seq(0, 0.6, by = 0.2),
                            #                    seq(0, 0.3, by = 0.1), seq(0, 0.15, by = 0.05)),
                            #          labels = list(seq(0, 1.2, by = 0.4), seq(0, 0.6, by = 0.2),
                            #                        seq(0, 0.3, by = 0.1), seq(0, 0.15, by = 0.05)))),
              xlab = list("Plant choice", cex = 1.2), 
              #ylab = list(expression("Leaving rate "(individuals^{"  "}*hr^{-1})*"        Attraction rate "(individuals^{"  "}*hr^{-1})),
              ylab = list(expression("Leaving rate "(hr^{-1})*"                       Attraction rate "(hr^{-1})),
                          cex = 1.2),
              layout = c(1,2), as.table = TRUE, strip = TRUE, pch = c(19, 1, 17, 2),
              type = 'p', cex = 1.2, col = "black", 
              key = list(x = 0.55, y = 0.55, corner = c(0,0),
                         text = list(lab = levels(latticegroups)), 
                         points = list(pch = c(19, 1, 17, 2), col = "black")), 
              prepanel = prepanel.ci,                      
              panel = function(x, y, ...) {                
                panel.abline(v = unique(as.numeric(x)),  
                             col = "white")              
                panel.superpose(x, y, ...)               
              },                                          
              panel.groups = panel.ci))                    
  # ltext(165, 9, "2009", cex = 1.3)
  # ltext(380, 9, "2010", cex = 1.3)
  # mtext(c("A", "B"), side = 3, cex = 1.3, adj = rep(0.06, 2), padj = c(-1.4, 15))
  # mtext(c("C", "D"), side = 3, cex = 1.3, adj = rep(0.64, 2), padj = c(-1.4, 15))
dev.off()






#####################################################################################
#### Equilibrial probabilities for stink bug data
#####################################################################################

## Functions to calculate equilibrium probabilities, P1* and P2* from parameter estimates
# Equilibrium equations from Equation A12 in Appendix A of Zeilinger et al. (2014)
ProbResults <- lapply(unique(paramData$week.trt), function(x) ProbEqmc(paramData[paramData$week.trt == x,]))
  
# Add week.trt combination to each element of the list and combine into one data.frame
for(i in 1:length(ProbResults)){
  ProbResults[[i]]$week.trt <- levels(cmmData$week.trt)[i]
} 
ProbResults <- ProbResults %>% rbindlist() %>% as.data.frame()
# Restructure data set
names(ProbResults) <- c("state", "median", "cil", "ciu", "week.trt")
ProbResults$median <- factor2numeric(ProbResults$median)
ProbResults$cil <- factor2numeric(ProbResults$cil)
ProbResults$ciu <- factor2numeric(ProbResults$ciu)
ProbResults

# Add columns for weeks and treatments, and remove 12.2 week trials
ProbResults <- dplyr::filter(ProbResults, week.trt != "12.2R" & week.trt != "12.2S")
ProbResults$week <- ProbResults$week.trt %>% str_extract(., "[0-9]+") %>% as.numeric()
ProbResults$trt <- ProbResults$week.trt %>% str_extract(., "[aA-zZ]+")
ProbResults



