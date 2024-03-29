#### Assumptions tests for consumers movement model for PdR1 preference-transmission experiment

rm(list = ls())
# libraries
my.packages <- c("openxlsx", "tidyr", "dplyr", "ggplot2", "data.table",
                 "lattice", "optimx", "bbmle", "numDeriv", "stringr",
                 "survival")
lapply(my.packages, require, character.only = TRUE)

## Input functions for P1 and P2 equations, 
## negative log likelihood (NLL) functions, 
## gradient functions, and AIC function
source("R_functions/consumer_movement_model_source_functions.R")
# factor2numeric() function
source("R_functions/factor2numeric.R")
# NLL gradient functions
source("R_functions/cmm_gradient_functions.R")


########################################################################################################################################
#### Importing data 
atdata <- read.xlsx("data/2017_data/PdR1_2017_preference-transmission_experiment_data.xlsx", sheet = "assumptions_test", detectDates = TRUE)
## Remove extraneous columns
## I can remove week, block, genotype, trt because they are all the same
## rep column distinguishes different trials, need to keep
atdata <- atdata %>% dplyr::select(rep, time_from_start_hr, xf_plant, test_plant, neutral_space, dead, missing)
str(atdata)
print.data.frame(atdata)


#########################################################################################################################################
#### Survival analysis for constant rate assumption
#########################################################################################################################################

#### Restructure data for survival analysis
## Create an empty data.frame for the survival data
nreps <- length(unique(atdata$rep)) # number of trials/replicates
survData <- data.frame(rep = rep(NA, nreps),
                       time0 = rep(NA, nreps),
                       moved0 = rep(NA, nreps),
                       choice1 = rep(NA, nreps),
                       moved_from1 = rep(NA, nreps),
                       time_from1 = rep(NA, nreps))
choices <- c("xf_plant", "test_plant", "neutral_space", "dead", "missing")

## For loop to convert to survival analysis data structure
for(i in 1:nreps){
  rep.i <- atdata[atdata$rep == i,]
  survData$rep[i] <- i
  survData$time0[i] <- rep.i %>% dplyr::filter(is.na(neutral_space)) %>% dplyr::select(time_from_start_hr) %>% min()
  choices.i <- rep.i %>% dplyr::filter(time_from_start_hr == survData$time0[i]) %>% dplyr::select_(.dots = choices) 
  survData$choice1[i] <- names(choices.i)[which(choices.i == 1)]
  survData$moved0[i] <- ifelse(is.na(survData$choice1[i]), 0, 1)
  rep.i2 <- rep.i %>% dplyr::filter(rep.i$time_from_start_hr >= survData$time0[i])
  choice.i2 <- rep.i2 %>% dplyr::select(all_of(survData$choice1[i]))
  survData$time_from1[i] <- ifelse(any(is.na(choice.i2)), min(rep.i2[which(is.na(choice.i2)), "time_from_start_hr"]), NA) 
  survData$moved_from1[i] <- ifelse(is.na(survData$time_from1[i]), 0, 1)
}

survData

## Fix the "missing" observations "by hand"
## Delete rep 5 because the bug died quickly
survData <- survData %>% dplyr::filter(rep != 5)
survData[survData$rep == 13,] <- c(13, 4, 1, "test_plant", 0, NA)
survData[survData$rep == 15,] <- c(15, 6, 1, "test_plant", 0, NA)
survData
## Convert time and moved columns to numeric
for(i in c(2,3,5,6)){
  survData[,i] <- as.numeric(survData[,i])
}
str(survData)

#### Settling time
# Obtaining Kaplan-Meier estimates of (S(t)) and time by plant choice for attraction rates
# choice1 is a variable that records the first choice of each individual consumer
# time0 is the time that each consumer moved from the neutral space to its first choice
# moved0 is a binary variable recording whether each consumer moved (1) or did not move (0)
setl.fit = survfit(Surv(time0, moved0) ~ choice1, data = survData) # Calculate K-M survival
sumsetl <- summary(setl.fit)
# Extract time, S(t), and choice from the survfit object
setl.dat <- data.frame("time" = sumsetl$time,
                       "surv" = sumsetl$surv,
                       "choice" = as.character(sumsetl$strata)) %>%
  mutate(time = as.numeric(time),
         surv = as.numeric(surv),
         rate = "Attraction")

## Leaving time
# Obtaining Kaplan-Meier estimates of (S(t)) and time by plant choice for leaving rates
# time1 is the time that each consumer its first choice
# moved1 records whether each consumer left the first choice (1) or did not leave (0)
leav.fit = survfit(Surv(time_from1, moved_from1) ~ choice1, data = survData) # Calculate K-M survival
sumleav <- summary(leav.fit)
# Extract time, S(t), and choice from the survfit object
leav.dat <- data.frame("time" = sumleav$time,
                       "surv" = sumleav$surv,
                       "choice" = as.character(sumleav$strata)) %>%
  mutate(time = as.numeric(time),
         surv = as.numeric(surv),
         rate = "Leaving")


#### Combine settling data and leaving data
survdat2 <- rbind(setl.dat, leav.dat)
survdat2 <- survdat2[survdat2$surv > 0,] # Remove K-M estimates of 0
survdat2$log.time <- log(survdat2$time) 
survdat2$log.surv <- log(-log(survdat2$surv))
survdat2

#### Plot K-M survival estimates vs. time
cmm_constant_rate_plot <- xyplot(log.surv ~ log.time|choice*rate, data = survdat2, 
                                 scales = list(alternating = FALSE, tck = c(1, 0), cex = 1.1),
                                 xlab = list("ln(time)", cex = 1.3),
                                 ylab = list("ln{-ln[S(t)]}", cex = 1.3), aspect = 1,
                                 layout = c(2,2), as.table = TRUE, strip = TRUE,
                                 pch = 16, col = "black", type = c('p', 'r'))

trellis.device(device = "tiff", file = "results/figures/2017_figures/cmm_assumptions_kaplan_meier_plots.tif")
  print(cmm_constant_rate_plot)
dev.off()


###############################################################################
#### Inspection of correlation matrices for consecutive choices assumption
###############################################################################


#### Contingency table
cchoiceData <- read.xlsx("data/2017_data/PdR1_2017_preference-transmission_experiment_data.xlsx", sheet = "choice_contingency_table", detectDates = TRUE)
str(cchoiceData)
(ctable <- with(cchoiceData, table(choice1, choice2)))
## Too few movement events for contingency table analysis


###############################################################################
#### Examine correlation matrices from CMM instead for independent choice assumption
#### 2016 MLE correlation matrices
corrMat16 <- readRDS("output/cmm_parameter_correlation_matrices_2016.rds")

#### Make correlation matrices look pretty for supplementary material
corrTableList <- lapply(corrMat16, extractCorrelationMatrix)

## Write all correlation matrices to a single Excel worksheet
wb <- createWorkbook()
addWorksheet(wb, sheetName = "corr_matrices_2016")
startRows <- seq(2,2+6*6,by=6)

for(i in 1:length(corrTableList)){
  writeData(wb, sheet = "corr_matrices_2016", x = corrTableList[[i]], startCol = 2, startRow = startRows[i])
  writeData(wb, sheet = "corr_matrices_2016", x = names(corrTableList)[i], startCol = 1, startRow = startRows[i])
}


#### How many correlations across week.genotype combinations are < 0.5?
## Note: low correlations are double counted but so are the totals, so they cancel each other out in the proportion
lowcor <- totals <- rep(0, length(corrTableList))
for(i in 1:length(corrTableList)){
  cor.i <- corrTableList[[i]][,-1]
  lowcor[i] <- sum(cor.i < 0.5)
  totals[i] <- sum(dim(cor.i)) - sum(cor.i == 1)
}
## Total proportion of parameter correlations that are < 0.5
sum(lowcor)/sum(totals)


###############################################################################
#### 2017 MLE correlation matrices

corrMat17 <- readRDS("output/cmm_parameter_correlation_matrices_2017.rds")

corrTableList17 <- lapply(corrMat17, extractCorrelationMatrix)

## Write all correlation matrices to a single Excel worksheet
addWorksheet(wb, sheetName = "corr_matrices_2017")
startRows <- seq(2,2+6*length(corrTableList17),by=6)

for(i in 1:length(corrTableList17)){
  writeData(wb, sheet = "corr_matrices_2017", x = corrTableList17[[i]], startCol = 2, startRow = startRows[i])
  writeData(wb, sheet = "corr_matrices_2017", x = names(corrTableList17)[i], startCol = 1, startRow = startRows[i])
}

saveWorkbook(wb, file = "results/cmm_assumptions_correlation_tables.xlsx", overwrite = TRUE)


#### How many correlations across week.genotype combinations are < 0.5?
lowcor <- totals <- rep(0, length(corrTableList17))
for(i in 1:length(corrTableList17)){
  cor.i <- corrTableList17[[i]][,-1]
  lowcor[i] <- sum(cor.i < 0.5, na.rm = TRUE)
  totals[i] <- sum(dim(cor.i)) - sum(is.na(cor.i)) - sum(cor.i == 1, na.rm = TRUE)
}
## Proportion of correlations that were below 0.5
sum(lowcor)/sum(totals)
