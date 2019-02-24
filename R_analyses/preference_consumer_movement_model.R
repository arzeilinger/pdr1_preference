#### ANALYSIS OF BGSS PREFERENCE ON PDR1 RESISTANT AND SUSCEPTIBLE INFECTED PLANTS
#### 2016 and 2017 data
#### USING CONSUMER MOVEMENT MODEL
# Code includes maximum likelihood estimation of 4 model variants, model selection,
# estimation of variance using normal approximation method
# model averaging, and calculation of equilibrial probabilities, Pj*

rm(list = ls())
# libraries
# loading dtplyr that replaces dplyr and data.table
my.packages <- c("openxlsx", "tidyr", "dplyr", "ggplot2", "data.table",
                 "lattice", "optimx", "bbmle", "numDeriv", "stringr",
                 "googlesheets", "RColorBrewer")
lapply(my.packages, require, character.only = TRUE)

## Input functions for P1 and P2 equations, 
## negative log likelihood (NLL) functions, 
## gradient functions, and AIC function
source("R_functions/consumer_movement_model_source_functions.R")
# Plotting functions using lattice
source("R_functions/lattice_plotting_functions.R")
# factor2numeric() function
source("R_functions/factor2numeric.R")
# NLL gradient functions
source("R_functions/cmm_gradient_functions.R")



###################################################################################################
###################################################################################################
#### 2016 experimental data
###################################################################################################
#### Importing data 

prefdata <- read.xlsx("data/2016_data/pdr1_preference_data.xlsx", sheet = "data")
str(prefdata)

# Separate leaf and symptom data from preference count data
leafdata <- prefdata[,13:ncol(prefdata)]
prefdata <- prefdata[,1:12]
prefdata$cage <- paste(prefdata$trt, prefdata$rep, sep = "")
head(prefdata)

# Checking out weirdness in S3 trial data
s3prefdata <- dplyr::filter(prefdata, week == 3, trt == "S")
dplyr::filter(s3prefdata, time_from_start_hr == 1 | time_from_start_hr == 2)


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
cmmData$N <- cmmData %>% with(., n1 + n2 + n3)
cmmData$week.trt <- factor(paste(cmmData$week, cmmData$trt, sep = ""))

r12Data <- dplyr::filter(cmmData, week == 12, trt == "R") 
dplyr::filter(cmmData, week == 3)


#### Fit models to each of the trt-week combinations
modelFits <- lapply(levels(cmmData$week.trt), function(x) optimizeCMM(dat = cmmData[cmmData$week.trt == x,], upperConstraint = 64, aiccN = 8))
names(modelFits) <- levels(cmmData$week.trt)
saveRDS(modelFits, file = "output/CMM_optimx_model_selection_output.rds")

modelFits <- readRDS("output/CMM_optimx_model_selection_output.rds")

#### Calculate variance-covariance and correlation matrices for each trt-week combination
matrices <- lapply(modelFits, getParCorrelations)
## Save vcov and corr matrices list
saveRDS(matrices, file = "output/cmm_parameter_correlation_matrices_2016.rds")

#### Extract and organize parameter estimates from all models
paramResults <- lapply(modelFits, mleTable)
names(paramResults) <- levels(cmmData$week.trt)

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

dplyr::filter(paramData, week.trt == "3S")

dplyr::filter(paramData, week.trt == "12S" | week.trt == "12R")
testpardat <- dplyr::filter(paramData, week.trt == "12R")
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

write.csv(plotPars, file = "results/pdr1_cmm_rate_parameter_estimates.csv", row.names = FALSE)

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
#tiff(filename = "results/figures/pdr1_cmm_rate_parameter_plot.tif")
     # width = 76*2, height = 76, units = "mm", 
     # res = 600, compression = "lzw")
#  plot.new()
parameter_plot <-  with(plotPars,
                        xyplot(estimate ~ dummyx|rate, groups = latticegroups,
                               ly = cil, uy = ciu,
                               scales = list(col = 1, alternating = 1, tck = c(1, 0), cex = 1.1, relation = "free",
                                             x = list(limits = c(0, 13), at = c(3,8,12),
                                                      labels = list(c("","", ""), c(3,8,12))),
                                             y = list(limits = list(c(0, 3), c(0, 0.6)),
                                                      at = list(seq(0, 3, by = 1), seq(0, 0.6, by = 0.2)),
                                                      labels = list(seq(0, 3, by = 1), seq(0, 0.6, by = 0.2)))),
                               xlab = list("Weeks post-inoculation", cex = 1.2), 
                               ylab = list(expression("Leaving rate "(hr^{-1})*"                          Attraction rate "(hr^{-1})),
                                           cex = 1.2),
                               layout = c(1,2), as.table = TRUE, strip = FALSE, pch = c(19, 1, 17, 2),
                               type = 'p', cex = 1.2, col = "black", 
                               key = list(x = 0.4, y = 0.85, corner = c(0,0),
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

trellis.device(device = "tiff", file = "results/figures/pdr1_cmm_rate_parameter_plot_test.tif")
  print(parameter_plot)
dev.off()






#####################################################################################
#### Equilibrial probabilities 
#####################################################################################

sd <- sqrt(test$variance)

p1sim <- simulateData(estimates[1], sd[1])
p2sim <- simulateData(estimates[2], sd[2])
mu1sim <- simulateData(estimates[3], sd[3])
mu2sim <- simulateData(estimates[4], sd[4])
paramSim <- cbind(p1sim, p2sim, mu1sim, mu2sim)

P1eq.func(as.numeric(paramSim[50,]))
P1eqsim <- apply(paramSim, 1, P1eq.func)

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

# Add columns for weeks and treatments, and remove 12.2 week trials
plotProbs <- dplyr::filter(ProbResults, week.trt != "12.2R" & week.trt != "12.2S")
plotProbs$week <- plotProbs$week.trt %>% str_extract(., "[0-9]+") %>% as.numeric()
plotProbs$trt <- plotProbs$week.trt %>% str_extract(., "[aA-zZ]+")
# Split parameter into rate and choice columns, and replace values with more meaningful terms
plotProbs$choice <- plotProbs$state %>% str_extract(., "[0-9]+") %>% as.numeric() %>%
  gsub(1, "infected", ., fixed = TRUE) %>%
  gsub(2, "Xf-free", ., fixed = TRUE)
plotProbs$latticegroups <- with(plotProbs, paste(trt, choice, sep=" ")) %>% factor()
plotProbs

write.csv(plotProbs, file = "results/pdr1_cmm_equilibrial_probabilities.csv", row.names = FALSE)

# Create dummy x variable to space out points
adj <- c(-0.5, -0.25, 0.25, 0.5)
plotProbs$dummyx <- 0
for(i in 1:length(levels(plotProbs$latticegroups))){
  group.i <- levels(plotProbs$latticegroups)[i]
  data.i <- which(plotProbs$latticegroups == group.i)
  adj.i <- adj[i]
  plotProbs$dummyx[data.i] <- plotProbs[plotProbs$latticegroups == group.i,"week"] + adj.i
}


# Rate parameter plots together
#tiff(filename = "results/figures/pdr1_cmm_rate_parameter_plot.tif")
# width = 76*2, height = 76, units = "mm", 
# res = 600, compression = "lzw")
#  plot.new()
probabilities_plot <-  with(plotProbs,
                          xyplot(median ~ dummyx, groups = latticegroups,
                                 ly = cil, uy = ciu,
                                 scales = list(col = 1, alternating = 1, tck = c(1, 0), cex = 1.1, relation = "free"),
                                               # x = list(limits = c(0, 13), at = c(3,8,12),
                                               #          labels = list(c("","", ""), c(3,8,12))),
                                               # y = list(limits = list(c(0, 3), c(0, 0.6)),
                                               #          at = list(seq(0, 3, by = 1), seq(0, 0.6, by = 0.2)),
                                               #          labels = list(seq(0, 3, by = 1), seq(0, 0.6, by = 0.2)))),
                                 xlab = list("Weeks post-inoculation", cex = 1.2), 
                                 ylab = list("Equilibrial probability", cex = 1.2),
                                 layout = c(1,1), as.table = TRUE, strip = FALSE, pch = c(19, 1, 17, 2),
                                 type = 'p', cex = 1.2, col = "black", 
                                 key = list(x = 0.4, y = 0.85, corner = c(0,0),
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

trellis.device(device = "tiff", file = "results/figures/pdr1_cmm_equilibrial_probabilities_plot.tif")
  print(probabilities_plot)
dev.off()



##########################################################################################################
#### Estimating rate parameters for each preference cage
#### To be used in logistic regression model for transmission
#### Only using the best estimates for each cage; not sure how to incorporate variance


# Checking out weirdness in S3 trial data
s3prefdata <- dplyr::filter(prefdata, week == 3, trt == "S")
dplyr::filter(s3prefdata, time_from_start_hr == 1 | time_from_start_hr == 2)
s3prefdata

# Add identifier column for each cage
cmmDataCage <- data.frame("t" = prefdata$time_from_start_hr,
                          "n1" = prefdata$source_plant,
                          "n2" = prefdata$test_plant,
                          "n3" = prefdata$neutral_space,
                          "N" = with(prefdata, source_plant + test_plant + neutral_space),
                          "week.cage" = factor(paste(prefdata$week, prefdata$cage, sep = "")))
length(levels(cmmDataCage$week.cage))
head(cmmDataCage)


#### Fit models to each of the trt-week combinations
modelFitsCage <- lapply(levels(cmmDataCage$week.cage), function(x) optimizeCMM(dat = cmmDataCage[cmmDataCage$week.cage == x,], upperConstraint = 8, aiccN = 8))
names(modelFitsCage) <- levels(cmmDataCage$week.cage)
saveRDS(modelFitsCage, file = "output/CMM_optimx_model_selection_output_per_cage.rds")

#modelFits <- readRDS("output/CMM_optimx_model_selection_output_per_cage.rds")

#### Extract and organize parameter estimates from all models
paramResultsCage <- lapply(modelFitsCage, function(x) tryCatch(mleTable(x), error = function(e) NA))
paramResultsCage <- lapply(modelFitsCage, mleTable)

#### Average results for all good models
paramAverageCage <- lapply(paramResultsCage, function(x) tryCatch(averageModels(x, dAIC.threshold = 7), error = function(e) NA))

paramAverageCage2 <- lapply(paramAverageCage, 
                            function(x) if(is.na(x)) data.frame("parameter" = c("p1", "p2", "mu1", "mu2"), 
                                                                "estimate" = rep(NA, 4),
                                                                "variance" = rep(NA, 4)) else x)

# Add week.trt combination to each element of the list and combine into one data.frame
for(i in 1:length(paramAverageCage)){
  paramAverageCage[[i]]$week.cage <- levels(cmmDataCage$week.cage)[i]
} 
paramDataCage <- rbindlist(paramAverageCage) %>% as.data.frame()
paramDataCage$estimate <- factor2numeric(paramDataCage$estimate)
paramDataCage$variance <- factor2numeric(paramDataCage$variance)

saveRDS(paramDataCage, file = "output/CMM_rate_parameters_per_cage.rds")



###################################################################################################
###################################################################################################
#### 2017 experimental data
###################################################################################################

#### Importing data 
# Import from local .xlsx file
# prefdata <- read.xlsx("data/2017_data/PdR1_2017_preference-transmission_experiment_data.xlsx", sheet = "BGSS_count_data", detectDates = TRUE)
# Import from Googlesheets
pdr1DataURL <- gs_url("https://docs.google.com/spreadsheets/d/14uJLfRL6mPrdf4qABeGeip5ZkryXmMKkan3mJHeK13k/edit?usp=sharing",
                      visibility = "private")
prefdataGS <- gs_read(pdr1DataURL, ws = "BGSS_count_data")
prefdata <- prefdataGS
str(prefdata)

# Separate leaf and symptom data from preference count data
leafdata <- prefdata[,15:ncol(prefdata)]
prefdata <- prefdata[,1:14]
prefdata$genotype <- factor(prefdata$genotype)
head(prefdata)


## Max number of BGSS for each trial is 8
# How many observations have >8 total bugs?
prefdata <- prefdata %>% dplyr::mutate(total_bgss = xf_plant + test_plant + neutral_space + dead + missing)
prefdata %>% dplyr::filter(total_bgss > 8)
# 2-2-007S-2 trials had 9 BGSS in them, by accident. Otherwise, data look good.


## Check howm many bugs I have in the freezer
# This should be roughly the same as the number of bugs at the 4d time point
freezerBugs <- prefdata %>% dplyr::filter(time_from_start_hr == 96) %>% dplyr::mutate(total_live_bugs = xf_plant + test_plant + neutral_space)
(totalFreezerBugs <- sum(freezerBugs$total_live_bugs)) 
# Total freezer bugs = 917
# Number of extraction plates
totalFreezerBugs/96


##############################################################################################################
#### Estimating attraction and leaving rates from Consumer Movement Model
#### Calculate total counts among cages for each week, genotype, and time point
# n1 = source plant (Xylella infected)
# n2 = test plant
# n3 = neutral space

cmmData <- prefdata %>% group_by(week, genotype, time_from_start_hr) %>% 
  summarise(n1 = sum(xf_plant, na.rm = TRUE),
            n2 = sum(test_plant, na.rm = TRUE),
            n3 = sum(neutral_space, na.rm = TRUE)) %>%
  as.data.frame()
cmmData$t <- cmmData$time_from_start_hr
cmmData$N <- cmmData %>% with(., n1 + n2 + n3)
cmmData$week.genotype <- factor(paste(cmmData$week, cmmData$genotype, sep = "-"))

# Check an example
data_5_102 <- dplyr::filter(cmmData, week == 5, genotype == "102") 
data_14_102 <- dplyr::filter(cmmData, week == 14, genotype == "102") 
data_14_094 <- dplyr::filter(cmmData, week == 14, genotype == "094") 

# Max N for any row of cmmData is 8*8 = 64
# Any rows of cmmData have more than N = 64 BGSS?
cmmData %>% dplyr::filter(N > 64)
cmmData %>% filter(week == 2 & genotype == "007")


#### Fit models to each of the week-genotype combinations
modelFits <- lapply(levels(cmmData$week.genotype), function(x) optimizeCMM(dat = cmmData[cmmData$week.genotype == x,], upperConstraint = 64, aiccN = 8))
names(modelFits) <- levels(cmmData$week.genotype)
saveRDS(modelFits, file = "output/CMM_optimx_model_selection_output_2017.rds")

modelFits <- readRDS("output/CMM_optimx_model_selection_output_2017.rds")

#### Calculate variance-covariance and correlation matrices for each trt-week combination
matrices <- lapply(modelFits, getParCorrelations)

#### Extract and organize parameter estimates from all models
paramResults <- lapply(modelFits, mleTable)
names(paramResults) <- levels(cmmData$week.genotype)


## Some variances are negative, need to figure out why
exop <- modelFits[[14]]$op.list$fixed
exhess <- hessian(func = NLLlist$fixed, x = as.numeric(exop[,grep("p",names(exop))]))
diag(solve(exhess))

exop2 <- modelFits[[14]]$op.list$mu.choice
exhess2 <- hessian(func = NLLlist$mu.choice, x = as.numeric(exop2[,grep("p",names(exop2))]))
diag(solve(exhess2))


#### Average results for all good models
paramAverage <- lapply(paramResults, function(x) averageModels(x, dAIC.threshold = 7))

# Add week.genotype combination to each element of the list and combine into one data.frame
for(i in 1:length(paramAverage)){
  paramAverage[[i]]$week.genotype <- levels(cmmData$week.genotype)[i]
} 
paramData <- rbindlist(paramAverage) %>% as.data.frame()
paramData$estimate <- factor2numeric(paramData$estimate)
paramData$variance <- factor2numeric(paramData$variance)

#### Calculate SE and 95% CI
# Some variances are negative, not sure why, but need to fix.
# Comment on ResearchGate forum (https://www.researchgate.net/post/In_R_how_to_estimate_confidence_intervals_from_the_Hessian_matrix) suggested taking absolute value of of variances
# Do this for now, but probably need to bootstrap/jackknife
# Also need to re-think model averaging of parameter estimates overall
paramData$se <- sqrt(abs(paramData$variance))
paramData$cil <- with(paramData, estimate - 1.96*se)
paramData$ciu <- with(paramData, estimate + 1.96*se)


#####################################################################################
#### Plotting rate parameters
# Structure data.frame for plotting
plotPars <- paramData
# Split week.trt into week and trt columns
week.genotype.split <- tstrsplit(plotPars$week.genotype, split = "-")
plotPars$week <- week.genotype.split[[1]] %>% as.numeric()
plotPars$genotype <- week.genotype.split[[2]] %>% factor()
# Split parameter into rate and choice columns, and replace values with more meaningful terms
plotPars$rate <- plotPars$parameter %>% str_extract(., "[aA-zZ]+") %>% 
  gsub("p", "attraction", ., fixed = TRUE) %>% 
  gsub("mu", "leaving", ., fixed = TRUE)
plotPars$choice <- plotPars$parameter %>% str_extract(., "[0-9]+") %>% as.numeric() %>%
  gsub(1, "infected", ., fixed = TRUE) %>%
  gsub(2, "Xf-free", ., fixed = TRUE) %>% factor()
plotPars$latticegroups <- with(plotPars, paste(genotype, choice, sep=" ")) %>% factor()
plotPars <- plotPars %>% arrange(., week)
plotPars

write.csv(plotPars, file = "results/pdr1_cmm_rate_parameter_estimates_2017.csv", row.names = FALSE)

# Create dummy x variable to space out points
# adj <- c(-0.5, -0.25, 0.25, 0.5)
# plotPars$dummyx <- 0
# for(i in 1:length(levels(plotPars$latticegroups))){
#   group.i <- levels(plotPars$latticegroups)[i]
#   data.i <- which(plotPars$latticegroups == group.i)
#   adj.i <- adj[i]
#   plotPars$dummyx[data.i] <- plotPars[plotPars$latticegroups == group.i,"week"] + adj.i
# }


#### Lattice plot -- currently doesn't work
# # Rate parameter plots together
# #tiff(filename = "results/figures/pdr1_cmm_rate_parameter_plot.tif")
# # width = 76*2, height = 76, units = "mm", 
# # res = 600, compression = "lzw")
# #  plot.new()
# parameter_plot <-  with(plotPars,
#                         xyplot(estimate ~ dummyx|rate*genotype, groups = choice,
#                                ly = cil, uy = ciu,
#                                scales = list(col = 1, alternating = 1, tck = c(1, 0), cex = 1.1, relation = "free",
#                                              x = list(limits = c(0, 15), at = c(2,5,8,14),
#                                                       labels = list(c("","","",""), c(2,5,8,14),
#                                                                     c("","","",""), c(2,5,8,14),
#                                                                     c("","","",""), c(2,5,8,14),
#                                                                     c("","","",""), c(2,5,8,14)))),
#                                              # y = list(limits = list(c(0, 3), c(0, 0.6)),
#                                              #          at = list(seq(0, 3, by = 1), seq(0, 0.6, by = 0.2)),
#                                              #          labels = list(seq(0, 3, by = 1), seq(0, 0.6, by = 0.2)))),
#                                xlab = list("Weeks post-inoculation", cex = 1.2), 
#                                ylab = list(expression("Leaving rate "(hr^{-1})*"                          Attraction rate "(hr^{-1})),
#                                            cex = 1.2),
#                                layout = c(4,2), as.table = TRUE, strip = FALSE, pch = c(19, 1),
#                                type = 'p', cex = 1.2, col = "black", 
#                                key = list(x = 0.4, y = 0.85, corner = c(0,0),
#                                           text = list(lab = levels(choice)), 
#                                           points = list(pch = c(19, 1), col = "black")), 
#                                prepanel = prepanel.ci,                      
#                                panel = function(x, y, ...) {                
#                                  panel.abline(v = unique(as.numeric(x)),  
#                                               col = "white")              
#                                  panel.superpose(x, y, ...)               
#                                },                                          
#                                panel.groups = panel.ci)) 
# parameter_plot
# # ltext(165, 9, "2009", cex = 1.3)
# # ltext(380, 9, "2010", cex = 1.3)
# # mtext(c("A", "B"), side = 3, cex = 1.3, adj = rep(0.06, 2), padj = c(-1.4, 15))
# # mtext(c("C", "D"), side = 3, cex = 1.3, adj = rep(0.64, 2), padj = c(-1.4, 15))
# 
# trellis.device(device = "tiff", file = "results/figures/2017_figures/pdr1_cmm_rate_parameter_plot_2017.tif")
# print(parameter_plot)
# dev.off()




#### Plotting using ggplot2
plotPars <- read.csv("results/pdr1_cmm_rate_parameter_estimates_2017.csv", header = TRUE)

# Prepend zeroes to genotype labels and make a factor
plotPars$genotype <- plotPars$genotype %>% formatC(., width = 3, format = "d", flag = "0") %>% factor()

# Color palette
my_palette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
my_palette = c(brewer.pal(5, "Blues")[c(4)], brewer.pal(5, "Set1")[c(3)])

# Plot
parameter_plot2 <- ggplot(data=plotPars, aes(x=week, y=estimate, group = choice, color = choice)) +
  #geom_line(aes(linetype=inoc.time, colour = trt), size=1.25) +
  geom_point(data = plotPars, size=3.5, position = position_dodge(width = 0.9)) +
  geom_errorbar(data = plotPars, aes(ymax=ciu, ymin=cil), width=0.2, position = position_dodge(width = 0.9)) +
  facet_grid(rate~genotype, scales = "free_y") +
  geom_hline(linetype = 2, yintercept = 0) +
  scale_color_manual(values = my_palette) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = unique(plotPars$week)) + 
  scale_y_continuous(name = "Rate (per hour)") +
                     #limits = c(0,10)) +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

parameter_plot2

ggsave("results/figures/2017_figures/pdr1_cmm_rate_parameter_plot_2017.jpg", plot = parameter_plot2,
       width = 14, height = 7, units = "in")




#### Plot raw numbers of bugs over time
plotCounts <- cmmData %>% mutate(xf_plant = n1/N, test_plant = n2/N) %>%
  dplyr::select(-n1, -n2, -n3, -t, -N, -week.genotype) %>% 
  gather(., key = "choice", value = "proportion", xf_plant, test_plant)

## Change choice names to "infected" and "xf-free"
plotCounts$choice <- ifelse(plotCounts$choice == "xf_plant", "infected", "xf_free")

bgss_counts_plot <- ggplot(data=plotCounts, aes(x=time_from_start_hr, y=proportion, group = choice)) +
  geom_line(aes(colour = choice), size=1.25) +
  #geom_point(aes(colour = choice), size=3.5, position = position_dodge(width = 0.9)) +
  #geom_errorbar(aes(ymax=ciu, ymin=cil), width=0.2, position = position_dodge(width = 0.9)) +
  facet_grid(week~genotype, scales = "fixed") +
  scale_color_manual(values = my_palette) +
  scale_x_continuous(name = "Time from start (hr)") + 
  scale_y_continuous(name = "Proportion of BGSS") +
  #limits = c(0,10)) +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

bgss_counts_plot

ggsave("results/figures/2017_figures/pdr1_bgss_counts_time-series_plot_2017.jpg", plot = bgss_counts_plot,
       width = 14, height = 7, units = "in")

#### Example of a single BGSS count plot
plotCounts1 <- plotCounts %>% dplyr::filter(week == 2 & genotype == "092")

bgss_counts_plot1 <- ggplot(data=plotCounts1, aes(x=time_from_start_hr, y=proportion, group = choice)) +
  geom_line(aes(colour = choice), size=1.25) +
  #geom_point(aes(colour = choice), size=3.5, position = position_dodge(width = 0.9)) +
  #geom_errorbar(aes(ymax=ciu, ymin=cil), width=0.2, position = position_dodge(width = 0.9)) +
  scale_color_manual(values = my_palette) +
  scale_x_continuous(name = "Time from start (hr)") + 
  scale_y_continuous(name = "Proportion of BGSS", limits = c(0,1)) +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank())
bgss_counts_plot1

ggsave("results/figures/2017_figures/pdr1_bgss_counts_time-series_single_plot_2017.jpg", plot = bgss_counts_plot1,
       width = 14, height = 7, units = "in")



##########################################################################################################
#### Estimating rate parameters for each preference cage
#### To be used in logistic regression model for transmission
#### Only using the best estimates for each cage; not sure how to incorporate variance

# Add identifier column for each cage
cmmDataCage <- with(prefdata, data.frame("t" = time_from_start_hr,
                                         "n1" = xf_plant,
                                         "n2" = test_plant,
                                         "n3" = neutral_space,
                                         "N" = xf_plant + test_plant + neutral_space,
                                         "cage" = factor(paste(week, block, genotype, trt, rep, sep = "-"))))
# Remove rows that are NA 
cmmDataCage <- cmmDataCage %>% dplyr::filter(!is.na(N))
length(levels(cmmDataCage$cage))
head(cmmDataCage)


#### Fit models to each of the trt-week combinations
modelFitsCage <- lapply(levels(cmmDataCage$cage), function(x) optimizeCMM(dat = cmmDataCage[cmmDataCage$cage == x,], upperConstraint = 8, aiccN = 8))
names(modelFitsCage) <- levels(cmmDataCage$cage)
saveRDS(modelFitsCage, file = "output/CMM_2017_optimx_model_selection_output_per_cage.rds")

#modelFits <- readRDS("output/CMM_optimx_model_selection_output_per_cage.rds")

#### Extract and organize parameter estimates from all models
paramResultsCage <- lapply(modelFitsCage, function(x) tryCatch(mleTable(x), error = function(e) NA))
paramResultsCage <- lapply(modelFitsCage, mleTable)

#### Average results for all good models
paramAverageCage <- lapply(paramResultsCage, function(x) tryCatch(averageModels(x, dAIC.threshold = 7), error = function(e) NA))

paramAverageCage2 <- lapply(paramAverageCage, 
                            function(x) if(is.na(x)) data.frame("parameter" = c("p1", "p2", "mu1", "mu2"), 
                                                                "estimate" = rep(NA, 4),
                                                                "variance" = rep(NA, 4)) else x)

# Add cage identifier to each element of the list and combine into one data.frame
for(i in 1:length(paramAverageCage)){
  paramAverageCage[[i]]$cage <- levels(cmmDataCage$cage)[i]
} 
paramDataCage <- rbindlist(paramAverageCage) %>% as.data.frame()
paramDataCage$estimate <- factor2numeric(paramDataCage$estimate)
paramDataCage$variance <- factor2numeric(paramDataCage$variance)

saveRDS(paramDataCage, file = "output/CMM_2017_rate_parameters_per_cage.rds")

#### Check some of the results
check1 <- paramDataCage %>% dplyr::filter(grepl("2-2-102R", cage))
check2 <- paramDataCage %>% dplyr::filter(grepl("14-2-007S", cage))
                                          