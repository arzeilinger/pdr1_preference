#### ANALYSIS OF BGSS PREFERENCE ON PDR1 RESISTANT AND SUSCEPTIBLE INFECTED PLANTS
#### 2016 and 2017 data
#### USING CONSUMER MOVEMENT MODEL
# Code includes maximum likelihood estimation of 4 model variants, model selection,
# estimation of variance using normal approximation method and model averaging

rm(list = ls())
# libraries
# loading dtplyr that replaces dplyr and data.table
my.packages <- c("openxlsx", "tidyr", "dplyr", "ggplot2", "data.table",
                 "lattice", "optimx", "BB", "bbmle", "numDeriv", "stringr",
                 "googlesheets", "RColorBrewer", "cowplot")
lapply(my.packages, require, character.only = TRUE)

## Input functions for P1 and P2 equations to consumer movement model (CMM)
## negative log likelihood (NLL) functions, 
## gradient functions, and AIC function
source("R_functions/consumer_movement_model_source_functions.R")
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

## Calculate number of dead BGSS
lastobs <- prefdata %>% dplyr::filter(time_from_start_hr == 192 & week != 12.2)
with(lastobs, table(week, trt))
totaldead <- sum(lastobs$dead, na.rm = TRUE)
totalbugs <- 48*8
## Proportion of bugs that died
totaldead/totalbugs
  

#### Calculate total counts among cages for each week, genotype, and time point
# Note: remove week == 12.2 as these trials were designed differently
# n1 = source plant (Xylella infected)
# n2 = test plant
# n3 = neutral space

cmmData <- prefdata %>% dplyr::filter(week != 12.2) %>%
  group_by(week, trt, time_from_start_hr) %>% 
  summarise(n1 = sum(source_plant),
            n2 = sum(test_plant),
            n3 = sum(neutral_space)) %>%
  as.data.frame()
cmmData$t <- cmmData$time_from_start_hr
cmmData$N <- cmmData %>% with(., n1 + n2 + n3)
cmmData$week.trt <- factor(paste(cmmData$week, cmmData$trt, sep = ""))

# Check subset of data
dplyr::filter(cmmData, week == 3)

#### Fit models to each of the trt-week combinations
## opimizeCMM() fits 4 model variants to each trt-week combination, calculates AIC, and returns optimx() output and AIC table
## Note: since revisiting my code after publication, optimizeCMM() no longer works properly...
## I suspect it has something to do with updated packages but I was unable to fix it with the time I had.
modelFits <- lapply(levels(cmmData$week.trt), function(x) optimizeCMM(dat = cmmData[cmmData$week.trt == x,], upperConstraint = 64, aiccN = 8))
names(modelFits) <- levels(cmmData$week.trt)

## Load previously saved MLE results
#modelFits <- readRDS("output/CMM_optimx_model_selection_output_2016.rds")

#### Extract model selection tables
selectTables16 <- lapply(1:length(modelFits), function(x) modelFits[[x]]$modelSelect)
names(selectTables16) <- names(modelFits)


#### Calculate variance-covariance and correlation matrices for each trt-week combination
matrices <- lapply(modelFits, getParCorrelations)

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
paramData$estimate <- as.numeric(paramData$estimate)
paramData$variance <- as.numeric(paramData$variance)

# Calculate SE and 95% CI
paramData$se <- sqrt(paramData$variance)
paramData$cil <- with(paramData, estimate - 1.96*se)
paramData$ciu <- with(paramData, estimate + 1.96*se)


#####################################################################################
#### Plotting rate parameters
#### Plotting rate parameters
# Structure data.frame for plotting
plotPars16 <- paramData
# Split parameter into rate and choice columns, and replace values with more meaningful terms
plotPars16$rate <- plotPars16$parameter %>% str_extract(., "[aA-zZ]+") %>% 
  gsub("p", "attraction", ., fixed = TRUE) %>% 
  gsub("mu", "leaving", ., fixed = TRUE)
plotPars16$choice <- plotPars16$parameter %>% str_extract(., "[0-9]+") %>% as.numeric() %>%
  gsub(1, "infected", ., fixed = TRUE) %>%
  gsub(2, "Xf-free", ., fixed = TRUE) %>% factor()
plotPars16$week <- plotPars16$week.trt %>% str_extract(., "[0-9]+") %>% as.numeric()
plotPars16$trt <- plotPars16$week.trt %>% str_extract(., "[aA-zZ]+") %>% factor()
plotPars16 <- plotPars16 %>% mutate(Trt = factor(ifelse(trt == "R", "Resistant", "Susceptible")))
plotPars16$Trt <- with(plotPars16, factor(Trt, levels(Trt)[c(2,1)]))


#### BW plot of CMM rates for paper
cmmPlot16 <- ggplot(data=plotPars16, aes(x=week, y=estimate, group = choice, shape = choice)) +
  ## Closed circles = infected source plant
  ## Open circles = Xf-free test plant
  geom_point(data = plotPars16, size=2, position = position_dodge(width = 0.9), colour = "black") +
  geom_errorbar(data = plotPars16, aes(ymax=ciu, ymin=cil), width=0.2, position = position_dodge(width = 0.9), colour = "black") +
  facet_grid(rate~Trt, scales = "free_y") +
  geom_hline(linetype = 2, yintercept = 0) +
  scale_shape_manual(values = c(16,1)) +
  scale_x_continuous(name = "", 
                     breaks = unique(plotPars16$week)) + 
  scale_y_continuous(name = "Rate (per hour)") +
  theme_bw(base_size=12) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") 

cmmPlot16

ggsave(file = "results/figures/2016_figures/pdr1_cmm_rate_parameter_plot_2016.jpg", plot = cmmPlot16,
       width = 14, height = 14, units = "in")


#####################################################################################
#### Plot raw numbers of bugs over time
plotCounts16 <- cmmData %>% mutate(xf_plant = n1/N, test_plant = n2/N) %>%
  dplyr::select(-n1, -n2, -n3, -t, -N) %>% 
  gather(., key = "choice", value = "proportion", xf_plant, test_plant)

## Change choice names to "infected" and "xf-free"
plotCounts16$choice <- ifelse(plotCounts16$choice == "xf_plant", "infected", "xf_free")

bgss_counts_plot16 <- ggplot(data=plotCounts16, aes(x=time_from_start_hr, y=proportion, group = choice)) +
  geom_line(aes(colour = choice), size=1.25) +
  facet_grid(trt~week, scales = "fixed") +
  scale_color_manual(values = my_palette) +
  scale_x_continuous(name = "Time from start (hr)") + 
  scale_y_continuous(name = "Proportion of vectors") +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank(),
        strip.background = element_blank()) 

bgss_counts_plot16

ggsave("results/figures/2016_figures/pdr1_bgss_counts_time-series_plot_2016.jpg", plot = bgss_counts_plot16,
       width = 14, height = 7, units = "in")




##########################################################################################################
#### Estimating rate parameters for each preference cage
#### To be used in SEM model for transmission
#### Only using the best estimates for each cage, ignoring variance

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

#### Load correct MLE results
modelFitsCage <- readRDS("output/CMM_optimx_model_selection_output_per_cage.rds")

#### Extract and organize parameter estimates from all models
## Note: profile likelihood to calculate variances fails numerous times. Ignore these because variances aren't used for per-cage analyses
paramResultsCage <- lapply(modelFitsCage, function(x) tryCatch(mleTable(x), error = function(e) NA))

#### Average results for all good models
paramAverageCage <- lapply(paramResultsCage, function(x) tryCatch(averageModels(x, dAIC.threshold = 7), error = function(e) NA))

paramAverageCage2 <- lapply(paramAverageCage, 
                            function(x) if(is.na(x)) data.frame("parameter" = c("p1", "p2", "mu1", "mu2"), 
                                                                "estimate" = rep(NA, 4),
                                                                "variance" = rep(NA, 4)) else x)

# Add week.trt combination to each element of the list and combine into one data.frame
for(i in 1:length(paramAverageCage2)){
  paramAverageCage2[[i]]$week.cage <- levels(cmmDataCage$week.cage)[i]
} 
paramDataCage2 <- rbindlist(paramAverageCage2) %>% as.data.frame()
paramDataCage2$estimate <- factor2numeric(paramDataCage2$estimate)
paramDataCage2$variance <- factor2numeric(paramDataCage2$variance)

saveRDS(paramDataCage, file = "output/CMM_rate_parameters_per_cage.rds")



###################################################################################################
###################################################################################################
#### 2017 experimental data
###################################################################################################

#### Importing data 
# Import from local .xlsx file
prefdata <- read.xlsx("data/2017_data/PdR1_2017_preference-transmission_experiment_data.xlsx", sheet = "BGSS_count_data", detectDates = TRUE)
str(prefdata)

# Separate leaf and symptom data from preference count data
leafdata <- prefdata[,15:ncol(prefdata)]
prefdata <- prefdata[,1:14]
prefdata$genotype <- factor(prefdata$genotype)
head(prefdata)


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

# Max N for any row of cmmData is 8*8 = 64
# Any rows of cmmData have more than N = 64 BGSS?
cmmData %>% dplyr::filter(N > 64)


#### Fit models to each of the week-genotype combinations
modelFits <- lapply(levels(cmmData$week.genotype), function(x) optimizeCMM(dat = cmmData[cmmData$week.genotype == x,], upperConstraint = 64, aiccN = 8))
names(modelFits) <- levels(cmmData$week.genotype)

## Load previously saved MLE results
#modelFits <- readRDS("output/CMM_optimx_model_selection_output_2017.rds")

#### Extract model selection tables
selectTables17 <- lapply(1:length(modelFits), function(x) modelFits[[x]]$modelSelect)
names(selectTables17) <- names(modelFits)

#### Calculate variance-covariance and correlation matrices for each trt-week combination
matrices <- lapply(modelFits, getParCorrelations)
names(matrices) <- levels(cmmData$week.genotype)

## Save vcov and corr matrices list
#saveRDS(matrices, file = "output/cmm_parameter_correlation_matrices_2017.rds")

#### Extract and organize parameter estimates from all models
paramResults <- lapply(modelFits, mleTable)
names(paramResults) <- levels(cmmData$week.genotype)


#### Average results for all good models
paramAverage <- lapply(paramResults, function(x) averageModels(x, dAIC.threshold = 7))

# Add week.genotype combination to each element of the list and combine into one data.frame
for(i in 1:length(paramAverage)){
  paramAverage[[i]]$week.genotype <- levels(cmmData$week.genotype)[i]
} 
paramData <- rbindlist(paramAverage) %>% as.data.frame()
paramData$estimate <- as.numeric(paramData$estimate)
paramData$variance <- as.numeric(paramData$variance)

#### Calculate SE and 95% CI
paramData$se <- sqrt(paramData$variance)
paramData$cil <- with(paramData, estimate - 1.96*se)
paramData$ciu <- with(paramData, estimate + 1.96*se)


#####################################################################################
#### Plotting rate parameters
# Structure data.frame for plotting
plotPars17 <- paramData
# Split week.trt into week and trt columns
week.genotype.split <- tstrsplit(plotPars17$week.genotype, split = "-")
plotPars17$week <- week.genotype.split[[1]] %>% as.numeric()
plotPars17$genotype <- week.genotype.split[[2]]
# Split parameter into rate and choice columns, and replace values with more meaningful terms
plotPars17$rate <- plotPars17$parameter %>% str_extract(., "[aA-zZ]+") %>% 
  gsub("p", "attraction", ., fixed = TRUE) %>% 
  gsub("mu", "leaving", ., fixed = TRUE)
plotPars17$choice <- plotPars17$parameter %>% str_extract(., "[0-9]+") %>% as.numeric() %>%
  gsub(1, "infected", ., fixed = TRUE) %>%
  gsub(2, "Xf-free", ., fixed = TRUE) %>% factor()
plotPars17$latticegroups <- with(plotPars17, paste(genotype, choice, sep=" ")) %>% factor()
plotPars17 <- plotPars17 %>% arrange(., week)
str(plotPars17)
head(plotPars17)

#write.csv(plotPars17, file = "results/pdr1_cmm_rate_parameter_estimates_2017.csv", row.names = FALSE)


#### Load previously saved data for plotting
#plotPars17 <- read.csv("results/pdr1_cmm_rate_parameter_estimates_2017.csv", header = TRUE)

# Prepend zeroes to genotype labels and make a factor
plotPars17$genotype <- plotPars17$genotype %>% formatC(., width = 3, format = "d", flag = "0") %>% factor()
str(plotPars17)


#### BW plot of CMM rates for paper
cmmPlot17 <- ggplot(data=plotPars17, aes(x=week, y=estimate, group = choice, shape = choice)) +
  ## Closed circles = infected source plant
  ## Open circles = Xf-free test plant
  geom_point(data = plotPars17, size=2, position = position_dodge(width = 0.9), colour = "black") +
  geom_errorbar(data = plotPars17, aes(ymax=ciu, ymin=cil), width=0.2, position = position_dodge(width = 0.9), colour = "black") +
  facet_grid(rate~genotype, scales = "free_y") +
  geom_hline(linetype = 2, yintercept = 0) +
  scale_shape_manual(values = c(16,1)) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = unique(plotPars17$week)) + 
  scale_y_continuous(name = "Rate (per hour)") +
  #limits = c(0,10)) +
  theme_bw(base_size=12) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") 

cmmPlot17

ggsave("results/figures/2017_figures/pdr1_cmm_rate_parameter_BW_plot_2017.jpg", plot = cmmPlot17,
       width = 14, height = 7, units = "in")



######################################################################################################
#### Plot raw counts of bugs over time
plotCounts17 <- cmmData %>% mutate(xf_plant = n1/N, test_plant = n2/N) %>%
  dplyr::select(-n1, -n2, -n3, -t, -N, -week.genotype) %>% 
  gather(., key = "choice", value = "proportion", xf_plant, test_plant)

## Change choice names to "infected" and "xf-free"
plotCounts17$choice <- ifelse(plotCounts17$choice == "xf_plant", "infected", "xf_free")

bgss_counts_plot17 <- ggplot(data=plotCounts17, aes(x=time_from_start_hr, y=proportion, group = choice)) +
  geom_line(aes(colour = choice), size=1.25) +
  #geom_point(aes(colour = choice), size=3.5, position = position_dodge(width = 0.9)) +
  #geom_errorbar(aes(ymax=ciu, ymin=cil), width=0.2, position = position_dodge(width = 0.9)) +
  facet_grid(genotype~week, scales = "fixed") +
  scale_color_manual(values = my_palette) +
  scale_x_continuous(name = "Time from start (hr)") + 
  scale_y_continuous(name = "Proportion of vectors") +
  #limits = c(0,10)) +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank(),
        strip.background = element_blank()) 

bgss_counts_plot17

ggsave("results/figures/2017_figures/pdr1_bgss_counts_time-series_plot_2017.jpg", plot = bgss_counts_plot17,
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

## Load previously saved MLE results
modelFitsCage <- readRDS("output/CMM_optimx_model_selection_output_per_cage.rds")

#### Extract and organize parameter estimates from all models
paramResultsCage <- lapply(modelFitsCage, function(x) tryCatch(mleTable(x), error = function(e) NA))

#### Average results for all good models
paramAverageCage <- lapply(paramResultsCage, function(x) tryCatch(averageModels(x, dAIC.threshold = 7), error = function(e) NA))

paramAverageCage2 <- lapply(paramAverageCage, 
                            function(x) if(is.na(x)) data.frame("parameter" = c("p1", "p2", "mu1", "mu2"), 
                                                                "estimate" = rep(NA, 4),
                                                                "variance" = rep(NA, 4)) else x)

# Add cage identifier to each element of the list and combine into one data.frame
for(i in 1:length(paramAverageCage2)){
  paramAverageCage2[[i]]$cage <- levels(cmmDataCage$cage)[i]
} 
paramDataCage2 <- rbindlist(paramAverageCage2) %>% as.data.frame()
paramDataCage2$estimate <- factor2numeric(paramDataCage2$estimate)
paramDataCage2$variance <- factor2numeric(paramDataCage2$variance)

saveRDS(paramDataCage, file = "output/CMM_2017_rate_parameters_per_cage.rds")


#####################################################################################################################
#### Saving plots and tables for both years for paper

#### Saving CMM rate plots from both years
cmm_rate_figure <- plot_grid(cmmPlot16, cmmPlot17,
                             ncol = 1, nrow = 2,
                             labels = c("(A)", "(B)"), label_size = 10)
cmm_rate_figure

ggsave(filename = "results/figures/cmm_rates_both_years_figure.tiff",
       plot = cmm_rate_figure,
       width = 10, height = 10, units = "cm", dpi = 600, compression = "lzw")


#### Saving model selection tables for both years

## Make the selections tables prettier
selectTables <- c(selectTables16, selectTables17)
for(i in 1:length(selectTables)){
  selectTables[[i]][,2:3] <- round(selectTables[[i]][2:3], digits = 2)
  selectTables[[i]]$model <- with(selectTables[[i]], ifelse(model == "choice", "Free Choice",
                                                            ifelse(model == "p.choice", "Free Attraction",
                                                                   ifelse(model == "mu.choice", "Free Leaving",
                                                                          "Fixed"))))
}

## Write all correlation matrices to a single Excel worksheet
wb <- createWorkbook()

#### model selection tables for both years
addWorksheet(wb, sheetName = "selection_tables")
startRows <- seq(2,2+6*length(selectTables),by=6)
for(i in 1:length(selectTables)){
  writeData(wb, sheet = "selection_tables", x = selectTables[[i]], startCol = 2, startRow = startRows[i])
  writeData(wb, sheet = "selection_tables", x = names(selectTables)[i], startCol = 1, startRow = startRows[i])
}

saveWorkbook(wb, file = "results/cmm_model_selection_tables.xlsx", overwrite = TRUE)



#### Time series of raw BGSS counts for both years
bgss_count_figure <- plot_grid(bgss_counts_plot16, bgss_counts_plot17,
                             ncol = 1, nrow = 2,
                             labels = c("a", "b"), label_size = 10)
bgss_count_figure

ggsave(filename = "results/figures/bgss_counts_both_years_figure.tiff",
       plot = bgss_count_figure,
       width = 14, height = 14, units = "cm", dpi = 300, compression = "lzw")
## Don't use this version. It's ugly. Just use the separate figures for the supplementary material



#### Comparing average movement rates between years
## 2016 rates
avgRates16 <- plotPars16 %>% group_by(rate) %>% summarise(mean = mean(estimate))
avgRates17 <- plotPars17 %>% group_by(rate) %>% summarise(mean = mean(estimate))                                          

