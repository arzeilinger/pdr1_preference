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
                 "googlesheets", "RColorBrewer", "cowplot")
lapply(my.packages, require, character.only = TRUE)

## Input functions for P1 and P2 equations, 
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
#write.csv(prefdata, "output/pdr1_preference_2016_data.csv", row.names = FALSE)
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

r12Data <- dplyr::filter(cmmData, week == 12, trt == "R") 
dplyr::filter(cmmData, week == 3)


#### Fit models to each of the trt-week combinations
modelFits <- lapply(levels(cmmData$week.trt), function(x) optimizeCMM(dat = cmmData[cmmData$week.trt == x,], upperConstraint = 64, aiccN = 8))
names(modelFits) <- levels(cmmData$week.trt)
saveRDS(modelFits, file = "output/CMM_optimx_model_selection_output_2016.rds")

modelFits <- readRDS("output/CMM_optimx_model_selection_output_2016.rds")

#### Extract model selection tables
selectTables16 <- lapply(1:length(modelFits), function(x) modelFits[[x]]$modelSelect)
names(selectTables16) <- names(modelFits)


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

# dplyr::filter(paramData, week.trt == "3S")
# 
# dplyr::filter(paramData, week.trt == "12S" | week.trt == "12R")
# testpardat <- dplyr::filter(paramData, week.trt == "12R")
# testpardat


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
plotPars16

write.csv(plotPars16, file = "results/pdr1_cmm_rate_parameter_estimates_2016.csv", row.names = FALSE)

# # Create dummy x variable to space out points
# adj <- c(-0.5, -0.25, 0.25, 0.5)
# plotPars$dummyx <- 0
# for(i in 1:length(levels(plotPars$latticegroups))){
#   group.i <- levels(plotPars$latticegroups)[i]
#   data.i <- which(plotPars$latticegroups == group.i)
#   adj.i <- adj[i]
#   plotPars$dummyx[data.i] <- plotPars[plotPars$latticegroups == group.i,"week"] + adj.i
# }


# Rate parameter plots together
#tiff(filename = "results/figures/pdr1_cmm_rate_parameter_plot.tif")
     # width = 76*2, height = 76, units = "mm", 
     # res = 600, compression = "lzw")
#  plot.new()
# parameter_plot <-  with(plotPars,
#                         xyplot(estimate ~ dummyx|rate, groups = latticegroups,
#                                ly = cil, uy = ciu,
#                                scales = list(col = 1, alternating = 1, tck = c(1, 0), cex = 1.1, relation = "free",
#                                              x = list(limits = c(0, 13), at = c(3,8,12),
#                                                       labels = list(c("","", ""), c(3,8,12))),
#                                              y = list(limits = list(c(0, 3), c(0, 0.6)),
#                                                       at = list(seq(0, 3, by = 1), seq(0, 0.6, by = 0.2)),
#                                                       labels = list(seq(0, 3, by = 1), seq(0, 0.6, by = 0.2)))),
#                                xlab = list("Weeks post-inoculation", cex = 1.2), 
#                                ylab = list(expression("Leaving rate "(hr^{-1})*"                          Attraction rate "(hr^{-1})),
#                                            cex = 1.2),
#                                layout = c(1,2), as.table = TRUE, strip = FALSE, pch = c(19, 1, 17, 2),
#                                type = 'p', cex = 1.2, col = "black", 
#                                key = list(x = 0.4, y = 0.85, corner = c(0,0),
#                                           text = list(lab = levels(latticegroups)), 
#                                           points = list(pch = c(19, 1, 17, 2), col = "black")), 
#                                prepanel = prepanel.ci,                      
#                                panel = function(x, y, ...) {                
#                                  panel.abline(v = unique(as.numeric(x)),  
#                                               col = "white")              
#                                  panel.superpose(x, y, ...)               
#                                },                                          
#                                panel.groups = panel.ci))                    
#                 # ltext(165, 9, "2009", cex = 1.3)
#                 # ltext(380, 9, "2010", cex = 1.3)
#                 # mtext(c("A", "B"), side = 3, cex = 1.3, adj = rep(0.06, 2), padj = c(-1.4, 15))
#                 # mtext(c("C", "D"), side = 3, cex = 1.3, adj = rep(0.64, 2), padj = c(-1.4, 15))
# 
# trellis.device(device = "tiff", file = "results/figures/pdr1_cmm_rate_parameter_plot_test.tif")
#   print(parameter_plot)
# dev.off()

#### Load rate parameter data for plotting
plotPars16 <- read.csv("results/pdr1_cmm_rate_parameter_estimates_2016.csv") %>% arrange(week)
plotPars16 <- plotPars16 %>% mutate(Trt = factor(ifelse(trt == "R", "Resistant", "Susceptible")))
plotPars16$Trt <- with(plotPars16, factor(Trt, levels(Trt)[c(2,1)]))


# Color palette
my_palette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
my_palette = c(brewer.pal(5, "Blues")[c(4)], brewer.pal(5, "Set1")[c(3)])

#### Color Plot
cmmPlot16_color <- ggplot(data=plotPars16, aes(x=week, y=estimate, group = choice, color = choice)) +
  #geom_line(aes(linetype=inoc.time, colour = trt), size=1.25) +
  geom_point(data = plotPars16, size=3.5, position = position_dodge(width = 0.9)) +
  geom_errorbar(data = plotPars16, aes(ymax=ciu, ymin=cil), width=0.2, position = position_dodge(width = 0.9)) +
  facet_grid(rate~Trt, scales = "free_y") +
  geom_hline(linetype = 2, yintercept = 0) +
  scale_color_manual(values = my_palette) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = unique(plotPars16$week)) + 
  scale_y_continuous(name = "Rate (per hour)") +
  #limits = c(0,10)) +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank(),
        strip.background = element_blank()) 

cmmPlot16_color

ggsave("results/figures/2016_figures/pdr1_cmm_rate_parameter_plot_2016.jpg", plot = cmmPlot16_color,
       width = 14, height = 7, units = "in")


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
  #limits = c(0,10)) +
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
  #geom_point(aes(colour = choice), size=3.5, position = position_dodge(width = 0.9)) +
  #geom_errorbar(aes(ymax=ciu, ymin=cil), width=0.2, position = position_dodge(width = 0.9)) +
  facet_grid(trt~week, scales = "fixed") +
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

bgss_counts_plot16

ggsave("results/figures/2016_figures/pdr1_bgss_counts_time-series_plot_2016.jpg", plot = bgss_counts_plot16,
       width = 14, height = 7, units = "in")




##########################################################################################################
#### Estimating rate parameters for each preference cage
#### To be used in logistic regression model for transmission
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
saveRDS(modelFitsCage, file = "output/CMM_optimx_model_selection_output_per_cage.rds")

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
#write.csv(prefdata, file = "output/pdr1_preference_2017_data.csv", row.names = FALSE)
# Import from Googlesheets
# pdr1DataURL <- gs_url("https://docs.google.com/spreadsheets/d/14uJLfRL6mPrdf4qABeGeip5ZkryXmMKkan3mJHeK13k/edit?usp=sharing",
#                       visibility = "private")
# prefdataGS <- gs_read(pdr1DataURL, ws = "BGSS_count_data")
# prefdata <- prefdataGS
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
data_14_102 <- dplyr::filter(cmmData, week == 14, genotype == "102") 
data_14_094 <- dplyr::filter(cmmData, week == 14, genotype == "094") 
data_14_092 <- dplyr::filter(cmmData, week == 14, genotype == "092") 

# Max N for any row of cmmData is 8*8 = 64
# Any rows of cmmData have more than N = 64 BGSS?
cmmData %>% dplyr::filter(N > 64)
cmmData %>% filter(week == 2 & genotype == "007")


#### Fit models to each of the week-genotype combinations
modelFits <- lapply(levels(cmmData$week.genotype), function(x) optimizeCMM(dat = cmmData[cmmData$week.genotype == x,], upperConstraint = 64, aiccN = 8))
names(modelFits) <- levels(cmmData$week.genotype)
saveRDS(modelFits, file = "output/CMM_optimx_model_selection_output_2017.rds")

modelFits <- readRDS("output/CMM_optimx_model_selection_output_2017.rds")

#### Extract model selection tables
selectTables17 <- lapply(1:length(modelFits), function(x) modelFits[[x]]$modelSelect)
names(selectTables17) <- names(modelFits)

#### Calculate variance-covariance and correlation matrices for each trt-week combination
matrices <- lapply(modelFits, getParCorrelations)
names(matrices) <- levels(cmmData$week.genotype)

## Save vcov and corr matrices list
saveRDS(matrices, file = "output/cmm_parameter_correlation_matrices_2017.rds")

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
paramData$estimate <- factor2numeric(paramData$estimate)
paramData$variance <- factor2numeric(paramData$variance)

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
plotPars17$genotype <- week.genotype.split[[2]] %>% factor()
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

write.csv(plotPars17, file = "results/pdr1_cmm_rate_parameter_estimates_2017.csv", row.names = FALSE)


#### Plotting using ggplot2
plotPars17 <- read.csv("results/pdr1_cmm_rate_parameter_estimates_2017.csv", header = TRUE)

# Prepend zeroes to genotype labels and make a factor
plotPars17$genotype <- plotPars17$genotype %>% formatC(., width = 3, format = "d", flag = "0") %>% factor()

# Color palette
my_palette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
my_palette = c(brewer.pal(5, "Blues")[c(4)], brewer.pal(5, "Set1")[c(3)])

#### Color Plot
parameter_plot2 <- ggplot(data=plotPars17, aes(x=week, y=estimate, group = choice, color = choice)) +
  #geom_line(aes(linetype=inoc.time, colour = trt), size=1.25) +
  geom_point(data = plotPars17, size=3.5, position = position_dodge(width = 0.9)) +
  geom_errorbar(data = plotPars17, aes(ymax=ciu, ymin=cil), width=0.2, position = position_dodge(width = 0.9)) +
  facet_grid(rate~genotype, scales = "free_y") +
  geom_hline(linetype = 2, yintercept = 0) +
  scale_color_manual(values = my_palette) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = unique(plotPars17$week)) + 
  scale_y_continuous(name = "Rate (per hour)") +
                     #limits = c(0,10)) +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank(),
        strip.background = element_blank()) 

parameter_plot2

ggsave("results/figures/2017_figures/pdr1_cmm_rate_parameter_plot_2017.jpg", plot = parameter_plot2,
       width = 14, height = 7, units = "in")


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

#### Example of a single BGSS count plot
plotCounts1 <- plotCounts %>% dplyr::filter(week == 2 & genotype == "092")

bgss_counts_plot1 <- ggplot(data=plotCounts1, aes(x=time_from_start_hr, y=proportion, group = choice)) +
  geom_line(aes(colour = choice), size=1.25) +
  geom_point(aes(colour = choice), size=3.5) +
  scale_color_manual(values = my_palette) +
  scale_x_continuous(name = "Time from start (hr)") + 
  scale_y_continuous(name = "Proportion of BGSS", limits = c(0,1)) +
  scale_shape_manual(values = 1) +
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

modelFitsCage <- readRDS("output/CMM_optimx_model_selection_output_per_cage.rds")

#### Extract and organize parameter estimates from all models
paramResultsCage <- lapply(modelFitsCage, function(x) tryCatch(mleTable(x), error = function(e) NA))
#paramResultsCage <- lapply(modelFitsCage, mleTable)

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
       width = 15, height = 15, units = "cm", dpi = 600, compression = "lzw")


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
## Don't use this version. It's ugly. Just use the separate figures.



#### Comparing average movement rates between years
## 2016 rates
avgRates16 <- plotPars16 %>% group_by(rate) %>% summarise(mean = mean(estimate))
avgRates17 <- plotPars17 %>% group_by(rate) %>% summarise(mean = mean(estimate))                                          

