#### MUNGING, CLEANING, AND COMBINING DATA SETS 
#### FOR PDR1 PREFERENCE-TRANSMISSION EXPERIMENTS 2016 AND 2017

rm(list = ls())
# Load packages
my.packages <- c("tidyr", "dplyr", "data.table", "openxlsx", "MASS", "googlesheets")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/factor2numeric.R")


##############################################################################################################
##############################################################################################################
#### 2016 data
##############################################################################################################

#### Munging data

# Preference data
prefdata <- read.xlsx("data/2016_data/pdr1_preference_data.xlsx", sheet = "data")
str(prefdata)


# Separate leaf and symptom data from preference count data
leafdata <- prefdata[,c("week", "trt", "rep", "n_leaves_test", "n_leaves_source_start", 
                        "n_leaves_source_end", "n_pd_leaves", "n_ms_petioles", "pd_index")] %>% dplyr::filter(., !is.na(n_pd_leaves))
leafdata$n_leaves_test <- factor2numeric(leafdata$n_leaves_test)
leafdata$n_leaves_source_start <- factor2numeric(leafdata$n_leaves_source_start)
leafdata$n_leaves_source_end <- factor2numeric(leafdata$n_leaves_source_end)
# Some trials I only counted leaves at beginning, other trials only at end; combine total leaf data
leafdata$n_leaves_source <- leafdata %>% with(., cbind(n_leaves_source_start, n_leaves_source_end)) %>% rowMeans(., na.rm = TRUE)
leafdata$n_leaves_source[is.nan(leafdata$n_leaves_source)] <- NA
str(leafdata)


################################################################################################################
#### Constructing PD indices

#### PD index of Rashed et al. 2013:
# 0 = asymptomatic
# 1 = 1-2 scorched leaves
# 2 = 3-4 scorched leaves
# 3 = all leaves scorched; a few matchstick petioles
# 4 = all leaves heavily scorched and many matchstick petioles
# 5 = only a few leaves at the end of the cane

table(leafdata$n_pd_leaves, leafdata$n_ms_petioles)

#### Alternative PD index: 
prop_pd_leaves <- leafdata %>% with(., ifelse(n_leaves_source == 0, 1, n_pd_leaves/n_leaves_source))
leafdata$pd_index2 <- leafdata %>% with(., prop_pd_leaves + n_ms_petioles)
plot(x = leafdata$pd_index, y = leafdata$pd_index2)


####################################################################################################################
#### Import and combine culturing data, leaf data, and preference data
# Import culturing data
culturedata <- read.xlsx("data/2016_data/pdr1_culturing_data.xlsx", sheet = "Infection data")
# Remove repeat culturing data
culturedata <- culturedata %>% dplyr::filter(., notes != "repeat culture" | is.na(notes))
str(culturedata)
culturedata

# Import preference data: estimates of rate parameters from CM model for each cage
# choice 1 = source plants, choice 2 = test plants
# Ignoring variance around parameters.
# TO DO: check on convergence issues
paramDataCage <- readRDS("output/CMM_rate_parameters_per_cage.rds")
# Reshape paramDataCage
paramDataCage <- paramDataCage[,c("parameter", "estimate", "week.cage")] %>% spread(., key = parameter, value = estimate)

# Merge culturing and leaf data on week-trt-rep combination
transdata <- inner_join(culturedata, leafdata, by = c("week", "trt", "rep")) 

# Merge with preference data
transdata$week.cage <- transdata %>% with(., paste(week, trt, rep, sep=""))
transdata <- left_join(transdata, paramDataCage, by = "week.cage")
str(transdata)

saveRDS(transdata, file = "output/pdr1_transmission_preference_dataset.rds")


#############################################################################################################
#### Import vector acquisition data and combine with transmission-preference dataset
acqData <- readRDS("output/pdr1_2016_vector_cfu_from_qpcr.rds")
# Remove qPCR serial dilution and NTC rows
acqData <- acqData %>% dplyr::filter(!is.na(insect_code))
str(acqData)
# Remove one crazy outlier!
acqData <- acqData %>% dplyr::filter(!(tube_code == 117 & wellNumber == 60))
# Simplify data set
acqDataVector <- acqData %>% dplyr::select(week, trt, rep, tube_code, cfu) %>% group_by(week, trt, rep, tube_code) %>% summarise(vectorcfu = mean(cfu))
acqDataVector$week <- as.numeric(acqDataVector$week)
acqDataVector$rep <- as.numeric(acqDataVector$rep)
str(acqDataVector)

#### Constructing vector infection index
# In each trial, I have Xf pops for each of the vectors. I could include in transmission model:
# total Xf population among all vectors
# proportion of vectors infected (b/c some cages have fewer than 8 vectors)
# an evenness index of Xf pops among vectors
# Matt's transmission parameters paper might have something to say about this
# multiple possibilities can be evaluated using AIC

## Average duplicates for each sample
acqDataCage <- acqDataVector %>% group_by(week, trt, rep) %>% summarise(cagecfu = mean(vectorcfu),
                                                                        sdcfu = sd(vectorcfu),
                                                                        logCagecfu = mean(log10(vectorcfu + 1)),
                                                                        propVectorInfected = sum(vectorcfu > 0, na.rm = TRUE)/length(vectorcfu[!is.na(vectorcfu)]))
print.data.frame(acqDataCage)

## Merge with acquisition data at cage level with transmission-preference data set
acqDataCage$week.cage <- with(acqDataCage, paste(week, trt, rep, sep=""))

transdata <- readRDS("output/pdr1_transmission_preference_dataset.rds")

transdata <- left_join(transdata, acqDataCage, by = c("week.cage", "week", "trt", "rep"))
transdata

#### Save final data set, including culturing/transmission/symptom data, CMM parameter estimates, and vector acquisition/qPCR data
saveRDS(transdata, file = "output/pdr1_transmission_preference_dataset.rds")

