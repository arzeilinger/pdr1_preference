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
acqDataVector <- acqData %>% dplyr::select(week, trt, rep, tube_code, cfu, vectorInfected) %>% 
  group_by(week, trt, rep, tube_code) %>% 
  summarise(vectorcfu = mean(cfu),
            vectorInfectious = ifelse(any(vectorInfected == 1), 1, 0))
acqDataVector$week <- as.numeric(acqDataVector$week)
acqDataVector$rep <- as.numeric(acqDataVector$rep)
acqDataVector$week.cage <- factor(with(acqDataVector, paste(week, trt, rep, sep = "")))
acqDataVector$vectorcfu <- as.integer(acqDataVector$vectorcfu)
str(acqDataVector)

## Save acquisition data at the vector level for acquisition analysis
saveRDS(acqDataVector, file = "output/pdr1_2016_vector_acquisition_dataset.rds")


#### Cage-level acquisition data
## Average duplicates for each sample
acqDataCage <- acqDataVector %>% group_by(week, trt, rep) %>% summarise(cagecfu = mean(vectorcfu, na.rm = TRUE),
                                                                        sdcfu = sd(vectorcfu, na.rm = TRUE),
                                                                        logCagecfu = mean(log10(vectorcfu + 1), na.rm = TRUE),
                                                                        totalVectorInfectious = sum(vectorInfectious, na.rm = TRUE),
                                                                        propVectorInfectious = totalVectorInfectious/length(vectorInfectious[!is.na(vectorInfectious)]))
print.data.frame(acqDataCage)

## Merge acquisition data at cage level with transmission-preference data set
acqDataCage$week.cage <- with(acqDataCage, paste(week, trt, rep, sep=""))

transdata <- readRDS("output/pdr1_transmission_preference_dataset.rds")

transdata <- left_join(transdata, acqDataCage, by = c("week.cage", "week", "trt", "rep"))
transdata



#### Model selection on PD symptom index
pdMod1 <- glm(test.plant.infection ~ week*trt + log10(source.cfu.per.g+1) + p1 + p2 + mu1 + mu2 + pd_index, data = transdata, family = "binomial")
pdMod2 <- glm(test.plant.infection ~ week*trt + log10(source.cfu.per.g+1) + p1 + p2 + mu1 + mu2 + pd_index2, data = transdata, family = "binomial")
# I think quasibinomial distribution might be better but then it doesn't calculate an AIC value. Need to look into this.
AICtab(pdMod1, pdMod2, base = TRUE)
plot(pdMod2)
summary(pdMod2)
# PD indices are essentially the same; go with Arash's index
## Remove pd_index2 column; pd_index shows similar results and was the same used by Arash/Krivanek papers
transdata <- transdata %>% dplyr::select(-pd_index2)


#### Save final data set, including culturing/transmission/symptom data, CMM parameter estimates, and vector acquisition/qPCR data
saveRDS(transdata, file = "output/pdr1_transmission_preference_dataset.rds")



##############################################################################################################
##############################################################################################################
#### 2017 data
##############################################################################################################
#### Importing and munging data
# Importing from local .xlsx file
# transdata <- read.xlsx("data/2017_data/PdR1_2017_preference-transmission_experiment_data.xlsx", sheet = "test_plant_culturing", detectDates = TRUE)
# Import from Googlesheets
pdr1DataURL <- gs_url("https://docs.google.com/spreadsheets/d/14uJLfRL6mPrdf4qABeGeip5ZkryXmMKkan3mJHeK13k/edit?usp=sharing",
                      visibility = "private")
transdataGS <- gs_read(pdr1DataURL, ws = "test_plant_culturing")
transdata17 <- transdataGS

# Remove Control plant samples
transdata17 <- transdata17 %>% dplyr::filter(!grepl("CTRL", genotype))

# Combine test_plant_infection_1 and test_plant_infection_2 columns
# In all cases, when I re-cultured a sample (test_plant_infection_2), the results were the same as the 1st time or more reliable
# So go with test_plant_infection_2 results when they are available
transdata17$test_plant_infection <- with(transdata17, ifelse(!is.na(test_plant_infection_2), test_plant_infection_2, test_plant_infection_1)) %>% as.integer()                                      

# Fix column classes
transdata17$genotype <- factor(transdata17$genotype)
transdata17$trt <- factor(transdata17$trt)
transdata17$block <- factor(transdata17$block)

str(transdata17)
summary(transdata17)


################################################################################################################
#### Join qPCR data and transmission data

#### Import qPCR vector acquisition data
## Calling object vectorData3 to keep with namings from "qpcr_calculations_pdr1_2017.R" script
vectorData3 <- readRDS("output/vector_acquisition_data_pdr1_2017.rds")

transVectorData <- transdata17 %>% dplyr::select(week, block, genotype, trt, rep, PD_symptoms_index, test_plant_infection, notes) %>%
  left_join(., vectorData3, by = c("genotype", "trt", "rep", "week", "block"))

## Comparing transmission and acquisition data
transVectorData %>% dplyr::filter(!is.na(propInfectious)) %>% arrange(week, block) %>% print.data.frame 
transVectorData %>% arrange(week, block) %>% tail(n = 30) %>% print.data.frame()

# Note: Trial 8-2-094R4 has a positive plant and all 8 bugs were negative. Need to double check this
# Also, 14-1-102R4 has a positive plant and 0 positive bugs, but only 1 bug tested.

saveRDS(transVectorData, file = "output/pdr1_transmission_preference_dataset_2017.rds")
str(transVectorData)

## Quick plot of proportion infectious over time
qplot(x = week, y = propInfectious, data = transVectorData)

## Calculate propInfected for test plants and mean propInfectious for vectors for each week-genotype combination
transVectorSummary <- transVectorData %>% dplyr::filter(!is.na(propInfectious)) %>% 
  group_by(week, genotype, trt) %>% summarise(propPlantInfected = sum(test_plant_infection)/length(test_plant_infection),
                                              meanPropInfectious = mean(propInfectious, na.rm = TRUE),
                                              nPropInfectious = length(propInfectious),
                                              sePropInfectious = sd(propInfectious, na.rm = TRUE)/sqrt(nPropInfectious)) 
transVectorSummary



################################################################################################################
#### Join source plant culturing data to transmission/acquisition data

#### Import culturing data
# Import from local .xlsx file
#sourcedata <- read.xlsx("data/2017_data/PdR1_2017_preference-transmission_experiment_data.xlsx", sheet = "source_plant_culturing", detectDates = TRUE)

# Import from Googlesheets
pdr1DataURL <- gs_url("https://docs.google.com/spreadsheets/d/14uJLfRL6mPrdf4qABeGeip5ZkryXmMKkan3mJHeK13k/edit?usp=sharing",
                      visibility = "private")
sourcedataGS <- gs_read(pdr1DataURL, ws = "source_plant_culturing")
sourcedata <- sourcedataGS
str(sourcedata)
summary(sourcedata)

# Convert classes for columns
sourcedata$block <- factor(sourcedata$block)
sourcedata$genotype <- factor(sourcedata$genotype)
sourcedata$trt <- factor(sourcedata$trt)
sourcedata$xf_cfu_per_g <- as.numeric(sourcedata$xf_cfu_per_g)

# Remove NAs and Control samples
print.data.frame(sourcedata[is.na(sourcedata$xf_cfu_per_g), ])
sourcedata <- sourcedata %>% dplyr::filter(!is.na(xf_cfu_per_g))


# #### Which plants were negative at least once?
# # make a plant ID
# sourcedata$plantID <- with(sourcedata, paste(week, block, genotype, trt, rep, sep = "-"))
# 
# negativeSamples <- sourcedata %>% filter(., xf_cfu_d0 == 0 & xf_cfu_d1 == 0 & xf_cfu_d2 == 0 & genotype != "CAB-") %>%
#   dplyr::select(., plantID, times_cultured, starts_with("xf_cfu"), notes) %>% arrange(., plantID, times_cultured)
# 
# # get all times cultured for each plant that was negative once
# negativePlants <- sourcedata[sourcedata$plantID %in% negativeSamples$plantID,] %>% 
#   dplyr::select(., plantID, times_cultured, starts_with("xf_cfu"), notes) %>% arrange(., plantID, times_cultured)


#### For plants that were cultured multiple times, get the best estimate of Xf pop
# Remove unnecessary columns and use spread() to get each plant on a single row
sourcedata <- sourcedata %>% mutate(plantID = factor(paste(week, block, genotype, rep, sep = "-"))) %>%
  mutate(times_cultured = paste("culture", times_cultured, sep = "_"))
sourcedata2 <- sourcedata %>% dplyr::select(-xf_plant_sample_mass, -xf_cfu_d0, -xf_cfu_d1, -xf_cfu_d2, -xf_plant_date_cultured, -notes, -notes2)
sourcedata2 <- sourcedata2 %>% spread(., key = times_cultured, value = xf_cfu_per_g)
print.data.frame(sourcedata2)
length(unique(sourcedata2$plantID)) == nrow(sourcedata2) # Check that spread() correctly prdouced a unique row for each plant ID; if so should be "TRUE"
# Make sure that the latest re-culture was always the best
recultures <- sourcedata2 %>% dplyr::filter(!is.na(culture_2))
recultures
# Six plants were initially negative or NA and were positive upon re-culturing
# Things are complicated, latest re-culture wasn't always the best. Use complicated ifelse statements
# Make the final xf pop estimates
# IMPORTANT: If plant initially tested negative (culture_1 = 0) but later culturing tested positive, I set xfpop = 100, 
# xfpop = 100 is an arbitrary value below the threshold of detection (which is around 700 - 800).
sourcedata2$xfpop <- with(sourcedata2, 
                          ifelse(culture_1 == 0 & (culture_2 == 0 | is.na(culture_2)) & (culture_3 == 0 | is.na(culture_3)), 0,
                                 ifelse((culture_1 == 0 | is.na(culture_1)) & (culture_2 > 0 | culture_3 > 0), 100,
                                        culture_1)))
print.data.frame(sourcedata2)

#### Which trials are missing from the final data set?
with(sourcedata2, table(week, genotype))
## Missing: 14-1-092S-1 and -2, 14-1-094R-4
## All were positive but highly contaminated or too dense to count on first culturing, not cultured again.

#### Some summary stats on source plant infections
# How many plants were "true" negatives?
sum(sourcedata2$xfpop == 0)
# 10 plants were "true" negatives
sourcedata2 %>% dplyr::filter(xfpop == 0)
# All true negative plants are Resistant lines; 9/10 are 102R


#### Load cleaned transmission/acquisition data set
transVectorData <- readRDS("output/pdr1_transmission_preference_dataset_2017.rds")
## Join with source culturing data, transVCData = trans, Vector, Culturing data
transVCData <- left_join(transVectorData, sourcedata2, by = c("week", "block", "genotype", "trt", "rep")) %>%
  dplyr::select(-culture_1, -culture_2, -culture_3)

with(transVectorData, table(week, genotype))
str(transVCData)


#######################################################################################################################
#### Join preference data to transmission/acquisition/culturing data

#### Read in preference rate estimates on per cage level
## choice1 = source (Xylella infected) plant
## choice2 = test plant
paramDataCage <- readRDS("output/CMM_2017_rate_parameters_per_cage.rds")

## Strip cage identifier to multiple columns
cage.strip <- tstrsplit(paramDataCage$cage, split = "-")
paramDataCage$week <- cage.strip[[1]]
paramDataCage$block <- factor(cage.strip[[2]])
paramDataCage$genotype <- cage.strip[[3]]
paramDataCage$trt <- factor(cage.strip[[4]])
paramDataCage$rep <- cage.strip[[5]]
## Add another "Rep" column, corresponding to "Rep" in phenolic dataset
paramDataCage$Rep2 <- with(paramDataCage, ifelse(block == 1, rep, 
                                                 ifelse(rep == 1, 5,
                                                        ifelse(rep == 2, 6,
                                                               ifelse(rep == 3, 7, 
                                                                      ifelse(rep == 4, 8, NA)))))) %>% as.numeric()
## Spread dataset so that each row is a cage; need to drop variances first
paramDataCage2 <- paramDataCage %>% dplyr::select(-variance, -cage) %>% spread(parameter, estimate)
paramDataCage2$week <- as.numeric(paramDataCage2$week)
paramDataCage2$rep <- as.numeric(paramDataCage2$rep)

str(paramDataCage2)
with(paramDataCage2, table(week, genotype))
## All trials are accounted for

#### Join with transmission/acquisition/culturing data set
transVCPdata <- left_join(transVCData, paramDataCage2, by = c("week", "block", "genotype", "trt", "rep"))
transVCPdata$genotype <- factor(transVCPdata$genotype)
str(transVCPdata)
summary(transVCPdata)

saveRDS(transVCPdata, file = "output/complete_2017_transmission-preference_dataset.rds")



#######################################################################################################################
#### Join phenolic/chemistry data to transmission/acquisition/culturing/preference data

#### Load transmission/acquisition/culturing/preference data set if not already loaded
transVCPdata <- readRDS("output/complete_2017_transmission-preference_dataset.rds")

## Filter phenolic data set to only Treatment == "Both" as these correspond to the xf_plants or source plants in the trials
## Chemistry data includes phenolics and volatiles, but continue to call it "phenData"
## Chemistry data does not include volatiles from stems because there were too many missing observations
phenData <- readRDS("output/full_phenolic_and_volatile_data_pdr1_2017.rds")
phenData <- phenData %>% dplyr::filter(Treatment == "Both")
phenData$Rep <- as.numeric(phenData$Rep)
# Make "Res" groups capitalized
phenData$Res <- with(phenData, ifelse(Res == "r", "R", "S"))
# Add a genotype column for merging with other data sets
phenData$genotype <- with(phenData, ifelse(Res == "R", "094", "092"))


#### Merge CMM preference parameter data set, phenolic data set, and transmission data set
phenPrefTransData <- transVCPdata %>% 
  left_join(., phenData, by = c("week" = "Week", "genotype", "trt" = "Res", "Rep2" = "Rep")) %>% 
  arrange(week)
str(phenPrefTransData)
summary(phenPrefTransData)

## A lot of NAs in chemical data were produced; what's going on?
## Because genotypes 007 and 102 are still in the data set
## Remove genotypes 007 and 102 because we never measured phenolics from them
phenPrefTransData <- phenPrefTransData %>% dplyr::filter(genotype == "092" | genotype == "094")
with(phenPrefTransData, table(genotype, week))
summary(phenPrefTransData)

phenPrefTransData %>% dplyr::filter(is.na(total.phenolics.stem)) %>% dplyr::select(1:5) %>% print.data.frame
phenPrefTransData %>% dplyr::filter(is.na(total.volatiles)) %>% dplyr::select(1:5) %>% print.data.frame
## week 2 block 1 was never measured, 3 other points missing from phenolic data and 5 others missing from volatile data

## Check missing reps again
with(phenPrefTransData, table(rep, week, genotype))

## Save final data set
saveRDS(phenPrefTransData, file = "output/chemistry_preference_transmission_dataset.rds")
