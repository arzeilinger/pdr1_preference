#### Analysis of PdR1 preference-transmission experiment 2017

rm(list = ls())
# Load packages
my.packages <- c("tidyr", "dplyr", "data.table", "openxlsx", "ggplot2",
                 "MASS", "logistf", "multcomp", "bbmle", "lme4", "googlesheets")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/factor2numeric.R")


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
# Things are complicated, latest re-culture wasn't always the best. Use complicated ifelse statements
# Make the final xf pop estimates
sourcedata2$xfpop <- with(sourcedata2, 
                          ifelse(culture_1 == 0 & (culture_2 == 0 | is.na(culture_2)) & (culture_3 == 0 | is.na(culture_3)), 0,
                                 ifelse((culture_1 == 0 | is.na(culture_1)) & (culture_2 > 0 | culture_3 > 0), 100,
                                        culture_1)))
print.data.frame(sourcedata2)

