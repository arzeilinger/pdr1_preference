#### Analysis of PdR1 preference-transmission experiment 2017

rm(list = ls())
# Load packages
my.packages <- c("tidyr", "dplyr", "data.table", "openxlsx", "ggplot2",
                 "MASS", "logistf", "multcomp", "bbmle", "lme4")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/factor2numeric.R")

# Function to look at whole df_tbl object
printTibble <- function(dftbl){
  return(print(dftbl, n = nrow(dftbl)))
}

#### Import culturing data
sourcedata <- read.xlsx("data/2017_data/PdR1_2017_preference-transmission_experiment_data.xlsx", sheet = "source_plant_culturing", detectDates = TRUE)
str(sourcedata)
summary(sourcedata)

# Remove 

# Convert classes for columns
sourcedata$block <- factor(sourcedata$block)
sourcedata$genotype <- factor(sourcedata$genotype)
sourcedata$trt <- factor(sourcedata$trt)
sourcedata$xf_cfu_per_g <- as.numeric(sourcedata$xf_cfu_per_g)

sourcedata[is.na(sourcedata$xf_cfu_per_g),]

# Which plants were negative at least once?
# make a plant ID
sourcedata$plantID <- with(sourcedata, paste(week, block, genotype, trt, rep, sep = "-"))

negativeSamples <- sourcedata %>% filter(., xf_cfu_d0 == 0 & xf_cfu_d1 == 0 & xf_cfu_d2 == 0 & genotype != "CAB-") %>%
  dplyr::select(., plantID, times_cultured, starts_with("xf_cfu"), notes) %>% arrange(., plantID, times_cultured)

# get all times cultured for each plant that was negative once
negativePlants <- sourcedata[sourcedata$plantID %in% negativeSamples$plantID,] %>% 
  dplyr::select(., plantID, times_cultured, starts_with("xf_cfu"), notes) %>% arrange(., plantID, times_cultured)


