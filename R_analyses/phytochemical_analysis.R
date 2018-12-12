#### PdR1 2017 PHYTO-CHEMICAL DATA
#### From Chris Wallis USDA-ARS

rm(list = ls())
# Load packages
my.packages <- c("tidyr", "dplyr", "data.table", "openxlsx", "ggplot2",
                 "MASS", "bbmle", "glmmLasso", "lmmen")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/factor2numeric.R")

# Import merged full phenolic data set
phenData <- readRDS("output/full_phenolic_data_pdr1_2017.rds")


#### Import original data files from Chris Wallis and merge
## Woody stem phenolics data set
stemPhenData <- read.xlsx("data/2017_data/BGSSandXfPhenolics.xlsx", sheet = "finalRawWoodPhenolics")
str(stemPhenData)
## Leaf phenolics data set
leafPhenData <- read.xlsx("data/2017_data/BGSSandXfPhenolics.xlsx", sheet = "finalLeafPhenolics")
str(leafPhenData)

# The "Number" column seems to match between leaf and stem data sets
stemPhenData %>% filter(Number == 200)
leafPhenData %>% filter(Number == 200)
nrow(stemPhenData)
nrow(leafPhenData)

#### Merge leaf and stem data sets
phenData <- full_join(stemPhenData, leafPhenData, by = c("Number", "Week", "Date", "Treatment", "Res", "Rep"), suffix = c(".stem", ".leaf"))
str(phenData)
# Save merged data set
saveRDS(phenData, "output/full_phenolic_data_pdr1_2017.rds")



