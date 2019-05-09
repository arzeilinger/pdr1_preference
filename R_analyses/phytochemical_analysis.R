#### PdR1 2017 PHYTO-CHEMICAL DATA
#### From Chris Wallis USDA-ARS

rm(list = ls())
# Load packages
my.packages <- c("tidyr", "dplyr", "data.table", "openxlsx", "ggplot2")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/factor2numeric.R")
source("R_functions/source_https.R")

## Source number2words function directly from GitHub
## The script includes a few examples, so returns those examples
source_https("https://gist.githubusercontent.com/EconometricsBySimulation/3282defc1557fd0614f2/raw/2c1a502b032aafa01c88b24a61b278f75c3e1333/numbers2words.R")


#### Import original data files from Chris Wallis and merge
## Woody stem phenolics data set
stemPhenData <- read.xlsx("data/2017_data/XylellaBGSS_final_phenolics_2019-04-19.xlsx", sheet = "Wood", detectDates = TRUE)
str(stemPhenData)
## Leaf phenolics data set
leafPhenData <- read.xlsx("data/2017_data/XylellaBGSS_final_phenolics_2019-04-19.xlsx", sheet = "Leaf", detectDates = TRUE)
str(leafPhenData)

# The "Number" column seems to match between leaf and stem data sets
stemPhenData %>% filter(Number == 200) %>% dplyr::select(1:9)
leafPhenData %>% filter(Number == 200) %>% dplyr::select(1:9)
nrow(stemPhenData)
nrow(leafPhenData)

#### Merge leaf and stem data sets
## Need to change "BothPost" Treatment level in stemPhenData to "Both"
table(stemPhenData$Treatment)
stemPhenData$Treatment[stemPhenData$Treatment == "BothPost"] <- "Both"
table(stemPhenData$Treatment)

phenData <- full_join(stemPhenData, leafPhenData, by = c("Number", "Week", "Date", "Treatment", "Res", "Rep"), suffix = c(".stem", ".leaf"))
## Need to fix some of the Treatment codes to correctly join data sets
str(phenData)
# Save merged data set
saveRDS(phenData, "output/full_phenolic_data_pdr1_2017.rds")


###############################################################################################
#### Import and clean up volatile data

#### Load compiled and cleaned phenolic data
phenData <- readRDS("output/full_phenolic_data_pdr1_2017.rds")
str(phenData)
with(phenData, table(Week, Res, Treatment))



#### Load foliar and stem volatile data
leafvolData <- read.xlsx("data/2017_data/Xylella_BGSS_Study_Volatiles_2019-05-09.xlsx", sheet = "Foliar Volatiles")
stemvolData <- read.xlsx("data/2017_data/Xylella_BGSS_Study_Volatiles_2019-05-09.xlsx", sheet = "Wood Volatiles")
str(leafvolData)  
## Chris said that some observations are missing for stem volatiles because there wasn't enough tissue
summary(stemvolData) ## Doesn't show NAs, observations were probably not entered into spreadsheet because all NA
with(stemvolData, table(Week, Res, Treatment))
with(leafvolData, table(Week, Res, Treatment))
## Compare column names
cbind(names(phenData)[1:10], names(leafvolData)[1:10], names(stemvolData)[1:10])
## Looks like the important columns are all the same


#### Merge phenolic and volatile data
## Don't include stem volatile data because of too many missing values
chemData <- full_join(phenData, leafvolData, by = c("Number", "Date", "Week", "Treatment", "Res", "Rep"))
with(chemData, table(Week, Res, Treatment))
## Fixing column names
## Need to strip out parentheses from column names
names(chemData) <- gsub(".(", ".", names(chemData), fixed = TRUE)
names(chemData) <- gsub(")", "", names(chemData), fixed = TRUE)
## Fix total leaf and total stem phenolics
names(chemData)[names(chemData) == "Total.Phenolics.ppm"] <- "total.phenolics.stem"
names(chemData)[names(chemData) == "total.phenolics.ppm"] <- "total.phenolics.leaf"
## Strip spaces and hyphens out of column names
names(chemData) <- gsub(" ", "", names(chemData), fixed = TRUE)
names(chemData) <- gsub("-", "", names(chemData), fixed = TRUE)
names(chemData)
## Need to convert leading numerics to words within column names 
for(i in 1:length(names(chemData))){
  if(grepl("^[[:digit:]]", names(chemData)[i])){
    splitName <- strsplit(names(chemData)[i], split = "")[[1]] 
    word1 <- number2words(splitName[1]) # Convert leading number to a word
    newName <- paste(c(word1, splitName[-1]), collapse = "") # Replace leading number with the word
    names(chemData)[i] <- newName
  }
}
names(chemData)
## Removed all the leading numbers in column names?
any(grep("^[[:digit:]]", names(chemData)))

#### Save data set of all chemistry data
saveRDS(chemData, "output/full_phenolic_and_volatile_data_pdr1_2017.rds")


#### I'm pretty sure one observation of Res = s, Week = 5, Treatment = Both should be Res = r....
## Try to figure out which one
## Try to look at residuals for all compounds
probObs <- chemData %>% dplyr::filter(Week == 5 & Res == "s" & Treatment == "Both") %>% 
  dplyr::select(10:ncol(chemData))
chemMeans <- colMeans(probObs)
repResiduals <- matrix(ncol = ncol(probObs), nrow = nrow(probObs))
for(i in 1:ncol(repResiduals)){
  repResiduals[,i] <- (probObs[,i] - chemMeans[i])^2
}
sumResiduals <- rowSums(repResiduals)
hist(sumResiduals)
which(sumResiduals == max(sumResiduals))
## Total squared residuals of Rep 9 is much larger than other replicates across all other reps
## Rep 9 might be the incorrectly labeled rep

