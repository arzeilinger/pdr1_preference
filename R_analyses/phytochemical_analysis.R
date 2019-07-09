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

#### Update chemical compound list
## Import sheet with updated names
newChemNames <- read.xlsx("data/2017_data/BGSS Xylella Compound Table Final.xlsx")
## Replace spaces in names with periods
newChemNames$Old.Name <- gsub(" ", ".", newChemNames$Old.Name)
newChemNames$Final.ID <- gsub(" ", ".", newChemNames$Final.ID)
str(newChemNames)

## Update leaf compounds
newChemNamesLeaf <- newChemNames %>% dplyr::filter(Tissue == "leaf") %>% dplyr::select(Old.Name, Final.ID)
namesToChangeLeaf <- names(leafPhenData)[names(leafPhenData) %in% newChemNamesLeaf$Old.Name]
oldNamesPhenLeaf <- names(leafPhenData)
for(i in 1:length(oldNamesPhenLeaf)){
  if(oldNamesPhenLeaf[i] %in% newChemNamesLeaf$Old.Name){
    oldName.i <- oldNamesPhenLeaf[i]
    newName.i <- newChemNamesLeaf[newChemNamesLeaf$Old.Name == oldName.i, "Final.ID"]
    oldNamesPhenLeaf[i] <- newName.i
  }
}
cbind(oldNamesPhenLeaf, names(leafPhenData))
names(leafPhenData) <- oldNamesPhenLeaf


## Update stem compounds
newChemNamesStem <- newChemNames %>% dplyr::filter(Tissue == "stem") %>% dplyr::select(Old.Name, Final.ID)
namesToChangeStem <- names(stemPhenData)[names(stemPhenData) %in% newChemNamesStem$Old.Name]
oldNamesPhenStem <- names(stemPhenData)
for(i in 1:length(oldNamesPhenStem)){
  if(oldNamesPhenStem[i] %in% newChemNamesStem$Old.Name){
    oldName.i <- oldNamesPhenStem[i]
    newName.i <- newChemNamesStem[newChemNamesStem$Old.Name == oldName.i, "Final.ID"]
    oldNamesPhenStem[i] <- newName.i
  }
}
cbind(oldNamesPhenStem, names(stemPhenData))
names(stemPhenData) <- oldNamesPhenStem
str(stemPhenData)

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



# #### Load foliar and stem volatile data
# leafvolData <- read.xlsx("data/2017_data/Xylella_BGSS_Study_Volatiles_2019-05-09.xlsx", sheet = "Foliar Volatiles")
# stemvolData <- read.xlsx("data/2017_data/Xylella_BGSS_Study_Volatiles_2019-05-09.xlsx", sheet = "Wood Volatiles")
# str(leafvolData)  
# ## Chris said that some observations are missing for stem volatiles because there wasn't enough tissue
# summary(stemvolData) ## Doesn't show NAs, observations were probably not entered into spreadsheet because all NA
# with(stemvolData, table(Week, Res, Treatment))
# with(leafvolData, table(Week, Res, Treatment))
# ## Compare column names
# cbind(names(phenData)[1:10], names(leafvolData)[1:10], names(stemvolData)[1:10])
# ## Looks like the important columns are all the same
# 
# 
# #### Merge phenolic and volatile data
# ## Don't include stem volatile data because of too many missing values
# chemData <- full_join(phenData, leafvolData, by = c("Number", "Date", "Week", "Treatment", "Res", "Rep"))
# with(chemData, table(Week, Res, Treatment))


#### Don't include volatile data
chemData <- phenData
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
saveRDS(chemData, "output/full_phenolic_data_pdr1_2017.rds")

chemData <- readRDS("output/full_phenolic_data_pdr1_2017.rds")
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

#### Check distribution of residuals from R plants
compObs <- chemData %>% dplyr::filter(Week == 5 & Res == "r" & Treatment == "Both") %>% 
  dplyr::select(10:ncol(chemData))
compMeans <- colMeans(compObs)
compResiduals <- matrix(ncol = ncol(compObs), nrow = nrow(compObs))
for(i in 1:ncol(compResiduals)){
  compResiduals[,i] <- (compObs[,i] - compMeans[i])^2
}
sumCompResiduals <- rowSums(compResiduals, na.rm = TRUE)
hist(sumCompResiduals)
which(sumCompResiduals == max(sumCompResiduals))


#### Examining which other replicates are missing
missingReps <- chemData %>% dplyr::filter(Treatment == "Both")
with(missingReps, table(Rep, Week, Res))
## I think Rep 9 Susceptible should be Rep 5 Resistant


# #### Rename Rep 9 Susceptible to Rep 5 Resistant
for(i in 1:nrow(chemData)){
  if(chemData[i, "Rep"] == 9){
    chemData[i, "Rep"] <- 5
    chemData[i, "Res"] <- "r"
  }
}


#### Remove Rep 9
chemData <- chemData %>% dplyr::filter(Rep != 9)
## Check missing reps again
missingReps <- chemData %>% dplyr::filter(Treatment == "Both")
with(missingReps, table(Rep, Week, Res))

#### Save chemistry data
saveRDS(chemData, "output/full_phenolic_data_pdr1_2017.rds")


#### How many pre-trial samples do I have?
chemData <- readRDS("output/full_phenolic_data_pdr1_2017.rds")
chemData %>% dplyr::filter(Treatment == "BothPre") %>% with(., table(Week, Res))
## Not many

## How does pre- and post-trial compare?
prePostTotalPhenolics <- chemData %>% dplyr::filter(Treatment == "BothPre" | Treatment == "Both") %>%
  dplyr::select(Week, Res, Rep, Treatment, total.phenolics.leaf) %>%
  spread(key = Treatment, value = total.phenolics.leaf)
