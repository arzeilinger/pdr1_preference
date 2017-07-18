#### Analysis of mortality in BGSS cages

rm(list = ls())
# libraries
my.packages <- c("openxlsx", "tidyr", "dplyr", "data.table", "googlesheets")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/pdr1_functions.R")

##############################################################
#### Extract colony data directly from google sheets
##############################################################

# Register BGSS colony google sheets
gsbgss <- gs_url("https://docs.google.com/spreadsheets/d/1ar4BDAHprFBn08PG9Es636We7e4_RJwbHIaX9cxDwnw/edit#gid=968476353",
                 visibility = "private")

# Download sheets for each active colony
# Need to change the vector of integers to reflect the index of active sheets
# Making multiple sheet requests requires 8 second wait between requests, to not get shutdown by Google API
# So this takes a while
cages <- lapply(1:11, gsExtract) 

cages <- cages %>% rbindlist() %>% as.data.frame()
names(cages) <- c("date", "adults", "nymphs", "notes", "sheetNum")

# Remove rows with no data
cages <- cages %>% dplyr::filter(., !(is.na(adults) & is.na(nymphs)))


########################################################################################
#### Get current number of Adults and Nymphs
########################################################################################

# Need to update the date of the most current colony counts (i.e., date plants were changed)
currentNumbers <- cages %>% dplyr::filter(., date == "2017-07-17")
colSums(currentNumbers[,c("adults", "nymphs")])


#########################################################################################
#### Assess mortality rate of adults from peak abundance to most recent date
#########################################################################################

# Make into list of cages again
cageList <- cages %>% split(., cages$sheetNum)

mortalityRates <- sapply(cageList, mortality, simplify = TRUE) %>% t() %>% as.data.frame()
names(mortalityRates) <- c("rate", "adults", "nymphs")
mortalityRates

# Average mortality rate
meanMortality <- mean(mortalityRates$rate, na.rm = TRUE)
# Total current number of adults
totalAdults <- sum(mortalityRates$adults)
# Total current number of nymphs
totalNymphs <- sum(mortalityRates$nymphs, na.rm = TRUE)

# How many do I expect to die over the next 9 weeks?
totaldead <- -meanMortality*9
remaining <- totalAdults - totaldead 
# I need more bugs:
512 - remaining

