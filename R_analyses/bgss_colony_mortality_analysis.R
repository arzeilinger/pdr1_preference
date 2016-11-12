#### Analysis of mortality in BGSS cages

rm(list = ls())
# libraries
my.packages <- c("xlsx", "tidyr", "dplyr", "data.table", "googlesheets")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/pdr1_functions.R")

##############################################################
#### Extract colony data directly from google sheets
##############################################################

# Register BGSS colony google sheets
gsbgss <- gs_url("https://docs.google.com/spreadsheets/d/1g6EuSFKA7AM3VzmF80dKix6ED3vBIle3UF7h-RFHuWo/edit#gid=968476353",
                 visibility = "private")

cages <- lapply(1:6, gsExtract) %>% rbindlist() %>% as.data.frame()
names(cages) <- c("date", "adults", "nymphs", "notes", "sheetNum", "weeks")

# Code works but API requests are limited to 6 at a time, so can't download all sheets

##############################################################
#### Extract colony data from .xlsx file
##############################################################

cages <- lapply(1:6, excelExtract) %>% rbindlist() %>% as.data.frame()
names(cages) <- c("date", "adults", "nymphs", "notes", "sheetNum", "weeks")


#####################################################################
##### Clean up data
# Remove rows with no data
cages <- dplyr::filter(cages, !(is.na(adults) & is.na(nymphs)))
# Make adults and nymphs columns numeric
cages$adults <- as.numeric(levels(cages$adults))[cages$adults]
cages$nymphs <- as.numeric(levels(cages$nymphs))[cages$nymphs]
cages$weeks <- as.numeric(cages$weeks)

# Make into list of cages again
cageList <- cages %>% split(., cages$sheetNum)

########################################################################################
#### Get current number of Adults and Nymphs
########################################################################################

currentNumbers <- sapply(cageList, currentFunc, simplify = TRUE) %>% t() %>% as.data.frame()

# Total adult count in colony cages
colSums(currentNumbers, na.rm = TRUE)
# V1 = adult count
# V2 = nymph count


#########################################################################################
#### Assess mortality rate of adults from peak abundance to most recent date
#########################################################################################

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


#### Extra code

# Modifying dates. Doesn't work well, dates are a mess
dates <- strsplit(as.character(cages$date), split = "^")
dates2 <- lapply(1:length(dates), function(x) ifelse(dates[[x]][1] == "^", 
                                                     paste0(dates[[x]][-1], collapse = ""),
                                                     paste0(dates[[x]], collapse = "")))

