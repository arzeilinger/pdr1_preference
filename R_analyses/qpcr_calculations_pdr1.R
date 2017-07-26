#### CALCULATING CFU FROM QPCR DATA
#### Ct values are estimated from LinRegQPCR program

rm(list = ls())
# libraries
# loading dtplyr that replaces dplyr and data.table
my.packages <- c("openxlsx", "tidyr", "dplyr", "ggplot2")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/qpcrFunctions.R")

# Directory for qPCR data
qpcrDir <- "data/2016_data/qpcr_data/"

#### Import serial dilution
serial_dilution <- read.xlsx(paste(qpcrDir, "serial_dilution.xlsx", sep = ""), sheet = "dilutions_R")

#### First qpcr run
#### Import plate setup and transform plate setups
ps419 <- read.csv(paste(qpcrDir, "041917_PD0059_plate_setup.csv", sep=""), header = TRUE) 
ps419 <- transformPlateSetupcsv(ps419)

#### Import qpcr data from LinRegPCR
ct419 <- readLinReg(file = "Adam_041917_PD0059_linreg.xlsx", dir = qpcrDir)
# Merge sample names from plate setup with data from LinRegPCR
ct419 <- ps419 %>% dplyr::select(., wellNumber, sample) %>% 
  left_join(., ct419, by = "wellNumber") %>% 
  dplyr::filter(., !is.na(sample))


#### Calculate CFUs from standard curve
cfulist419 <- calculateCFU(qpcrdata = ct419, serial_dilution = serial_dilution, getModel = TRUE)
scurve419 <- cfulist419[[1]]
mod419 <- summary(cfulist419[[2]])
cfu419 <- cfulist419[[3]]
# Checking results and standard curve coefficients


####################################################################################################
#### Processing qPCR data for all plates
####################################################################################################

# Define vector of plate setup files and linregPCR output files
plateNames <- c("041917_PD0059", "042017_PD0059_2", "042017_PD0059", "042117_PD0059", "042417_PD0059", "042517_PD0059_2", "042517_PD0059")
psfiles <- paste(plateNames, "plate_setup.csv", sep = "_")
outputfiles <- paste("Adam", plateNames, "linreg.xlsx", sep = "_")

# Import and transform plate setups
plateSetupList <- lapply(psfiles, function(x) read.csv(paste(qpcrDir, x, sep=""), header = TRUE) %>% transformPlateSetupcsv())

# Import and merge linregPCR output files
outputList <- lapply(1:length(plateNames), function(x){
  readLinReg(file = outputfiles[[x]], dir = qpcrDir) %>%
    left_join(., plateSetupList[[x]][,c("wellNumber", "sample")], by = "wellNumber") %>%
    dplyr::filter(., !is.na(sample))
})

# Calculate CFUs for each plate using the standard curves
cfuList <- lapply(outputList, function(x) calculateCFU(qpcrdata = x, serial_dilution = serial_dilution, getModel = TRUE))
# get the standard curve regression models out
sModels <- lapply(1:length(cfuList), function(x) cfuList[[x]][[2]] %>% summary())
# extract CFU data and make into a single data.frame
vectorcfu <- lapply(1:length(cfuList), function(x) cfuList[[x]][[3]]) %>% rbindlist() %>% as.data.frame() 
saveRDS(vectorcfu, file = "output/pdr1_2016_vector_cfu_from_qpcr.rds")

#### The treatments (genotype-week) are not included in the plate setup, but are associated with a tube code
#### Merge cfu data sets with treatment codes, from a Google Sheet

require(googlesheets)

# Register BGSS colony google sheets
codesgs <- gs_url("https://docs.google.com/spreadsheets/d/1XX5HnGltAvShBskd5of5RNUZfhKrz35-y1WO00QlxpM/edit#gid=0",
                    visibility = "private")
tubeCodes <- gs_read(codesgs, ws = 1, range = "A1:C277")
# Split date and week information
tubeCodes$week <- tstrsplit(tubeCodes$treatment, " ")[[4]]
tubeCodes$week[tubeCodes$week == 6] <- 8  # Change week 6 to week 8
tubeCodes$treatment <- tstrsplit(tubeCodes$insect_code, "")[[1]]
tubeCodes$rep <- tstrsplit(tubeCodes$insect_code, "")[[2]]
tubeCodes$tube_code <- factor(tubeCodes$tube_code)

# Some of the sample names in the CFU data set include both the tube code and cage rep; I need to strip out cage rep
samplesplit <- tstrsplit(vectorcfu$sample, " ")
vectorData <- cbind(vectorcfu, samplesplit[[1]], samplesplit[[2]]) 
names(vectorData) <- c(names(vectorcfu), "tube_code", "sample2") 
# check what the NAs are for, then delete them
vectorData[is.na(vectorData$tube_code),]
vectorData <- vectorData %>% dplyr::filter(!is.na(tube_code))

# merge tubeCodes and vector cfu data
vectorData <- vectorData %>% full_join(., tubeCodes, by = "tube_code")
saveRDS(vectorData, file = "output/pdr1_2016_vector_cfu_from_qpcr.rds")

# check that the genotype-rep codes align between data sets, where they were included in the pcr sample names
nameCheck <- vectorData[!is.na(vectorData$sample2),]
with(nameCheck, sum(sample2 != paste("(", treatment, rep, ")", sep="")))
nameCheck[which(nameCheck$sample2 != paste("(", nameCheck$treatment, nameCheck$rep, ")", sep="")),]
# Need to check these samples

vectorData <- readRDS("output/pdr1_2016_vector_cfu_from_qpcr.rds")
# Look at discrepancies between qPCR duplicates
repCheck <- vectorData %>% group_by(tube_code) %>% summarise(mean = mean(Cq, na.rm = TRUE), 
                                                             min = min(Cq, na.rm = TRUE), 
                                                             max = max(Cq, na.rm = TRUE))
repCheck$diff <- repCheck$max - repCheck$min
# Select only those samples that have medium to large differences between duplicates
badReps <- repCheck[repCheck$diff > 5,]
badReps %>% print(., n = nrow(.))

# Samples to prioritize for re-running qPCR
badRepData <- vectorData %>% dplyr::filter(vectorData$tube_code %in% badReps$tube_code) %>% 
  dplyr::select(tube_code, insect_code, week, Cq, Quality_checks)
write.csv(badRepData, file = "output/qpcr_bad_replicates.csv", row.names = FALSE)
# Remove serial dilution and NTC rows
vectorData <- vectorData %>% dplyr::filter(!is.na(insect_code))

