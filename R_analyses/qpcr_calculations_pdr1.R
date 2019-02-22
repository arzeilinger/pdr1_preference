#### CALCULATING CFU FROM QPCR DATA
#### Ct values are estimated from LinRegQPCR program

rm(list = ls())
# libraries
# loading dtplyr that replaces dplyr and data.table
my.packages <- c("openxlsx", "tidyr", "dplyr", "ggplot2", "googlesheets", "modeest")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/qpcrFunctions.R")

# Directory for qPCR data
qpcrDir <- "data/2016_data/Jeffs_qpcr_data/"

#### Import serial dilution
serial_dilution <- read.xlsx(paste(qpcrDir, "serial_dilution.xlsx", sep = ""), sheet = "dilutions_R")

#### If only final CFU data set is needed
#### FINAL DATA SET
vectorData <- readRDS("output/pdr1_2016_vector_cfu_from_qpcr.rds")



##########################################################################################################
#### To estimate CFU for only one plate
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
vectorcfu %>% dplyr::select(sample, cfu) %>% tail()


###############################################################################################################
#### The treatments (genotype-week) are not included in the plate setup, but are associated with a tube code
#### Merge cfu data sets with treatment codes, from a Google Sheet


# Register BGSS colony google sheets
codesgs <- gs_url("https://docs.google.com/spreadsheets/d/1XX5HnGltAvShBskd5of5RNUZfhKrz35-y1WO00QlxpM/edit#gid=0",
                    visibility = "private")
tubeCodes <- gs_read(codesgs, ws = 1, range = "A1:C277")
saveRDS(tubeCodes, file = "data/Jeffs_tube_codes_raw_googlesheet.rds")
# Split date and week information
tubeCodes <- readRDS("data/Jeffs_tube_codes_raw_googlesheet.rds")
tubeCodes$week <- tstrsplit(tubeCodes$treatment, " ")[[4]]
tubeCodes$week[tubeCodes$week == 6] <- 8  # Change week 6 to week 8
tubeCodes$trt <- tstrsplit(tubeCodes$insect_code, "")[[1]]
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

#### FINAL DATA SET
vectorData <- readRDS("output/pdr1_2016_vector_cfu_from_qpcr.rds")



#####################################################################################################
#### Re-running 2016 PdR1 vector DNA extracts on qPCR

qpcrDir <- "data/2016_data/qpcr_data/"

#### Import plate setup 
## From googlesheets (Michael's template)
# tsgs_url <- gs_url("https://docs.google.com/spreadsheets/d/1Xi-S8UeRQWufjaeWKYzJYN590kHX0qOJ2kweCVVg8Kw/edit#gid=0",
#                   visibility = "private")
# plateSetup1 <- gs_read(tsgs_url, ws = 1, range = "A6:M14")  %>% transformPlateSetupcsv()

## From excel file
plateSetup <- read.xlsx(paste(qpcrDir, "2019-02-07_pdr1_qpcr_plate_setup.xlsx", sep = ""), 
                        sheet = 1, rows = c(6:14), cols = c(1:13)) %>% 
  transformPlateSetupcsv()

#### Read in lineRegPCR file
qpcrOutput <- readLinReg(file = "2019-02-07_pdr1.xlsx", dir = qpcrDir) %>%
  left_join(., plateSetup[,c("wellNumber", "sample")], by = "wellNumber") %>%
  dplyr::filter(., !is.na(sample))
qpcrOutput$Cq[qpcrOutput$Cq >= 41] <- 0

# Looking at results
(qpcrCheck <- qpcrOutput[,c("Cq", "N0", "sample", "Sample_Use", "Quality_checks")] %>% arrange(., Cq)) 
# I think Sanjeet switched 2PC (+ control) with NTC (- control) samples

#### Merge these re-run samples with the final 2016 data set
#### FINAL DATA SET
vectorData <- readRDS("output/pdr1_2016_vector_cfu_from_qpcr.rds")

## filter to just 2016 vectors; use the fact that the 2017 samples are the largest numbers
qpcrOutput <- qpcrOutput %>% mutate(sampleNum = as.numeric(sample)) %>% arrange(sampleNum)
rerun16 <- qpcrOutput %>% dplyr::filter(sampleNum < 800 & !is.na(sampleNum))

## Merge rerun samples data with tubeCodes data set
rerun16 <- rerun16 %>% 
  dplyr::select(-sampleNum) %>%
  left_join(., tubeCodes, by = c("sample" = "tube_code")) %>% 
  rename(tube_code = sample) # Switch name to tube_code to match vectorData
## rbind rerun samples with vectorData
vectorData2 <- list(vectorData, rerun16) %>% rbindlist(., fill = TRUE)

#### Clean up data set
## Define vector infection status based on Cq
vectorData2$vectorInfected <- ifelse(vectorData2$Cq > 0 & vectorData2$Cq < 41, 1, 0) 

## Save data set with rerun samples
saveRDS(vectorData2, file = "output/pdr1_2016_vector_cfu_from_qpcr.rds")



#################################################################################################################
## Data quality checkts
# check that the genotype-rep codes align between data sets, where they were included in the pcr sample names
nameCheck <- vectorData[!is.na(vectorData$sample2),]
with(nameCheck, sum(sample2 != paste("(", trt, rep, ")", sep="")))
nameCheck[which(nameCheck$sample2 != paste("(", nameCheck$trt, nameCheck$rep, ")", sep="")),]
# Need to check these samples
# I checked and the insect_code is correct, sample2 is incorrect

# Look at discrepancies between qPCR duplicates
#### Summarize duplicates. 
vectorData[vectorData$Cq >= 41, "Cq"] <- 0
repCheck <- vectorData %>% group_by(tube_code) %>% summarise(minCq = min(Cq, na.rm = TRUE),
                                                             maxCq = max(Cq, na.rm = TRUE),
                                                             maxDiff = maxCq - minCq,
                                                             # round Cq values to nearest 10 then calculate the mode (most common value)
                                                             # Not sure this modeCq is working correctly, need to check it, maybe try other mode functions
                                                             #modeCq = mfv(round(Cq, digits = -1))[1], 
                                                             nSample = length(Cq),
                                                             # Need to carryover experiment IDs
                                                             trt = first(trt),
                                                             treatment = first(treatment),
                                                             rep = first(rep),
                                                             week = first(week))

## Samples that have a large difference among duplicates and need to be run again
badReps <- repCheck %>% dplyr::filter(maxDiff > 5) %>% mutate("tube_code" = as.numeric(tube_code)) %>%
  dplyr::filter(!is.na(tube_code)) 
badReps %>% print.data.frame()
badReps %>% nrow()

write.csv(badReps, file = "output/qpcr_bad_replicates_2016.csv", row.names = FALSE)


#### Select samples to run again, to try to improve efficiency (2018-03-01)
badRepTubes <- unique(badRepData$tube_code) %>% as.numeric
badReps.i <- (which(badRepTubes == max(badRepTubes))-10):which(badRepTubes == max(badRepTubes))
badTubesSelected <- badRepTubes[badReps.i]

tubes <- vectorData %>% dplyr::filter(!is.na(insect_code)) 
tubes <- tubes$tube_code %>% as.numeric %>% unique 
goodTubes <- setdiff(tubes, badRepTubes)

# 21 samples to run qPCR again
selectTubes <- c(sample(goodTubes, 10), badTubesSelected)

# Check the info about these selected samples
selectInfo <- vectorData %>% dplyr::filter(tube_code %in% selectTubes)

# Save info
write.csv(selectInfo, "output/samples_for_qpcr_re-run_2018-03-02.csv", row.names = FALSE)
write.csv(selectTubes, "output/just_tube_codes_for_qpcr_re-run_2018_03_02.csv", row.names = FALSE)


#### Are all consequtive tube codes represented in the data set?
samplecodes <- vectorData %>% dplyr::filter(!grepl("D", vectorData$tube_code)) %>% dplyr::select(tube_code) 
samplecodes <- as.numeric(samplecodes$tube_code) %>% unique() %>% sort()
# Yes, so I can set up re-run plates with consequtive tube-code numbers 



