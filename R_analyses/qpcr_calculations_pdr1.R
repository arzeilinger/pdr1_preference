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
cfus <- lapply(1:length(cfuList), function(x) cfuList[[x]][[3]])

