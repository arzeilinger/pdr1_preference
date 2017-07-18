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
ps419 <- read.csv(paste(qpcrDir, "04_19_17 PD_0059_plate_setup.csv", sep=""), header = TRUE) 
ps419 <- transformPlateSetupcsv(ps419)

#### Import qpcr data from LinRegPCR
ct419 <- readLinReg(file = "Adam_41917_PD0059_linreg.xlsx", dir = qpcrDir)
# Merge sample names from plate setup with data from LinRegPCR
ct419 <- ps419 %>% dplyr::select(., wellNumber, sample) %>% 
  left_join(., ct419, by = "wellNumber") %>% 
  dplyr::filter(., !is.na(sample))


#### Calculate CFUs from standard curve
cfulist419 <- calculateCFU(qpcrdata = ct419, serial_dilution = serial_dilution, getModel = TRUE)
scurve419 <- cfulist419[[1]]
mod419 <- cfulist419[[2]]
cfu419 <- cfulist419[[3]]
# Checking results and standard curve coefficients

