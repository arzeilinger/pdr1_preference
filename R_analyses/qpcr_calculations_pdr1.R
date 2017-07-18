#### CALCULATING CFU FROM QPCR DATA
#### Ct values are estimated from LinRegQPCR program

rm(list = ls())
# libraries
# loading dtplyr that replaces dplyr and data.table
my.packages <- c("openxlsx", "tidyr", "dplyr", "ggplot2")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/qpcrFunctions.R")

# Directory for qPCR data
qpcrDir <- "data/qpcr_data/"

#### Import serial dilution
serial_dilution <- read.xlsx("data/qpcr_data/serial_dilution.xlsx", sheet = "dilutions_R")

#### First qpcr run
#### Import plate setup and transform plate setups
ps419 <- read.csv(paste(qpcrDir, "04_19_17 PD_0059_plate_setup.csv", sep=""), header = TRUE) 
ps419 <- transformPlateSetupcsv(ps419)

#### Import qpcr data from LinRegPCR
ct419 <- readLinReg(file = "Adam_41917_PD0059_linreg.xlsx", dir = qpcrDir)

