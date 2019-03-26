#### CALCULATING CFU FROM QPCR DATA
#### Ct values are estimated from LinRegQPCR program

rm(list = ls())
# libraries
my.packages <- c("openxlsx", "tidyr", "dplyr", "ggplot2", "googlesheets", "data.table", "modeest")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/qpcrFunctions.R")

# Directory for qPCR data
qpcrDir <- "data/2017_data/qpcr_data/"

#### Import serial dilution
#serial_dilution <- read.xlsx(paste(qpcrDir, "serial_dilution.xlsx", sep = ""), sheet = "dilutions_R")

#### If only final CFU data set is needed
#### FINAL DATA SET
#vectorData <- readRDS("output/pdr1_2016_vector_cfu_from_qpcr.rds")





####################################################################################################
#### Processing qPCR data for all plates
####################################################################################################

# Define vector of plate setup files and linregPCR output files
# plateNames <- c("2018-06-14_pdr1_2017", "2018-06-18_pdr1_2017", "2018-07-13_pdr1_2017", "2018-07-17_pdr1_2017", 
#                 "2018-07-23_pdr1_2017", "2018-07-30_pdr1_2017", "2018-08-24_pdr1_2017", "2018-08-31_pdr1_2017",
#                 "2018-09-05_pdr1_2017", "2018-09-07_pdr1_2017")
# psfiles <- paste(plateNames, "qpcr_plate_setup.xlsx", sep = "_")
# outputfiles <- paste(plateNames, ".xlsx", sep = "")

## Define vector of plate setup files and linregPCR output files
## Use list.files to get all the file names, rather than inputing the dates
qpcrFiles <- list.files(path = qpcrDir)
# Plate set up files
psfiles <- qpcrFiles[grep(qpcrFiles, pattern = "qpcr_plate_setup")]
# Files of outputs from LinRegPCR -- the qpcr results
outputfiles <- qpcrFiles[which(!grepl(qpcrFiles, pattern = "plate_setup") & 
                                 !grepl(qpcrFiles, pattern = "extraction") & 
                                 !grepl(qpcrFiles, pattern = "serial_dilution") & 
                                 grepl(qpcrFiles, pattern = ".xlsx"))]

# Import and transform plate setups
plateSetupList <- lapply(psfiles, function(x) read.xlsx(paste(qpcrDir, x, sep=""), sheet = 1, colNames = FALSE) %>% 
                           transformPlateSetupMV() %>%
                           dplyr::select(., wellNumber, well, sample))

# Import and merge linregPCR output files
outputList <- lapply(1:length(outputfiles), function(x){
  readLinReg(file = outputfiles[[x]], dir = qpcrDir) %>%
    left_join(., plateSetupList[[x]], by = "wellNumber") %>%
    dplyr::filter(., !is.na(sample)) %>%
    mutate("plateName" = outputfiles[x])
})

## Fix specific sample names in 2018-06-14 and 2018-08-31 plates
## These sample names should be numeric but are character because of a space in front. Not sure why.
outputList[[1]][outputList[[1]]$sample == " 9", "sample"] <- 9
outputList[[8]][outputList[[8]]$sample == " 92", "sample"] <- 92
outputList[[8]][outputList[[8]]$sample == " 96", "sample"] <- 96

#### Filter out 2016 samples from the last plate: 2019-02-07
## The 2017 samples have tube code numbers > 800; the 2016 numbers are < 300
outputList[[29]] <- outputList[[29]] %>% mutate(sample = as.numeric(sample)) %>% dplyr::filter(sample > 800)

#### Import codes for experimental treatments, and merge with tube codes
tubeCodes <- read.xlsx(paste(qpcrDir, "extraction_tube_codes_pdr1_2017.xlsx", sep = ""))
str(tubeCodes)

tubeCodes %>% dplyr::filter(week == 2 & block == 1 & genotype == "102" & rep == 4)


## Join tube codes and plate runs
qpcrData <- outputList %>%
  rbindlist() %>%
  mutate(sampleNum = as.numeric(sample)) %>%
  left_join(., tubeCodes, by = c("sampleNum" = "tube_code"))
## Set any Cq values above 40 to 0
qpcrData$Cq[qpcrData$Cq >= 41] <- 0


#### Summarize duplicates. 
vectorData <- qpcrData %>% group_by(sample) %>% summarise(minCq = min(Cq, na.rm = TRUE),
                                                          maxCq = max(Cq, na.rm = TRUE),
                                                          maxDiff = maxCq - minCq,
                                                          # round Cq values to nearest 10 then calculate the mode (most common value)
                                                          # Not sure this modeCq is working correctly, need to check it, maybe try other mode functions
                                                          modeCq = mfv(round(Cq, digits = -1))[1], 
                                                          nSample = length(Cq),
                                                          # Need to carryover experiment IDs
                                                          genotype = first(genotype),
                                                          trt = first(trt),
                                                          rep = first(rep),
                                                          week = first(week),
                                                          block = first(block),
                                                          plateName = last(plateName))

## Samples that have a large difference among duplicates and need to be run again
vectorData %>% dplyr::filter(maxDiff > 5) %>% mutate("sample" = as.numeric(sample)) %>%
  dplyr::filter(!is.na(sample)) %>% arrange(plateName) %>% tail(n = 20) #%>% print.data.frame()

qpcrData %>% dplyr::filter(sample == 901)

qpcrData %>% dplyr::filter(plateName == "2018-11-28_pdr1_2017.xlsx") %>% dplyr::select(well, sample, Cq, N0, Sample_Use, Quality_checks) %>% arrange(sample)

check <- qpcrData %>% dplyr::filter(week == 14 & block == 1 & genotype == "102" & rep == 4)
check


#### Editing results of sample 901 "by hand"
## Because of how I calculate the Cq mode for each sample, sample 901 gets a "0" for most frequent Cq, which is not correct
## Also, sample 901 was the only recovered vector from cage 14-1-102R4, and the test plant was positive
## So sample 901 should be "positive"
vectorData[vectorData$sample == 901, "modeCq"] <- 40

# Make a binary "infectious" variable for each bug
# For now just take the max Cq value and if 0 < Cq < 40 then call the bug "infected"
# Then calculate proportion of vectors infectious in each cage
vectorData2 <- vectorData %>% mutate(infectious = ifelse(modeCq > 0, 1, 0)) %>%
  # Calculate proportion of bugs infectious for each cage
  group_by(genotype, trt, rep, week, block) %>% summarise(nbugs = length(infectious),
                                                          totalInfectious = sum(infectious),
                                                          propInfectious = totalInfectious/nbugs)

## Remove qPCR control and PA (Preference Assumption test) samples
vectorData3 <- vectorData2 %>% dplyr::filter(!is.na(genotype) & !grepl("PA", genotype))
# Fix column classes
vectorData3$block <- factor(vectorData3$block)
vectorData3$genotype <- factor(vectorData3$genotype)
vectorData3$trt <- factor(vectorData3$trt)

#### Save acquisition data
saveRDS(vectorData3, file = "output/vector_acquisition_data_pdr1_2017.rds")

#### Importing and munging transmission data
# Importing from local .xlsx file
# transdata <- read.xlsx("data/2017_data/PdR1_2017_preference-transmission_experiment_data.xlsx", sheet = "test_plant_culturing", detectDates = TRUE)
# Import from Googlesheets
pdr1DataURL <- gs_url("https://docs.google.com/spreadsheets/d/14uJLfRL6mPrdf4qABeGeip5ZkryXmMKkan3mJHeK13k/edit?usp=sharing",
                      visibility = "private")
transdataGS <- gs_read(pdr1DataURL, ws = "test_plant_culturing")

transdata <- transdataGS

# Remove Control plant samples
transdata <- transdata %>% dplyr::filter(!grepl("CTRL", genotype))

# Combine test_plant_infection_1 and test_plant_infection_2 columns
# In all cases, when I re-cultured a sample (test_plant_infection_2), the results were the same as the 1st time or more reliable
# So go with test_plant_infection_2 results when they are available
transdata$test_plant_infection <- with(transdata, ifelse(!is.na(test_plant_infection_2), test_plant_infection_2, test_plant_infection_1)) %>% as.integer()                                      

# Fix column classes
transdata$genotype <- factor(transdata$genotype)
transdata$trt <- factor(transdata$trt)
transdata$block <- factor(transdata$block)

str(transdata)
summary(transdata)




#### Plot time series of vector infections
vectorPlot <- ggplot(data=transVectorSummary, aes(x=week, y=meanPropInfectious)) +
  geom_line(aes(linetype=genotype, colour = trt), size=1.25) +
  geom_point(aes(shape=genotype, colour = trt), size=3.5) +
  geom_errorbar(aes(ymax=meanPropInfectious+sePropInfectious, ymin=meanPropInfectious-sePropInfectious), width=0.2) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = c(2,5,8,14), limits = c(1,14)) + 
  scale_y_continuous(name = "Mean proportion of infectious vectors",
                     limits = c(0,1)) +
  # ylab("% insects on source plant") + 
  # ylim(c(0,100)) +
  # xlab("Weeks post inoculation") +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

vectorPlot

ggsave("results/figures/2017_figures/infectious_vectors_line_plot_2017.jpg", plot = vectorPlot,
       width = 7, height = 7, units = "in")


###############################################################################################################
###############################################################################################################
#### Analysis of acquisition data

#### Final aquisition data, at the cage level
vectorData3 <- readRDS("output/vector_acquisition_data_pdr1_2017.rds")
str(vectorData3)
summary(vectorData3)

#### Poisson regression with experimental treatments
acqMod <- glm(totalInfectious ~ block + week*genotype, data = vectorData3, family = "poisson")
plot(simulateResiduals(acqMod))
testOverdispersion(simulateResiduals(acqMod))
## Data are definitely overdispered
summary(acqMod)


## quasi-Possion regression
acqMod2 <- glm(totalInfectious ~ block + week*genotype, data = vectorData3, family = "quasipoisson")
plot(acqMod2)
summary(acqMod2)


##############################################################################################################
#### Fitting non-linear models to transmission data

# Load NLL functions for all models
source("R_functions/nll_functions_pdr1.R")


#### Initial parameters
## Initial parameters for linear model
a.l0 <- 0 # Y-intercept assumed at 0
b.l0 <- 0.5/14 # Looks like the line might hit 50% at 14 weeks, calculate initial slope from this

## Initial parameters for Holling Type IV
a0 <- 0.1 # Asymptotic propInfected
# Peak transmission is at 8 weeks, peak = -2b/c, b/c = -4, and c should be < 0
b0 <- 8
c0 <- -2

## Initial parameters for Ricker model
a.rick0 <- 0.25/5 # Initial slope of the line
b.rick0 <- 1/8 # x value for the peak

## Initial parameters for Logistic Growth model
b.logist0 <- 0.25/5/4 # Initial slope of the line divided by 4
a.logist0 <- -6*b.logist0 # Negative value of x at halfway to asymptote (inflection point) multiplied by b

## Initial parameters for Michaelis-Menten model
a.mm0 <- 0.5 # Asymptote
b.mm0 <- 6 # value of x when y = a/2


#### Fitting multiple non-linear models to Resistant and Susceptible trts separately

# Get data sets
vectorR <- vectorData3[vectorData3$trt == "R",]
vectorS <- vectorData3[vectorData3$trt == "S",]

#### Fitting resistant line data
vectorRresults <- optimizeVectorModels(vectorR)
(modelSelectR <- vectorRresults[[2]])
opListR <- vectorRresults[[1]]

# Get model predictions for plotting
bestModR <- opListR$holling4Op # Holling4 model is the best
bestcoefR <- as.list(coef(bestModR))

newDatR <- data.frame(week = seq(0, 14, length = 101))
newDatR$totalInfectious <- with(newDatR, (bestcoefR$a*week^2)/(bestcoefR$b + bestcoefR$c*week + week^2))

#### Fitting susceptible line data
vectorSresults <- optimizeVectorModels(vectorS)
(modelSelectS <- vectorSresults[[2]])
opListS <- vectorSresults[[1]]

# Get model predictions for plotting
bestModS <- opListS$holling4Op # Ricker model is the best
bestcoefS <- as.list(coef(bestModS))

newDatS <- data.frame(week = seq(0, 14, length = 101))
newDatS$totalInfectious <- with(newDatS, (bestcoefS$a*week^2)/(bestcoefS$b + bestcoefS$c*week + week^2))


#### Plotting results
## Calculate mean infectiousness for each trt and week
vectorSummary2 <- vectorData3 %>% dplyr::filter(!is.na(totalInfectious)) %>% 
  group_by(week, trt) %>% summarise(meanTotalInfectious = mean(totalInfectious, na.rm = TRUE),
                                    nTotalInfectious = length(totalInfectious),
                                    seTotalInfectious = sd(totalInfectious, na.rm = TRUE)/sqrt(nTotalInfectious)) 
vectorSummary2

# Plot S and R genotypes summarised together
# Use Holling 4 model for S and for R
# Vector acquisition plot
vectorplotNL <- ggplot(data=vectorSummary, aes(x=week, y=meanTotalInfectious)) +
  #geom_line(aes(linetype=genotype, colour = trt), size=1.25) +
  geom_point(aes(shape = trt, colour = trt), size=3.5) +
  geom_errorbar(aes(ymax=meanTotalInfectious+seTotalInfectious, ymin=meanTotalInfectious-seTotalInfectious), width=0.2) +
  # Defining the colors for the lines based on the default ggplot2 colors and ggplot_build()$data
  geom_smooth(data = newDatR, aes(x=week, y=totalInfectious), method = "loess", colour = "#F8766D", se = FALSE) +
  geom_smooth(data = newDatS, aes(x=week, y=totalInfectious), method = "loess", colour = "#00BFC4", se = FALSE) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = c(2,5,8,14)) + 
  scale_y_continuous(name = "Mean number of infectious vectors",
                     limits = c(-0.2,8.2)) +
  # ylab("% insects on source plant") + 
  # ylim(c(0,100)) +
  # xlab("Weeks post inoculation") +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

vectorplotNL

ggsave("results/figures/2017_figures/acquisition_non-linear_plot_2017.jpg", plot = vectorplotNL,
       width = 7, height = 7, units = "in")


###############################################################################################################
###############################################################################################################
#### Look at 4th attempt at standard curve (2018-08-31 qPCR plate)

# CFU/uL counts
serial_dilution <- read.xlsx(paste(qpcrDir, "serial_dilution_2017_4.xlsx", sep=""), sheet = "serial_dilution")

# plate that contains standard curve Cq
plate7 <- qpcrData %>% dplyr::filter(plateName == "2018-08-31_pdr1_2017") %>% 
  # Split standard curve sample names
  # Column "dilution" refers to D0 - D9 serial dilutions or non-SC samples; column "sctrt" refers to whether the dilutions were treated with PMA ("PMA") or not ("NT")
  mutate(dilution = tstrsplit(sample, split = "-")[[1]], sctrt = tstrsplit(sample, split = "-")[[2]]) %>%
  # Merge with CFU counts
  left_join(., serial_dilution, by = "dilution")

# Non-treated standard curve (not treated with PMA)
ntsc <- plate7 %>% dplyr::filter(sctrt == "NT") %>% arrange(dilution)
qplot(x = log10(cfu_per_ul), y = Cq, data = ntsc, ylim = c(0,40))
# Fit a regression line
ntscMod <- lm(Cq ~ log10(cfu_per_ul), data = ntsc[ntsc$cfu_per_ul > 0,])
summary(ntscMod)

# Standard curve treated with PMA
pmasc <- plate7 %>% dplyr::filter(sctrt == "PMA") %>% arrange(dilution)
qplot(x = sample, y = Cq, data = pmasc, ylim = c(0,40))

# Unknowns from plate 7
unknowns7 <- plate7 %>% dplyr::filter(is.na(sctrt)) %>% arrange(Cq)
hist(unknowns7$Cq)




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


#################################################################################################################
## Data quality checkts
# check that the genotype-rep codes align between data sets, where they were included in the pcr sample names
nameCheck <- vectorData[!is.na(vectorData$sample2),]
with(nameCheck, sum(sample2 != paste("(", trt, rep, ")", sep="")))
nameCheck[which(nameCheck$sample2 != paste("(", nameCheck$trt, nameCheck$rep, ")", sep="")),]
# Need to check these samples
# I checked and the insect_code is correct, sample2 is incorrect

# Look at discrepancies between qPCR duplicates
repCheck <- vectorData %>% group_by(tube_code) %>% summarise(mean = mean(Cq, na.rm = TRUE), 
                                                             min = min(Cq, na.rm = TRUE), 
                                                             max = max(Cq, na.rm = TRUE))
repCheck$diff <- repCheck$max - repCheck$min
hist(repCheck$diff)
# Select only those samples that have medium to large differences between duplicates
badReps <- repCheck[repCheck$diff > 5,]
badReps %>% print(., n = nrow(.))
# Samples to prioritize for re-running qPCR
badRepData <- vectorData %>% dplyr::filter(vectorData$tube_code %in% badReps$tube_code) %>% 
  dplyr::select(tube_code, insect_code, week, Cq, Quality_checks) %>%
  dplyr::filter(!is.na(insect_code))
badRepData %>% print.data.frame
length(unique(badRepData$tube_code))

write.csv(badRepData, file = "output/qpcr_bad_replicates.csv", row.names = FALSE)


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



#####################################################################################################
#### Re-running 2016 PdR1 vector DNA extracts on qPCR

qpcrDir <- "data/2016_data/qpcr_data/"

#### Import plate setup 
## From googlesheets (Michael's template)
# tsgs_url <- gs_url("https://docs.google.com/spreadsheets/d/1Xi-S8UeRQWufjaeWKYzJYN590kHX0qOJ2kweCVVg8Kw/edit#gid=0",
#                   visibility = "private")
# plateSetup1 <- gs_read(tsgs_url, ws = 1, range = "A6:M14")  %>% transformPlateSetupcsv()

## From excel file
plateSetup <- read.xlsx(paste(qpcrDir, "2018-04-23_pdr1_2016_rerun_plate_setup.xlsx", sep = ""), 
                        sheet = 1, rows = c(6:14), cols = c(1:13)) %>% 
  transformPlateSetupcsv()

#### Read in lineRegPCR file
qpcrOutput <- readLinReg(file = "2018-04-23_pdr1_2016_rerun.xlsx", dir = qpcrDir) %>%
  left_join(., plateSetup[,c("wellNumber", "sample")], by = "wellNumber") #%>%
  #dplyr::filter(., !is.na(sample))

# Looking at results
(qpcrCheck <- qpcrOutput[,c("Cq", "N0", "sample", "Sample_Use", "Quality_checks")] %>% arrange(., Cq)) 


cfures1 <- calculateCFU(qpcrOutput, serial_dilution = serial_dilution, getModel = TRUE)
cfures1  



# ##########################################################################################################
# #### To estimate CFU for only one plate
# #### First qpcr run (2018-06-14)
# #### Import plate setup and transform plate setups
# ps1 <- read.xlsx(paste(qpcrDir, "2018-06-14_pdr1_2017_qpcr_plate_setup.xlsx", sep=""), sheet = 1, colNames = FALSE) 
# ps1 <- transformPlateSetupMV(ps1)
# 
# #### Import qpcr data from LinRegPCR
# ct1 <- readLinReg(file = "2018-06-14_pdr1_2017.xlsx", dir = qpcrDir)
# # Merge sample names from plate setup with data from LinRegPCR
# ct1 <- ps1 %>% dplyr::select(., wellNumber, well, sample) %>% 
#   left_join(., ct1, by = "wellNumber") %>% 
#   dplyr::filter(., !is.na(sample))
# 
# 
# ## Extract and plot standard curve samples
# scdata <- ct1 %>% dplyr::filter(., grepl("D", sample))
# scdata <- scdata %>% mutate(., dilution = tstrsplit(sample, "-", keep = 1, type.convert = TRUE)[[1]]) %>% arrange(dilution)
# plot(x = as.numeric(factor(scdata$dilution)), y = scdata$Cq)
# 
# 
# ## Quality checks
# # Look at discrepancies between qPCR duplicates
# repCheck <- ct1 %>% group_by(sample) %>% summarise(mean = mean(Cq, na.rm = TRUE), 
#                                                    min = min(Cq, na.rm = TRUE), 
#                                                    max = max(Cq, na.rm = TRUE))
# repCheck$diff <- repCheck$max - repCheck$min
# hist(repCheck$diff)
# repCheck <- repCheck %>% arrange(diff)
# print(repCheck, n = nrow(repCheck))
# # For samples # 17, 18, 19, and 10 the duplicates are way off -- one amplified and one did not. Re-run these. 
# 
# # Histogram of Cq values from experimental samples
# unknowns <- ct1 %>% dplyr::filter(., !grepl("D", sample) & !grepl("PC", sample) & !grepl("NC", sample))
# hist(unknowns$Cq[unknowns$Cq > 0])
# 
# # #### Calculate CFUs from standard curve
# # cfulist419 <- calculateCFU(qpcrdata = ct419, serial_dilution = serial_dilution, getModel = TRUE)
# # scurve419 <- cfulist419[[1]]
# # mod419 <- summary(cfulist419[[2]])
# # cfu419 <- cfulist419[[3]]
# # # Checking results and standard curve coefficients
# 
# 
# 
# ##########################################################################################################
# ##### Second pdr1 2017 plate (2018-06-18)
# #### Import plate setup and transform plate setups
# ps2 <- read.xlsx(paste(qpcrDir, "2018-06-18_pdr1_2017_qpcr_plate_setup.xlsx", sep=""), sheet = 1, colNames = FALSE) 
# ps2 <- transformPlateSetupMV(ps2)
# 
# #### Import qpcr data from LinRegPCR
# ct2 <- readLinReg(file = "2018-06-18_pdr1_2017.xlsx", dir = qpcrDir)
# # Merge sample names from plate setup with data from LinRegPCR
# ct2 <- ps2 %>% dplyr::select(., wellNumber, well, sample) %>% 
#   left_join(., ct2, by = "wellNumber") %>% 
#   dplyr::filter(., !is.na(sample))
# 
# 
# ## Extract and plot standard curve samples
# scdata <- ct2 %>% dplyr::filter(., grepl("D", sample))
# scdata <- scdata %>% mutate(., dilution = tstrsplit(sample, "-", keep = 1, type.convert = TRUE)[[1]]) %>% arrange(dilution)
# plot(x = as.numeric(factor(scdata$dilution)), y = scdata$Cq,
#      xlab = "Serial Dilution", ylab = "Cq")
# 
# ## Import serial dilution and plot with Ct
# serial_dilution <- read.xlsx(paste(qpcrDir, "serial_dilution_2017_2.xlsx", sep = ""), sheet = "serial_dilution")
# # Join ct and cfu data
# scdata <- scdata %>% left_join(., serial_dilution, by = c("sample" = "dilution"))
# # plot cfu/ul and N0
# with(scdata, plot(x = log10(cfu_ul+1), y = log10(N0), 
#                   xlab = "CFU/uL (log10)", ylab = "N0 (log10"))
# with(scdata, abline(lm(log10(N0) ~ log10(cfu_ul+1))))
# # Look at the regression stats
# with(scdata, lm(log10(N0) ~ log10(cfu_ul+1))) %>% summary()
# 
# pos <- ct2 %>% dplyr::filter(Cq > 0 & Cq < 40)
# 
# 
# ## Quality checks
# # Look at discrepancies between qPCR duplicates
# repCheck <- ct2 %>% group_by(sample) %>% summarise(mean = mean(Cq, na.rm = TRUE), 
#                                                    min = min(Cq, na.rm = TRUE), 
#                                                    max = max(Cq, na.rm = TRUE))
# repCheck$diff <- repCheck$max - repCheck$min
# hist(repCheck$diff)
# repCheck <- repCheck %>% arrange(diff)
# print(repCheck, n = nrow(repCheck))
# # For samples # 18, 26, 29, 31, 32, 34, 40, 41, 45, 48 the duplicates are way off -- one amplified and one did not. Re-run these.
# # Also, 5NC, 4NC, and NTC have some contamination (in one duplicate).
# 
# # LinRegPCR's flags:
# flags <- ct2 %>% dplyr::filter(Cq > 0)
# flags
# # Most of the samples that where the duplicates were off also had LinReg flags -- amplification did not plateau -- suggesting that they are spurious positives
# 
# 
# # Histogram of Cq values from experimental samples
# unknowns <- ct2 %>% dplyr::filter(., !grepl("D", sample) & !grepl("PC", sample) & !grepl("NC", sample))
# hist(unknowns$Cq[unknowns$Cq > 0])
# hist(unknowns$Cq, breaks = seq(0,40,by=2))
# 
# 
# 
# 
# ##########################################################################################################
# ##### 3rd pdr1 2017 plate (2018-07-13)
# #### Import plate setup and transform plate setups
# ps3 <- read.xlsx(paste(qpcrDir, "2018-07-13_pdr1_2017_qpcr_plate_setup.xlsx", sep=""), sheet = 1, colNames = FALSE) %>% transformPlateSetupMV()
# 
# #### Import qpcr data from LinRegPCR
# ct3 <- readLinReg(file = "2018-07-13_pdr1_2017.xlsx", dir = qpcrDir)
# # Merge sample names from plate setup with data from LinRegPCR
# ct3 <- ps3 %>% dplyr::select(., wellNumber, well, sample) %>% 
#   left_join(., ct3, by = "wellNumber") %>% 
#   dplyr::filter(., !is.na(sample))
# 
# 
# ## Extract and plot standard curve samples
# scdata <- ct3 %>% dplyr::filter(., grepl("D", sample))
# scdata <- scdata %>% mutate(., dilution = tstrsplit(sample, "-", keep = 1, type.convert = TRUE)[[1]]) %>% arrange(dilution)
# plot(x = as.numeric(factor(scdata$dilution)), y = scdata$Cq,
#      xlab = "Serial Dilution", ylab = "Cq")
# 
# ## Import serial dilution and plot with Ct
# serial_dilution <- read.xlsx(paste(qpcrDir, "serial_dilution_2017_3.xlsx", sep = ""), sheet = "serial_dilution")
# # Join ct and cfu data
# scdata <- scdata %>% left_join(., serial_dilution, by = c("sample" = "dilution"))
# # plot cfu/ul and N0
# with(scdata, plot(x = log10(cfu_ul+1), y = log10(N0), 
#                   #ylim = c(0,40), 
#                   xlab = "CFU/uL (log10)", ylab = "N0 (log10)"))
# with(scdata, abline(lm(log10(N0) ~ log10(cfu_ul+1))))
# # Look at the regression stats
# with(scdata, lm(log10(N0) ~ log10(cfu_ul+1))) %>% summary()
# 
# pos <- ct3 %>% dplyr::filter(Cq > 0 & Cq < 40)
# 
# ## Quality checks
# # Look at discrepancies between qPCR duplicates
# repCheck <- ct3 %>% group_by(sample) %>% summarise(mean = mean(Cq, na.rm = TRUE), 
#                                                    min = min(Cq, na.rm = TRUE), 
#                                                    max = max(Cq, na.rm = TRUE))
# repCheck$diff <- repCheck$max - repCheck$min
# hist(repCheck$diff)
# repCheck <- repCheck %>% arrange(diff)
# print(repCheck, n = nrow(repCheck))
# # Duplicates look very clean
# 
# # LinRegPCR's flags:
# flags <- ct3 %>% dplyr::filter(Cq > 0)
# flags
# 
# 
# # Histogram of Cq values from experimental samples
# unknowns <- ct3 %>% dplyr::filter(., !grepl("D", sample) & !grepl("PC", sample) & !grepl("NC", sample))
# hist(unknowns$Cq[unknowns$Cq > 0])
# hist(unknowns$Cq, breaks = seq(0,40,by=2))
# 
# 
# 
# 
# ####################################################################################################
# #### Compiling plates run so far and combining with treatment data
# 
# ## Combine plate runs
# ctData <- rbind(ct1, ct2, ct3)
# 
# ## Import tube codes
# tubeCodes <- read.xlsx(paste(qpcrDir, "extraction_tube_codes_pdr1_2017.xlsx", sep = ""))
# str(tubeCodes)
# 
# ## Join tube codes and plate runs
# ctData <- ctData %>% 
#   dplyr::filter(., !grepl("D", sample) & !grepl("PC", sample) & !grepl("NC", sample) & !grepl("NTC", sample)) %>% 
#   mutate(sample = as.numeric(sample)) %>%
#   left_join(., tubeCodes, by = c("sample" = "tube_code"))
# 
# ## Combine duplicates. For now just take the max Cq value and if 0 < Cq < 40 then call the bug "infected"
# ## Then calculate proportion of vectors infectious in each cage
# vectorData <- ctData %>% group_by(sample) %>% summarise(maxCq = max(Cq, na.rm = TRUE),
#                                                         # Need to carryover experiment IDs
#                                                         genotype = first(genotype),
#                                                         trt = first(trt),
#                                                         rep = first(rep),
#                                                         week = first(week),
#                                                         block = first(block)) %>%
#   # Make a binary "infectious" variable for each bug
#   mutate(infectious = ifelse(maxCq > 0 & maxCq < 40, 1, 0)) %>%
#   # Calculate proportion of bugs infectious for each cage
#   group_by(genotype, trt, rep, week, block) %>% summarise(n = length(infectious),
#                                                           totalInfectious = sum(infectious),
#                                                           propInfectious = totalInfectious/n)
# 
# 
# #### Importing and munging transmission data
# # Importing from local .xlsx file
# # transdata <- read.xlsx("data/2017_data/PdR1_2017_preference-transmission_experiment_data.xlsx", sheet = "test_plant_culturing", detectDates = TRUE)
# # Import from Googlesheets
# pdr1DataURL <- gs_url("https://docs.google.com/spreadsheets/d/14uJLfRL6mPrdf4qABeGeip5ZkryXmMKkan3mJHeK13k/edit?usp=sharing",
#                       visibility = "private")
# transdataGS <- gs_read(pdr1DataURL, ws = "test_plant_culturing")
# transdata <- transdataGS
# 
# # Remove Control plant samples
# transdata <- transdata %>% dplyr::filter(!grepl("CTRL", genotype))
# 
# # Combine test_plant_infection_1 and test_plant_infection_2 columns
# # In all cases, when I re-cultured a sample (test_plant_infection_2), the results were the same as the 1st time or more reliable
# # So go with test_plant_infection_2 results when they are available
# transdata$test_plant_infection <- with(transdata, ifelse(!is.na(test_plant_infection_2), test_plant_infection_2, test_plant_infection_1)) %>% as.integer()                                      
# 
# # Fix column classes
# transdata$genotype <- factor(transdata$genotype)
# transdata$trt <- factor(transdata$trt)
# transdata$block <- factor(transdata$block)
# 
# str(transdata)
# summary(transdata)
# 
# 
# #### Join qPCR data and transmission data
# vectorData$block <- factor(vectorData$block)
# vectorData$genotype <- factor(vectorData$genotype)
# vectorData$trt <- factor(vectorData$trt)
# 
# transVectorData <- transdata %>% dplyr::select(week, block, genotype, trt, rep, test_plant_infection) %>%
#   left_join(., vectorData, by = c("genotype", "trt", "rep", "week", "block"))
# 
# ## Comparing transmission and acquisition data
# transVectorData %>% dplyr::filter(!is.na(propInfectious)) %>% print.data.frame

