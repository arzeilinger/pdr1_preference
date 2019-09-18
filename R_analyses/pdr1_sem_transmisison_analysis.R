#### Analysis of PdR1 transmission data using structural equation modeling

rm(list = ls())

my.packages <- c("tidyr", "dplyr", "data.table", "ggplot2", "piecewiseSEM")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/selectSEM.R")

#### For CMM parameters
## choice1 = source (Xylella infected) plant
## choice2 = test plant

#### Some (controversial) rules of thumb for assessing absolute model fit
## Fisher's C small and p-value from chi-squared > 0.05


##########################################################################################################
#### 2017 SEM analysis
##########################################################################################################

#### Load transmission, acquisition, culturing, preference (and phenolic, when it's ready) data set
transVCPdata <- readRDS("output/complete_2017_transmission-preference_dataset.rds")
str(transVCPdata)
summary(transVCPdata)


#### Select only relevant variables and clean data up
#### Remove mu2 outlier
semData17 <- transVCPdata %>% 
  mutate(acquisition = asin(sqrt(propInfectious)),
         ## Trt code: R = 1, S = 0; so that the variable indicates the presence/absence of PdR1 trait...
         ## and a positive relationship between trt and infection status indicates greater transmission from Resistant
         trtNumeric = 1 - (as.numeric(trt) - 1),
         XylellaPopulation = log(xfpop+1),
         pd_index = log(PD_symptoms_index+1)) %>%
  dplyr::select(test_plant_infection, trtNumeric, week, pd_index, acquisition, XylellaPopulation, mu1, mu2, p1, p2) %>%
  dplyr::filter(complete.cases(.)) %>%
  dplyr::filter(mu2 < 6.5) %>%
  dplyr::rename(transmission = test_plant_infection,
                resistanceTrait = trtNumeric)
str(semData17)

#### Fit 2017 data to full set of models and run model selection
modelSelect17 <- selectSEM(semData17)
(results17 <- modelSelect17[[1]])

## Best model is noPDModel; attraction and leaving rates without pd_index
modelSelect17[[2]]$noPDModel
