#### Analysis of PdR1 transmission data using structural equation modeling

rm(list = ls())

my.packages <- c("tidyr", "dplyr", "data.table", "ggplot2", "lavaan", "semPlot", "piecewiseSEM")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/selectSEM.R")

## For tutorial on lavaan package: http://lavaan.ugent.be/tutorial/index.html
## For tutorial on semPlot::semPaths: https://www.r-bloggers.com/ploting-sems-in-r-using-semplot/

#### For CMM parameters
## choice1 = source (Xylella infected) plant
## choice2 = test plant

#### Some (controversial) rules of thumb for assessing absolute model fit
#### For lavaan:
## CFI/TLI>0.95
## RMSEA<0.05
## SRMR<0.06
#### For piecewiseSEM:
## Fisher's C small and p-value from chi-squared > 0.05

##########################################################################################################
##########################################################################################################
#### SEM analysis with lavaan package
##########################################################################################################

#### Specify transmission model
## Best model based on biology/theory
transModel <- '# Regressions
                  transmission ~ acquisition + p2 + mu2
                  acquisition ~ p1 + mu1 + XylellaPopulation
                  p1 ~ pd_index
                  mu1 ~ pd_index
                  XylellaPopulation ~ resistanceTrait + week
                  pd_index ~ resistanceTrait + XylellaPopulation'

## Model including direct effects of attraction/leaving to/from source plants on transmission
transModel2 <- '# Regressions
                  transmission ~ acquisition + p2 + mu2 + p1 + mu1
                  acquisition ~ p1 + mu1 + XylellaPopulation
                  p1 ~ pd_index
                  mu1 ~ pd_index
                  XylellaPopulation ~ resistanceTrait
                  pd_index ~ resistanceTrait + XylellaPopulation'


## Model including latent variable "inoculation"
transModelInoc <- '# Regressions
                  transmission ~ acquisition + inoculation
                  acquisition ~ p1 + mu1 + XylellaPopulation
                  p1 ~ pd_index
                  mu1 ~ pd_index
                  XylellaPopulation ~ resistanceTrait
                  pd_index ~ resistanceTrait + XylellaPopulation
                  # Latent variable definition
                  inoculation =~ p2 + mu2
                  # Covariances
                  p1 ~~ mu2
                  p2 ~~ mu1'





##########################################################################################################
#### 2016 SEM analysis
##########################################################################################################

transdata <- readRDS("output/pdr1_transmission_preference_dataset.rds") %>% dplyr::filter(., week != 12.2)
transdata$source.cfu.per.g <- as.numeric(transdata$source.cfu.per.g)
transdata$test.plant.infection <- as.integer(transdata$test.plant.infection)
transdata$trt <- factor(transdata$trt)
str(transdata)

## Select only relevant variables
semData16 <- transdata %>% 
  mutate(log.source.cfu = log10(source.cfu.per.g+1),
         # Trt code: R = 1, S = 0; so a positive relationship between trt and infection status indicates greater trans from Resistant
         trtNumeric = 1 - (as.numeric(trt) - 1)) %>%
  dplyr::select(test.plant.infection, propVectorInfectious, log.source.cfu, pd_index, trtNumeric, week, p1, p2, mu1, mu2) %>%
  dplyr::filter(., complete.cases(.)) %>%
  ## Remove outlier for p1
  dplyr::filter(p1 < 7) %>%
  rename(transmission = test.plant.infection,
         acquisition = propVectorInfectious,
         resistanceTrait = trtNumeric,
         XylellaPopulation = log.source.cfu)
xData <- semData16 %>% 
  dplyr::select(-transmission) %>%
  scale(., scale = TRUE, center = TRUE)
semData16 <- cbind(semData16$transmission, xData) %>% 
  as.data.frame() %>%
  rename(transmission = V1)
semData16$transmission <- as.integer(semData16$transmission)
str(semData16)


fit16 <- sem(transModel, data = semData16, ordered = c("transmission", "pd_index"),
             optim.method = "nlminb", control = list(iter.max = 8000))


#### Compare model fits with covariances
fit16Inoc <- sem(transModelInoc, data = semData16, #ordered = c("transmission", "pd_index"), 
                optim.method = "BFGS")
## Model without inoculation/covariances is much better
AICtab(fit16, fit16Inoc, base = TRUE, sort = TRUE)


#### Compare model fits with direct links of p1 and mu1 to transmission
fit16mediation <- sem(transModel2, data = semData16, #ordered = c("transmission", "pd_index"), 
                      optim.method = "BFGS")
AICtab(fit16, fit16mediation, base = TRUE, sort = TRUE)
## No direct effect is best



#### Original model is best
## Checking model fit
summary(fit16, standardized = TRUE, fit.measures = TRUE)

## Plot best SEM model
semPaths(fit16, what = "std", fade = FALSE, residuals = FALSE)

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
  mutate(log.xfpop = log10(xfpop+1),
         # Trt code: R = 1, S = 0; so a positive relationship between trt and infection status indicates greater trans from Resistant
         trtNumeric = 1 - (as.numeric(trt) - 1)) %>%
  dplyr::select(test_plant_infection, trtNumeric, week, PD_symptoms_index, propInfectious, log.xfpop, mu1, mu2, p1, p2) %>%
  dplyr::filter(complete.cases(.)) %>%
  dplyr::filter(mu2 < 6.5) %>%
  rename(transmission = test_plant_infection,
         acquisition = propInfectious,
         resistanceTrait = trtNumeric,
         XylellaPopulation = log.xfpop,
         pd_index = PD_symptoms_index)
xData <- semData17 %>% 
  dplyr::select(-transmission) %>%
  scale(., scale = TRUE, center = TRUE)
semData17 <- cbind(semData17$transmission, xData) %>% 
  as.data.frame() %>%
  rename(transmission = V1)
str(semData17)

fit17 <- sem(transModel, data = semData17, ordered = c("transmission", "pd_index"),
             optim.method = "BFGS")


#### Compare model fits with covariances
fit17Inoc <- sem(transModelInoc, data = semData17, #ordered = c("transmission", "pd_index"), 
                 optim.method = "BFGS")
AICtab(fit17, fit17Inoc, base = TRUE, sort = TRUE)
## Model without inoculation/covariances is much better

#### Compare model fits with direct links between p1 and mu1 and transmission
fit17mediation <- sem(transModel2, data = semData17, #ordered = c("transmission", "pd_index"), 
                      optim.method = "BFGS")
AICtab(fit17, fit17mediation, base = TRUE, sort = TRUE)
## No direct effect is best
summary(fit17mediation, standardized = TRUE)


#### Results from best model
summary(fit17, standardized = TRUE, fit.measures = TRUE)



##########################################################################################################
##########################################################################################################
#### SEM analysis using piecewiseSEM
##########################################################################################################
##########################################################################################################
#### 2016 SEM analysis 
##########################################################################################################

transdata <- readRDS("output/pdr1_transmission_preference_dataset.rds") %>% dplyr::filter(., week != 12.2)
transdata$source.cfu.per.g <- as.numeric(transdata$source.cfu.per.g)
transdata$test.plant.infection <- as.integer(transdata$test.plant.infection)
transdata$trt <- factor(transdata$trt)
str(transdata)

## Select only relevant variables
semData16 <- transdata %>% 
  mutate(acquisition = asin(sqrt(propVectorInfectious)),
         # Trt code: R = 1, S = 0; so a positive relationship between trt and infection status indicates greater trans from Resistant
         trtNumeric = 1 - (as.numeric(trt) - 1),
         XylellaPopulation = log(source.cfu.per.g+1),
         pd_index = log(pd_index+1)) %>%
  dplyr::select(test.plant.infection, acquisition, XylellaPopulation, pd_index, trtNumeric, week, p1, p2, mu1, mu2) %>%
  dplyr::filter(., complete.cases(.)) %>%
  ## Remove outlier for p1
  dplyr::filter(p1 < 7) %>%
  rename(transmission = test.plant.infection,
         resistanceTrait = trtNumeric)
         #XylellaPopulation = source.cfu.per.g)
str(semData16)


#### Run model selection for 2016
modelSelect16 <- selectSEM(semData16)
(results16 <- modelSelect16[[1]])

## Best model is No Preference Model
modelSelect16[[2]]$attractModel2
coefTable16 <- modelSelect16[[2]]$attractModel2$coefficients

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
         # Trt code: R = 1, S = 0; so a positive relationship between trt and infection status indicates greater trans from Resistant
         trtNumeric = 1 - (as.numeric(trt) - 1),
         XylellaPopulation = log(xfpop+1),
         pd_index = log(PD_symptoms_index+1)) %>%
  dplyr::select(test_plant_infection, trtNumeric, week, pd_index, acquisition, XylellaPopulation, mu1, mu2, p1, p2) %>%
  dplyr::filter(complete.cases(.)) %>%
  dplyr::filter(mu2 < 6.5) %>%
  rename(transmission = test_plant_infection,
         resistanceTrait = trtNumeric)
         #XylellaPopulation = xfpop,
         #pd_index = PD_symptoms_index)
str(semData17)

modelSelect17 <- selectSEM(semData17)
(results17 <- modelSelect17[[1]])

## Best model is attractModel2; attraction rates without pd_index
modelSelect17[[2]]$attractModel2
