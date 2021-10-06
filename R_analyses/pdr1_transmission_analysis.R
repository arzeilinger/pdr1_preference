#### ANALYSIS OF PDR1 TRANSMISSION DATA

rm(list = ls())
# Load packages
my.packages <- c("data.table", "openxlsx", "MASS", "googlesheets",
                 "lme4", "ggplot2", "bbmle", "multcomp", "DHARMa", "cowplot", "broom",
                 "rms", "tidyr", "dplyr")
lapply(my.packages, require, character.only = TRUE)

## Load functions
source("R_functions/factor2numeric.R")
# Load Negative Log-Likelihood functions for all models
source("R_functions/nll_functions_pdr1.R")


## Specifications for ggplot2 text size
baseTextSize <- 16


######################################################################################################################
######################################################################################################################
#### Analysis of 2016 data
######################################################################################################################

#### Load combined and filtered transmission/preference data set
# Remove second 12-week set of trials as these were designed differently
transdata <- readRDS("output/pdr1_transmission_preference_dataset.rds") %>% dplyr::filter(., week != 12.2)
transdata$source.cfu.per.g <- as.numeric(transdata$source.cfu.per.g)
transdata$test.plant.infection <- as.integer(transdata$test.plant.infection)
transdata$trt <- factor(transdata$trt)
str(transdata)

######################################################################################################################
#### Analysis of PD symptom data
# Ordered logistic regression
# also called partial odds logistic regression
transdata$pd_index <- factor(transdata$pd_index, ordered = TRUE, levels = c("0", "1", "2", "3", "4", "5"))
transdataDD <- transdata %>% dplyr::select(pd_index, week, trt)
dd <- datadist(transdataDD)
options(datadist = "dd")
lrmMod16 <- lrm(pd_index ~ week*trt, data = transdataDD)
lrmMod16
summary(lrmMod16, conf.type = "simultaneous")

#############################################################################################################
#### source xf pop
## Data look overdispesed; considered negative binomial but quasi-poisson performs better
## Quasi-poisson
sourcexfQP <- glm(source.cfu.per.g ~ week*trt, data = transdata, family = "quasipoisson")
plot(sourcexfQP)
summary(sourcexfQP)
tidy(sourcexfQP, conf.int = TRUE)
## Even though the error variance is highly skewed, model fit seems best among different options


#######################################################################################################
#### Plotting
transSummary <- transdata %>% 
  mutate(pd_index_num = factor2numeric(pd_index)) %>%
  group_by(., week, trt) %>% 
  summarise(propInfected = sum(test.plant.infection, na.rm = TRUE)/sum(!is.na(test.plant.infection)),
            meanPD = mean(pd_index_num, na.rm = TRUE),
            sePD = sd(pd_index_num, na.rm = TRUE)/sqrt(sum(!is.na(pd_index_num))),
            meancfu = mean(log10(source.cfu.per.g+1)),
            secfu = sd(log10(source.cfu.per.g+1))/sqrt(sum(!is.na(source.cfu.per.g))))


#### Symptoms plot
symptom16Plot <- ggplot(data=transSummary, aes(x=week, y=meanPD)) +
  geom_line(aes(linetype=trt), size=1.25) +
  geom_point(aes(shape=trt), size=2.5) +
  geom_errorbar(aes(ymax=meanPD+sePD, ymin=meanPD-sePD), width=0.2) +
  scale_x_continuous(name = "", 
                     breaks = c(3,8,12), limits = c(0,12)) + 
  scale_y_continuous(name = "PD severity index",
                     breaks = seq(0,5,by=1), limits = c(0, 5)) +
  scale_shape_manual(values = c(16,1)) +
  theme_bw(base_size=baseTextSize) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none") 

symptom16Plot


#### Xf pops in source plant plot
sourcexf16plot <- ggplot(data=transSummary, aes(x=week, y=meancfu)) +
  geom_line(aes(linetype=trt), size=1.25) +
  geom_point(aes(shape=trt), size=2.5) +
  geom_errorbar(aes(ymax=meancfu+secfu, ymin=meancfu-secfu), width=0.2) +
  scale_x_continuous(name = "", 
                     breaks = c(3,8,12), limits = c(0,12)) + 
  scale_y_continuous(name = "Xf pop. size (log10, CFU/g)",
                     breaks = seq(2,10,by=2), limits = c(2,10)) +
  scale_shape_manual(values = c(16,1)) +
  theme_bw(base_size=baseTextSize) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none") 

sourcexf16plot



##############################################################################################################
#### Fitting non-linear models to acquisition data
##############################################################################################################


#### Initial parameters
## Initial parameters for linear model
a.l0 <- 0 # Y-intercept assumed at 0
b.l0 <- 0.6/12 # Looks like the line might hit 50% at 14 weeks, calculate initial slope from this

## Initial parameters for Holling Type IV
a0 <- 0.1 # Asymptotic totalInfectious
# Peak transmission is at 8 weeks, peak = -2b/c, b/c = -4, and c should be < 0
b0 <- 8
c0 <- -2

## Initial parameters for Ricker model
a.rick0 <- 0.5/3 # Initial slope of the line
b.rick0 <- 1/8 # x value for the peak

## Initial parameters for Logistic Growth model
b.logist0 <- 0.5/3/4 # Initial slope of the line divided by 4
a.logist0 <- -4*b.logist0 # Negative value of x at halfway to asymptote (inflection point) multiplied by b

## Initial parameters for Michaelis-Menten model
a.mm0 <- 1 # Asymptote
b.mm0 <- 4 # value of x when y = a/2


#### Fitting multiple non-linear models to Resistant and Susceptible trts separately
## Load acquisition data
acqDataVector <- readRDS("output/pdr1_2016_vector_acquisition_dataset.rds")
#acqDataVector %>% ungroup() %>% write.csv(., "output/pdr1_2016_vector_acquisition_dataset.csv", row.names = FALSE)
str(acqDataVector)

## Summarize data by week and trt
acqSummary16 <- acqDataVector %>% ungroup() %>%
  dplyr::filter(!is.na(vectorInfectious)) %>%
  group_by(trt, week) %>% 
  summarise(n = length(vectorInfectious),
            nInfected = sum(vectorInfectious),
            propInfected = nInfected/n)
## Split into Resistant and Sucsceptible data sets
acqR <- acqSummary16 %>% dplyr::filter(trt == "R")
acqS <- acqSummary16 %>% dplyr::filter(trt == "S")

#### Fitting resistant line 2016 data
vectorRresults16 <- optimizeTransModels(acqR, nbugs = acqR$n)
(modelSelectR16 <- vectorRresults16[[2]])
opListR <- vectorRresults16[[1]]

# Get model predictions for plotting
bestModR <- opListR$rickerOp # Ricker model is the best
bestcoefR <- as.list(coef(bestModR))

## Estimate confidence intervals for model parameter estimates using profile method and combine data
ciR16 <- confint(bestModR)
## Using profile method 
acqParamR16 <- data.frame(dataset = "R16",
                          model = "Ricker",
                          coef = names(coef(bestModR)),
                          estimate = coef(bestModR),
                          cil = ciR16[,1],
                          ciu = ciR16[,2])

newDatAcqR16 <- data.frame(week = seq(0, 12, length = 201))
newDatAcqR16$nInfected <- with(newDatAcqR16, bestcoefR$a*week*exp(-bestcoefR$b*week))


#### Fitting susceptible line 2016 data
vectorSresults16 <- optimizeTransModels(acqS, nbugs = acqS$n)
(modelSelectS16 <- vectorSresults16[[2]])
opListS <- vectorSresults16[[1]]

# Get model predictions for plotting
bestModS <- opListS$rickerOp # Ricker model is the best
bestcoefS <- as.list(coef(bestModS))

## Estimate confidence intervals for model parameter estimates using profile method and combine data
ciS16 <- confint(bestModS) 
acqParamS16 <- data.frame(dataset = "S16",
                          model = "Ricker",
                          coef = names(coef(bestModS)),
                          estimate = coef(bestModS),
                          cil = ciS16[,1],
                          ciu = ciS16[,2])

newDatAcqS16 <- data.frame(week = seq(0, 12, length = 201))
newDatAcqS16$nInfected <- with(newDatAcqS16, bestcoefS$a*week*exp(-bestcoefS$b*week))


#### Plotting results
## Calculate mean infectiousness for each trt and week
vectorSummary16 <- acqDataVector %>% ungroup() %>% 
  dplyr::filter(!is.na(vectorInfectious)) %>%
  group_by(trt, week, rep) %>% summarise(propInfectious = sum(vectorInfectious)/length(vectorInfectious)) %>%
  group_by(trt, week) %>% summarise(meanPropInfectious = mean(propInfectious),
                                    nPropInfectious = length(propInfectious),
                                    sePropInfectious = sd(propInfectious)/sqrt(nPropInfectious)) 
vectorSummary16

## Plot S and R genotypes summarised together
## Vector acquisition plot black and white
vectorplotNL16 <- ggplot(data=vectorSummary16, aes(x=week, y=meanPropInfectious)) +
  # R = Closed circles and solid line
  # S = Open circles and dashed line
  geom_point(aes(shape = trt), size=2) +
  geom_errorbar(aes(ymax=meanPropInfectious+sePropInfectious, ymin=meanPropInfectious-sePropInfectious), width=0.2) +
  geom_line(data = newDatAcqR16, aes(x = week, y = nInfected), linetype = 1, colour = "black", lwd = 1.25) +
  geom_line(data = newDatAcqS16, aes(x = week, y = nInfected), linetype = 2, colour = "black", lwd = 1.25) +
  scale_x_continuous(name = "", 
                     breaks = c(3,8,12), limits = c(0,12)) + 
  scale_y_continuous(name = "Prop. vectors infected",
                     breaks = seq(0,1,0.2), limits = c(-0.01,1.01)) +
  scale_shape_manual(values = c(16,1)) +
  theme_bw(base_size=baseTextSize) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none") 

vectorplotNL16



##############################################################################################################
#### Non-linear models of transmission
##############################################################################################################

#### Summarize transmission: mean transmission between trt
transdata %>% group_by(trt) %>% summarise(mean = mean(test.plant.infection, na.rm = TRUE))


# Setting up data
transSummarynl16 <- transdata %>% group_by(week, trt) %>% 
  summarise(n = length(test.plant.infection),
            nInfected = sum(test.plant.infection, na.rm = TRUE),
            propInfected = nInfected/n)
transSummarynl16


#### Initial parameters
## Initial parameters for linear model
a.l0 <- 1 # Y-intercept assumed at 0
b.l0 <- -0.5/8 # Looks like the line might hit 50% at 14 weeks, calculate initial slope from this

## Initial parameters for Holling Type IV
a0 <- 0.1 # Asymptotic propInfected
# Assume peak transmission is at 8 weeks, peak = -2b/c, b/c = -4, and c should be < 0
b0 <- 8
c0 <- -2

## Initial parameters for Ricker model
a.rick0 <- 0.75*(1/8)/exp(-1) # Estimate of a from estimate of b and peak transmission
b.rick0 <- 1/8 # x value for the peak

## Initial parameters for Logistic Growth model
b.logist0 <- (0.125-0.75)/(12-3)/4 # Initial slope of the line divided by 4
a.logist0 <- -6*b.logist0 # Negative value of x at halfway to asymptote (inflection point) multiplied by b

## Initial parameters for Michaelis-Menten model
a.mm0 <- 0.25 # Asymptote
b.mm0 <- 8 # value of x when y = a/2


##############################################################################################################
#### Fitting multiple non-linear models to Resistant and Susceptible trts separately

# Get data sets
transR <- transSummarynl16[transSummarynl16$trt == "R",]
transS <- transSummarynl16[transSummarynl16$trt == "S",]


#### Fitting resistant line data
transRresults16 <- optimizeTransModels(dat = transR, nbugs = 8)
(modelSelectR16 <- transRresults16[[2]])
opListR <- transRresults16[[1]]

# Get model predictions for plotting
bestModR <- opListR$rickerOp 
# Logistic, Ricker, and Linear are all equivalent; only Ricker makes sense if we assume dynamics start at origin
bestcoefR <- as.list(coef(bestModR))

## Estimate confidence intervals for model parameter estimates using profile method and combine data
ciR16 <- confint(bestModR)
paramR16 <- data.frame(dataset = "R16",
                       model = "Ricker",
                       coef = row.names(ciR16),
                       estimate = coef(bestModR),
                       cil = ciR16[,1],
                       ciu = ciR16[,2])

## Predict propInfected from model
newDatTransR16 <- data.frame(week = seq(0, 12, length = 201))
newDatTransR16$propInfected <- with(newDatTransR16, bestcoefR$a*week*exp(-bestcoefR$b*week))

#### Fitting susceptible line data
transSresults16 <- optimizeTransModels(dat = transS, nbugs = 8)
(modelSelectS16 <- transSresults16[[2]])
opListS <- transSresults16[[1]]

## Get model predictions for plotting
bestModS <- opListS$rickerOp # Ricker model is the best
bestcoefS <- as.list(coef(bestModS))

## Estimate confidence intervals for model parameter estimates using profile method and combine data
ciS16 <- confint(bestModS)
paramS16 <- data.frame(dataset = "S16",
                       model = "Ricker",
                       coef = row.names(ciS16),
                       estimate = coef(bestModS),
                       cil = ciS16[,1],
                       ciu = ciS16[,2])
nlParamTable <- rbind(paramR16, paramS16)

## Predict propInfected from model
newDatTransS16 <- data.frame(week = seq(0, 12, length = 201))
newDatTransS16$propInfected <- with(newDatTransS16, bestcoefS$a*week*exp(-bestcoefS$b*week))


#### Plotting results
# Plot S and R genotypes summarised together
# Transmission plot
transplotNL16 <- ggplot(data=transSummarynl16, aes(x=week, y=propInfected)) +
  # R = Circles and solid line
  # S = Triangles and dashed line
  geom_point(aes(shape = trt), size=2) +
  geom_line(data = newDatTransR16, aes(x = week, y = propInfected), linetype = 1, colour = "black", lwd = 1.25) +
  geom_line(data = newDatTransS16, aes(x = week, y = propInfected), linetype = 2, colour = "black", lwd = 1.25) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = c(3,8,12), limits = c(0,12)) + 
  scale_y_continuous(name = "Prop. test plants infected",
                     breaks = seq(0,1,by=0.2), limits = c(0,1)) +
  scale_shape_manual(values = c(16, 1)) +
  theme_bw(base_size=baseTextSize) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none")

transplotNL16



##############################################################################################################
##############################################################################################################
#### 2017 data
##############################################################################################################

transVCPdata <- readRDS("output/complete_2017_transmission-preference_dataset.rds")
str(transVCPdata)
with(transVCPdata, table(week, genotype))


##############################################################################################################
#### ## Analysis of PD symptoms using Partial Odds Logistic Regression
transVCPdata$PD_symptoms_index <- factor(transVCPdata$PD_symptoms_index, ordered = TRUE, levels = c("0", "1", "2", "3", "4", "5"))
pdData <- transVCPdata %>% dplyr::select(PD_symptoms_index, block, week, genotype) %>% dplyr::filter(complete.cases(.))
dd <- datadist(pdData)
options(datadist = "dd")
lrmMod17 <- lrm(PD_symptoms_index ~ block + week*genotype, data = pdData)
lrmMod17
summary(lrmMod17, conf.type = "simultaneous")

##############################################################################################################
#### Analysis of Xf pops in source plants
## Generalized linear model with quasipoisson distribution
poispopMod1 <- glm(xfpop ~ block + week*genotype, data = transVCPdata, family = "quasipoisson")
plot(poispopMod1)
# Residuals look about the same as lm() with sqrt()
summary(poispopMod1)
tidy(poispopMod1, conf.int = TRUE)

## Negative binomial GLM
sourcexfNB <- glm.nb(xfpop ~ block + week*genotype, data = transVCPdata)
summary(sourcexfNB)

#### NOTE ON ALTERNATIVE MODELS: Removing false negatives and all negatives had negligible effects on results.
#### These models had reduced significance of week main effect but were otherwise unchanged.
#### Log10 transformation works better if all negatives are removed
#### Results from negative binomial are qualitatively similar to quasipoisson ...
#### but warning suggests poor fitting. Quasipoisson also works better for both years.



##############################################################################################################
#### Plotting PD symptoms, Xf pops, transmission, and acquisition from 2017

#### Plotting PD symptoms
## Summary
pdSummary <- transVCPdata %>% group_by(week, genotype, trt) %>% 
  mutate(pd_index_num = factor2numeric(PD_symptoms_index)) %>%
  summarise(meanPD = mean(pd_index_num, na.rm = TRUE), 
            sePD = sd(pd_index_num, na.rm = TRUE)/sqrt(sum(!is.na(pd_index_num))))

## PD symptoms plot
symptom17plot <- ggplot(data=pdSummary, aes(x=week, y=meanPD)) +
  geom_line(aes(linetype=genotype), colour = "black", size=1.25) +
  geom_point(aes(shape=genotype), colour = "black", size=2.5) +
  geom_errorbar(aes(ymax=meanPD+sePD, ymin=meanPD-sePD), width=0.2) +
  scale_x_continuous(name = "", 
                     breaks = c(2,5,8,14), limits = c(0,14)) + 
  scale_y_continuous(name = "",
                     limits = c(0,5)) +
  # 007 = open circles, dashed line
  # 092 = open triangles, dashed line
  # 094 = closed circles, solid line
  # 102 = closed triangles, solid line
  scale_shape_manual(NULL, values = c(1,2,16,17)) +
  scale_linetype_manual(NULL, values = c(2,2,1,1)) +
  theme_bw(base_size=baseTextSize) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank(),
        legend.justification = c(0,0),
        legend.position = c(0.1, 0.3))

symptom17plot


#### Plotting source plant xf populations
## Summary
sourceSummary <- transVCPdata %>% 
  dplyr::filter(!is.na(xfpop)) %>%
  mutate(logxfpop = log10(xfpop+1), sqrtxfpop = sqrt(xfpop)) %>% 
  group_by(week, genotype, trt) %>% 
  summarise_at(c("logxfpop", "sqrtxfpop"), list(mean = ~mean(.), n = ~sum(.), se = ~sd(.)/sqrt(length(.))))

## Xf pops in source plant plot
sourcexf17plot <- ggplot(data=sourceSummary, aes(x=week, y=logxfpop_mean)) +
  geom_line(aes(linetype=genotype), colour = "black", size=1.25) +
  geom_point(aes(shape=genotype), colour = "black", size=2.5) +
  geom_errorbar(aes(ymax=logxfpop_mean+logxfpop_se, ymin=logxfpop_mean-logxfpop_se), width=0.2) +
  scale_x_continuous(name = "", 
                     breaks = c(2,5,8,14), limits = c(0,14)) + 
  scale_y_continuous(name = "",
                     breaks = seq(2,10,by=2), limits = c(2,10)) +
  # 007 = open circles, dashed line
  # 092 = open triangles, dashed line
  # 094 = closed circles, solid line
  # 102 = closed triangles, solid line
  scale_shape_manual(values = c(1,2,16,17)) +
  scale_linetype_manual(values = c(2,2,1,1)) +
  theme_bw(base_size=baseTextSize) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none") 

sourcexf17plot


##############################################################################################################
#### Fitting non-linear models to acquisition data
##############################################################################################################

## Remove trials where pre-screen plant was positive and test plant was positive
badTrials <- with(transVCPdata, which(grepl("PS plant positive", notes) & test_plant_infection == 1))
transdata2 <- transVCPdata[-badTrials,]
with(transdata2, table(week, genotype)) # Sample size by week and genotype


#### Initial parameters
## Initial parameters for linear model
a.l0 <- 0 # Y-intercept assumed at 0
b.l0 <- 0.5/8 # Looks like the line might hit 50% at 14 weeks, calculate initial slope from this

## Initial parameters for Holling Type IV
a0 <- 0.1 # Asymptotic totalInfectious
# Peak transmission is at 8 weeks, peak = -2b/c, b/c = -4, and c should be < 0
b0 <- 8
c0 <- -2

## Initial parameters for Ricker model
a.rick0 <- 0.1/2 # Initial slope of the line
b.rick0 <- 1/8 # x value for the peak

## Initial parameters for Logistic Growth model
b.logist0 <- 0.1/2/4 # Initial slope of the line divided by 4
a.logist0 <- -4*b.logist0 # Negative value of x at halfway to asymptote (inflection point) multiplied by b

## Initial parameters for Michaelis-Menten model
a.mm0 <- 1 # Asymptote
b.mm0 <- 8 # value of x when y = a/2


#### Fitting multiple non-linear models to Resistant and Susceptible trts separately

# Get data sets
acqSummary17 <- transdata2 %>% group_by(trt, week) %>% summarise(n = sum(nbugs),
                                                                 nInfected = sum(totalInfectious),
                                                                 propInfected = nInfected/n)
acqR <- acqSummary17 %>% dplyr::filter(trt == "R")
acqS <- acqSummary17 %>% dplyr::filter(trt == "S")

#### Fitting resistant line data
vectorRresults17 <- optimizeTransModels(acqR, nbugs = acqR$n)
(modelSelectR17 <- vectorRresults17[[2]])
opListR <- vectorRresults17[[1]]

# Get model predictions for plotting
bestModR <- opListR$holling4Op # Holling4 model is the best
bestcoefR <- as.list(coef(bestModR))

## Estimate confidence intervals for model parameter estimates using profile method and combine data
ciR17 <- confint(bestModR, method = "quad")
acqParamR17 <- data.frame(dataset = "R17",
                       model = "Holling_type_IV",
                       coef = names(coef(bestModR)),
                       estimate = coef(bestModR),
                       cil = ciR17[,1],
                       ciu = ciR17[,2])

newDatAcqR17 <- data.frame(week = seq(0, 14, length = 201))
newDatAcqR17$nInfected <- with(newDatAcqR17, (bestcoefR$a*week^2)/(bestcoefR$b + bestcoefR$c*week + week^2))


#### Fitting susceptible line data
vectorSresults17 <- optimizeTransModels(acqS, nbugs = acqS$n)
(modelSelectS17 <- vectorSresults17[[2]])
opListS <- vectorSresults17[[1]]

# Get model predictions for plotting
bestModS <- opListS$holling4Op # Ricker model is the best
bestcoefS <- as.list(coef(bestModS))

## Estimate confidence intervals for model parameter estimates using profile method and combine data
ciS17 <- confint(bestModS, method = "quad")
acqParamS17 <- data.frame(dataset = "S17",
                          model = "Holling_type_IV",
                          coef = names(coef(bestModS)),
                          estimate = coef(bestModS),
                          cil = ciS17[,1],
                          ciu = ciS17[,2])

newDatAcqS17 <- data.frame(week = seq(0, 14, length = 201))
newDatAcqS17$nInfected <- with(newDatAcqS17, (bestcoefS$a*week^2)/(bestcoefS$b + bestcoefS$c*week + week^2))


#### Plotting results
## Calculate mean infectiousness for each trt and week
vectorSummary2 <- transdata2 %>% dplyr::filter(!is.na(propInfectious)) %>% 
  group_by(trt, week) %>% summarise(meanpropInfectious = mean(propInfectious, na.rm = TRUE),
                                    npropInfectious = length(propInfectious),
                                    sepropInfectious = sd(propInfectious, na.rm = TRUE)/sqrt(npropInfectious)) 
vectorSummary2

# Plot S and R genotypes summarised together
# Use Holling 4 model for S and for R
# Vector acquisition plot
vectorplotNL17 <- ggplot(data=vectorSummary2, aes(x=week, y=meanpropInfectious)) +
  # R = Closed circles and solid line
  # S = Open circles and dashed line
  geom_point(aes(shape = trt), size=2) +
  geom_errorbar(aes(ymax=meanpropInfectious+sepropInfectious, ymin=meanpropInfectious-sepropInfectious), width=0.2) +
  geom_line(data = newDatAcqR17, aes(x = week, y = nInfected), linetype = 1, colour = "black", lwd = 1.25) +
  geom_line(data = newDatAcqS17, aes(x = week, y = nInfected), linetype = 2, colour = "black", lwd = 1.25) +
  scale_x_continuous(name = "", 
                     breaks = c(2,5,8,14), limits = c(0,14)) + 
  scale_y_continuous(name = "",
                     breaks = seq(0,1,0.2), limits = c(-0.01,1.01)) +
  scale_shape_manual(values = c(16,1)) +
  theme_bw(base_size=baseTextSize) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none") 

vectorplotNL17


#############################################################################################
#### Combining non-linear acquisition model results from both years for paper

#### Combine results from analyses
nlModelSelectionList <- list(modelSelectR16,
                             modelSelectS16,
                             modelSelectR17,
                             modelSelectS17)
names(nlModelSelectionList) <- c("Resistant-2016", "Susceptible-2016", "Resistant-2017", "Susceptible-2017")

## Save parameter estimates to one sheet and all model selection tables to a second sheet of the same workbook
wb <- createWorkbook()
addWorksheet(wb, sheetName = "nl_model_selection_tables")
startRows <- seq(1, length(nlModelSelectionList)*7)

for(i in 1:length(nlModelSelectionList)){
  modelSelectionDF <- data.frame(modelName = attr(nlModelSelectionList[[i]], "row.names"),
                                 AICc = round(nlModelSelectionList[[i]]$AICc, 2),
                                 dAICc = round(nlModelSelectionList[[i]]$dAICc, 2),
                                 df = nlModelSelectionList[[i]]$df)
  modelSelectionDF[,"modelName"] <- with(modelSelectionDF, ifelse(modelName == "rickerOp", "Ricker",
                                                                  ifelse(modelName == "holling4Op", "Holling Type IV",
                                                                         ifelse(modelName == "linearOp", "Linear",
                                                                                ifelse(modelName == "MMOp", "Michaelis-Menten",
                                                                                       "Logistic Growth")))))
  writeData(wb, sheet = "nl_model_selection_tables", x = modelSelectionDF, startCol = 2, startRow = i*7)
  writeData(wb, sheet = "nl_model_selection_tables", x = names(nlModelSelectionList)[i], startCol = 1, startRow = i*7)
}

## Parameter estimates
nlParamTable <- rbind(acqParamR16, acqParamS16, acqParamR17, acqParamS17)
nlParamTable <- nlParamTable %>% mutate(estimateCI = paste(round(estimate, 3), " [", round(cil, 3), ", ", round(ciu,3), "]", sep = "")) %>%
  mutate("Genotype-Year" = ifelse(dataset == "R16", "Resistant-2016",
                                  ifelse(dataset == "S16", "Susceptible-2016",
                                         ifelse(dataset == "R17", "Resistant-2017",
                                                "Susceptible-2017")))) %>%
  dplyr::select(-estimate, -cil, -ciu, -dataset)
addWorksheet(wb, sheetName = "nl_parameter_estimates")
writeData(wb, sheet = "nl_parameter_estimates", x = nlParamTable, startCol = 2, startRow = 2)

saveWorkbook(wb, file = "results/nl_model_vector_acquisition_results_tables_both_years.xlsx", overwrite = TRUE)



##############################################################################################################
#### Fitting non-linear models to transmission data
##############################################################################################################

#### Summarize transmission: mean transmission between trt
transVCPdata %>% group_by(trt) %>% summarise(mean = mean(test_plant_infection, na.rm = TRUE))


# Setting up data
transSummarynl17 <- transdata2 %>% group_by(week, trt) %>% 
  summarise(n = length(test_plant_infection),
            nInfected = sum(test_plant_infection),
            propInfected = nInfected/n)



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


##############################################################################################################
#### Fitting multiple non-linear models to Resistant and Susceptible trts separately

# Get data sets
transR <- transSummarynl17[transSummarynl17$trt == "R",]
transS <- transSummarynl17[transSummarynl17$trt == "S",]


#### Fitting resistant line data
transRresults17 <- optimizeTransModels(dat = transR, nbugs = 16)
(modelSelectR17 <- transRresults17[[2]])
opListR <- transRresults17[[1]]

# Get model predictions for plotting
bestModR <- opListR$holling4Op # Holling4 model is the best
bestcoefR <- as.list(coef(bestModR))

## Estimate confidence intervals for model parameter estimates using profile method and combine data
ciR17 <- confint(bestModR, method = "quad")

## Profile method doesn't work, try quadratic approximation
paramR17 <- data.frame(dataset = "R17",
                       model = "Holling_type_IV",
                       coef = names(coef(bestModR)),
                       estimate = coef(bestModR),
                       cil = ciR17[,1],
                       ciu = ciR17[,2])

## Predict propInfected from model
newDatTransR17 <- data.frame(week = seq(0, 14, length = 201))
newDatTransR17$propInfected <- with(newDatTransR17, (bestcoefR$a*week^2)/(bestcoefR$b + bestcoefR$c*week + week^2))


#### Fitting susceptible line data
transSresults17 <- optimizeTransModels(dat = transS, nbugs = 16)
(modelSelectS17 <- transSresults17[[2]])
opListS <- transSresults17[[1]]

# Get model predictions for plotting
bestModS <- opListS$rickerOp # Ricker model is the best
bestcoefS <- as.list(coef(bestModS))

## Estimate confidence intervals for model parameter estimates using profile method and combine data
ciS17 <- confint(bestModS)
paramS17 <- data.frame(dataset = "S17",
                       model = "Ricker",
                       coef = row.names(ciS17),
                       estimate = coef(bestModS),
                       cil = ciS17[,1],
                       ciu = ciS17[,2])

## Predict propInfected from model
newDatTransS17 <- data.frame(week = seq(0, 14, length = 201))
newDatTransS17$propInfected <- with(newDatTransS17, bestcoefS$a*week*exp(-bestcoefS$b*week))


#### Plotting results
# Plot S and R genotypes summarised together
# Use Ricker model for S and Holling 4 model for R
## 2017 transmission plot black and white
transplotNL17 <- ggplot(data=transSummarynl17, aes(x=week, y=propInfected)) +
  # R = Circles and solid line
  # S = Triangles and dashed line
  geom_point(aes(shape = trt), size=2) +
  geom_line(data = newDatTransR17, aes(x = week, y = propInfected), linetype = 1, colour = "black", lwd = 1) +
  geom_line(data = newDatTransS17, aes(x = week, y = propInfected), linetype = 2, colour = "black", lwd = 1) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = c(2,5,8,14), limits = c(0,14)) + 
  scale_y_continuous(name = "",
                     breaks = seq(0,1,by=0.2), limits = c(0,1)) +
  scale_shape_manual(values = c(16,1)) +
  theme_bw(base_size=baseTextSize) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none") 

transplotNL17



##############################################################################
##############################################################################
#### Combining figures for paper
##############################################################################


##############################################################################
#### Combining PD symptom, Xf pop, and NL transmission results for paper for each year
#### Giant 8-panel figure

giant_figure <- plot_grid(symptom16Plot, symptom17plot, sourcexf16plot, sourcexf17plot,
                          vectorplotNL16, vectorplotNL17, transplotNL16, transplotNL17,
                          align = "h", ncol = 2, nrow = 4, 
                          labels = paste("(", LETTERS[1:8], ")", sep = ""),
                          label_x = 0.22, label_y = 0.94,
                          label_size = 14)

ggsave(filename = "results/figures/pd_symptom_xf_pop_nl_transmission_figure.tiff",
       plot = giant_figure,
       width = 16, height = 32, units = "cm", dpi = 600, compression = "lzw")



#### Combine results from analyses for manuscript
nlModelSelectionList <- list(modelSelectR16,
                             modelSelectS16,
                             modelSelectR17,
                             modelSelectS17)
names(nlModelSelectionList) <- c("R16", "S16", "R17", "S17")

## Save parameter estimates to one sheet and all model selection tables to a second sheet of the same workbook
wb <- createWorkbook()
addWorksheet(wb, sheetName = "nl_model_selection_tables")
startRows <- seq(1, length(nlModelSelectionList)*7)

for(i in 1:length(nlModelSelectionList)){
  modelSelectionDF <- data.frame(modelName = attr(nlModelSelectionList[[i]], "row.names"),
                                 AICc = round(nlModelSelectionList[[i]]$AICc, 2),
                                 dAICc = round(nlModelSelectionList[[i]]$dAICc, 2),
                                 df = nlModelSelectionList[[i]]$df)
  writeData(wb, sheet = "nl_model_selection_tables", x = modelSelectionDF, startCol = 2, startRow = i*7)
  writeData(wb, sheet = "nl_model_selection_tables", x = names(nlModelSelectionList)[i], startCol = 1, startRow = i*7)
}

## Parameter estimates
nlParamTable <- rbind(paramR16, paramS16, paramR17, paramS17)
nlParamTable <- nlParamTable %>% mutate(estimateCI = paste(round(estimate, 3), " [", round(cil, 3), ", ", round(ciu,3), "]", sep = "")) %>%
  dplyr::select(-estimate, -cil, -ciu)
addWorksheet(wb, sheetName = "nl_parameter_estimates")
writeData(wb, sheet = "nl_parameter_estimates", x = nlParamTable, startCol = 2, startRow = 2)

saveWorkbook(wb, file = "results/nl_model_transmission_results_tables_both_years.xlsx", overwrite = TRUE)

