#### ANALYSIS OF PDR1 TRANSMISSION DATA

rm(list = ls())
# Load packages
my.packages <- c("tidyr", "dplyr", "data.table", "openxlsx", "MASS", "googlesheets",
                 "lme4", "ggplot2", "bbmle", "multcomp")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/factor2numeric.R")


######################################################################################################################
######################################################################################################################
#### Analysis of 2016 data
######################################################################################################################

#### Analysis of PD symptom data
# Remove second 12-week set of trials
transdata <- readRDS("output/pdr1_transmission_preference_dataset.rds") %>% dplyr::filter(., week != 12.2)
transdata$source.cfu.per.g <- as.numeric(transdata$source.cfu.per.g)
transdata$test.plant.infection <- as.integer(transdata$test.plant.infection)

#### symptom data
sympMod <- glm(pd_index ~ week*trt*log10(source.cfu.per.g+1), data = transdata, family = "poisson")
sympTrtMod <- glm(pd_index ~ week*trt, data = transdata, family = "poisson")
sympNoInterxnMod <- glm(pd_index ~ week + trt, data = transdata, family = "poisson")
AICctab(sympMod, sympTrtMod, sympNoInterxnMod, base = TRUE, delta = TRUE)
# sympTrtMod (with interaction) is best and quasipoisson distribution is best, AIC won't work with quasipoisson
sympTrtMod <- glm(pd_index ~ week*trt*log10(source.cfu.per.g+1), data = transdata, family = "quasipoisson")
plot(sympTrtMod)
summary(sympTrtMod)
# Linear model with transformed symptom data
boxcox(pd_index+1 ~ week*trt*log10(source.cfu.per.g+1), data = transdata)
sympLM <- lm(log(pd_index+1) ~ week*trt*log10(source.cfu.per.g+1), data = transdata)
plot(sympLM)
summary(sympLM)
# Ordered logistic regression
# also called partial odds logistic regression
transdata$pd_index <- factor(transdata$pd_index, ordered = TRUE, levels = c("0", "1", "2", "3", "4", "5"))
olrMod <- polr(pd_index ~ week*trt, data = transdata, Hess = TRUE, method = "logistic")
summary(olrMod)
# Calculate p values from t statistic
olrCoefs <- coef(summary(olrMod))
p <- pnorm(abs(olrCoefs[, "t value"]), lower.tail = FALSE)*2
(olrCoefs <- cbind(olrCoefs, "p-value" = p))
# Results: quasi-Poisson model, transformed LM model, and POLR model all give same result -> only week is significant positive 






#### Analyzing acquisition data at the per-vector level, with cage as a random effect
acqDataVector$cage <- with(acqDataVector, paste(week, trt, rep, sep = ""))
acqDataVector$vectorcfu <- floor(acqDataVector$vectorcfu)

# Analysis of CFU data
acqMod <- glmer(vectorcfu ~ week*trt + (1|cage),
                data = acqDataVector, family = "poisson")
summary(acqMod)

# Analysis of infection status
acqDataVector$vectorInfected <- ifelse(acqDataVector$vectorcfu > 0, 1, 0) # Turn CFUs into binomial presence/absence
infectedMod <- glmer(vectorInfected ~ week*trt + (1|cage),
                     data = acqDataVector, family = "binomial")
summary(infectedMod)


acqSummary <- acqDataCage %>% group_by(week, trt) %>% summarise(meancfu = mean(cagecfu),
                                                                secfu = sd(cagecfu)/sqrt(length(cagecfu[!is.na(cagecfu)])),
                                                                logMeancfu = mean(logCagecfu),
                                                                logSEcfu = sd(logCagecfu)/sqrt(length(logCagecfu[!is.na(logCagecfu)])),
                                                                meanPercInfected = mean(propVectorInfected)*100,
                                                                sePercInfected = sd(propVectorInfected)/sqrt(length(propVectorInfected[!is.na(propVectorInfected)]))*100)
acqSummary

#### Plotting
## Mean CFU
vectorxfplot <- ggplot(data=acqSummary, aes(x=week, y=meancfu)) +
  # geom_bar(position=position_dodge(), stat="identity", 
  #          aes(fill=trt)) +
  # geom_hline(aes(yintercept=50), linetype="dashed") +
  geom_line(aes(linetype=trt), size=1.25) +
  geom_point(aes(shape=trt), size=2.5) +
  geom_errorbar(aes(ymax=meancfu+secfu, ymin=meancfu-secfu), width=0.2) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = c(3,8,12)) + 
  scale_y_continuous(name = "Xylella CFU/vector") +
  # ylab("% insects on source plant") + 
  # ylim(c(0,100)) +
  # xlab("Weeks post inoculation") +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 
vectorxfplot
ggsave("results/figures/vector_xf_line_plot.jpg", plot = vectorxfplot,
       width = 7, height = 7, units = "in")


## Log mean CFU
logvectorxfplot <- ggplot(data=acqSummary, aes(x=week, y=logMeancfu)) +
  # geom_bar(position=position_dodge(), stat="identity", 
  #          aes(fill=trt)) +
  # geom_hline(aes(yintercept=50), linetype="dashed") +
  geom_line(aes(linetype=trt), size=1.25) +
  geom_point(aes(shape=trt), size=2.5) +
  geom_errorbar(aes(ymax=logMeancfu+logSEcfu, ymin=logMeancfu-logSEcfu), width=0.2) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = c(3,8,12)) + 
  scale_y_continuous(name = "Xylella CFU per vector (log10 transformed)") +
  # ylab("% insects on source plant") + 
  # ylim(c(0,100)) +
  # xlab("Weeks post inoculation") +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 
logvectorxfplot
ggsave("results/figures/vector_xf_line_plot_log.jpg", plot = logvectorxfplot,
       width = 7, height = 7, units = "in")


## Proportion of vectors infected
propInfectiousplot <- ggplot(acqSummary, aes(x=week, y=meanPercInfected)) +
  # geom_bar(position=position_dodge(), stat="identity", 
  #          aes(fill=trt)) +
  # geom_hline(aes(yintercept=50), linetype="dashed") +
  geom_line(aes(linetype=trt), size=1.25) +
  geom_point(aes(shape=trt), size=2.5) +
  geom_errorbar(aes(ymax=meanPercInfected+sePercInfected, ymin=meanPercInfected-sePercInfected), width=0.2) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = c(3,8,12)) + 
  scale_y_continuous(name = "Mean percent vectors positive for X. fastidiosa",
                     limits = c(0,100)) +
  # ylab("% insects on source plant") + 
  # ylim(c(0,100)) +
  # xlab("Weeks post inoculation") +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

propInfectiousplot
ggsave("results/figures/vector_prop_infectious_line_plot.jpg", plot = propInfectiousplot,
       width = 7, height = 7, units = "in")

#############################################################################################################
#### Transmission models
transdata <- readRDS("output/pdr1_transmission_preference_dataset.rds")
transdata$source.cfu.per.g <- as.numeric(transdata$source.cfu.per.g)
transdata$test.plant.infection <- as.integer(transdata$test.plant.infection)
str(transdata)

# Remove the second 12-week trials and NAs
transdata <- transdata %>% dplyr::filter(., week != 12.2, !is.na(pd_index) & !is.na(pd_index2))

#### Model selection on PD symptom index
pdMod1 <- glm(test.plant.infection ~ week*trt + log10(source.cfu.per.g+1) + p1 + p2 + mu1 + mu2 + pd_index, data = transdata, family = "binomial")
pdMod2 <- glm(test.plant.infection ~ week*trt + log10(source.cfu.per.g+1) + p1 + p2 + mu1 + mu2 + pd_index2, data = transdata, family = "binomial")
# I think quasibinomial distribution might be better but then it doesn't calculate an AIC value. Need to look into this.
AICtab(pdMod1, pdMod2, base = TRUE)
plot(pdMod2)
summary(pdMod2)
# PD indices are essentially the same; go with Arash's index

# Re-load data set to retain data points that are NA for pd indices
transdata <- readRDS("output/pdr1_transmission_preference_dataset.rds") %>% dplyr::filter(., week != 12.2)
transdata$source.cfu.per.g <- as.numeric(transdata$source.cfu.per.g)
transdata$test.plant.infection <- as.integer(transdata$test.plant.infection)


#### Model selection on preference rate parameters
FullModel <- glm(test.plant.infection ~ week*trt + log10(source.cfu.per.g+1) + p1 + p2 + mu1 + mu2 + pd_index, data = transdata, family = "binomial")
prefModChoice1 <- glm(test.plant.infection ~ week*trt + log10(source.cfu.per.g+1) + p1 + mu1 + pd_index, data = transdata, family = "binomial")
prefModp1 <- glm(test.plant.infection ~ week*trt + log10(source.cfu.per.g+1) + p1 + pd_index, data = transdata, family = "binomial")
prefModmu1 <- glm(test.plant.infection ~ week*trt + log10(source.cfu.per.g+1) + mu1 + pd_index, data = transdata, family = "binomial")
prefModnull <- glm(test.plant.infection ~ week*trt + log10(source.cfu.per.g+1) + pd_index, data = transdata, family = "binomial")

AICtab(FullModel, prefModChoice1, prefModp1, prefModmu1, prefModnull, base = TRUE)
# Including only mu1 model is best
summary(prefModmu1)


#### Model selection to include source xf pop and PD index
# start from prefModmu1
noSourceMod <- glm(test.plant.infection ~ week*trt*mu1*pd_index, data = transdata, family = "binomial")
noPDMod <- glm(test.plant.infection ~ week*trt*log10(source.cfu.per.g+1)*mu1, data = transdata, family = "binomial")
trtprefMod <- glm(test.plant.infection ~ week*trt*mu1, data = transdata, family = "binomial")
trtMod <- glm(test.plant.infection ~ week*trt, data = transdata, family = "binomial")
noIntrxnMod <- glm(test.plant.infection ~ week + trt, data = transdata, family = "binomial")

AICctab(noSourceMod, noPDMod, trtprefMod, trtMod, noIntrxnMod, base = TRUE)
# trtMod seems best
summary(noIntrxnMod)

# Contrast testing difference in transmission between genotypes at 12 weeks
transdata$week.trt <- with(transdata, factor(paste(week, trt, sep = ".")))
trtModMC <- glm(test.plant.infection ~ week.trt, data = transdata, family = "binomial")
transContrastTest <- glht(trtModMC, linfct = mcp(week.trt = "Tukey"))
summary(transContrastTest)




#### source xf pop
boxcox(source.cfu.per.g + 1 ~ week*trt, data = transdata, lambda = seq(-2, 2, by=0.5))

sourcexfMod <- lm(log10(source.cfu.per.g+1) ~ week*trt, data = transdata)
sourceNoInterxnMod <- lm(log10(source.cfu.per.g+1) ~ week + trt, data = transdata)
AICctab(sourcexfMod, sourceNoInterxnMod, base = TRUE)

plot(sourceNoInterxnMod)
summary(sourceNoInterxnMod)

# Multiple comparison testing difference in transmission between genotypes at 12 weeks
sourceModMC <- lm(log10(source.cfu.per.g+1) ~ week.trt, data = transdata)
transContrastTest <- glht(sourceModMC, linfct = mcp(week.trt = "Tukey"))
summary(transContrastTest)


## Proportion of source plants infected by week and trt
transdata$source.plant.infection <- ifelse(transdata$source.cfu.per.g == 0, 0, 1)
sourceProp <- transdata %>% group_by(., week, trt) %>% summarise(propInfected = sum(source.plant.infection, na.rm = TRUE)/sum(!is.na(source.plant.infection)))
sourceProp

sourcePropMod <- glm(source.plant.infection ~ week*trt, data = transdata, family = "binomial")
summary(sourcePropMod)
# Need to use Firth's correction
sourcePropMod <- logistf(source.plant.infection ~ week*trt, data = transdata)
summary(sourcePropMod)

#######################################################################################################
#### Plotting
transSummary <- transdata %>% group_by(., week, trt) %>% 
  summarise(percInfected = 100*sum(test.plant.infection, na.rm = TRUE)/sum(!is.na(test.plant.infection)),
            meanPD = mean(pd_index, na.rm = TRUE),
            sePD = sd(pd_index, na.rm = TRUE)/sqrt(sum(!is.na(pd_index))),
            meancfu = mean(log10(source.cfu.per.g+1)),
            secfu = sd(log10(source.cfu.per.g+1))/sqrt(sum(!is.na(source.cfu.per.g))))

# Transmission plot
transplot <- ggplot(data=transSummary, aes(x=week, y=percInfected)) +
  # geom_bar(position=position_dodge(), stat="identity", 
  #          aes(fill=trt)) +
  # geom_hline(aes(yintercept=50), linetype="dashed") +
  geom_line(aes(linetype=trt), size=1.25) +
  geom_point(aes(shape=trt), size=2.5) +
  #geom_errorbar(aes(ymax=meancfu+secfu, ymin=meancfu-secfu), width=0.2) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = c(3,8,12)) + 
  scale_y_continuous(name = "Percent test plants positive for X. fastidiosa",
                     limits = c(0,100)) +
  # ylab("% insects on source plant") + 
  # ylim(c(0,100)) +
  # xlab("Weeks post inoculation") +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

transplot
ggsave("results/figures/transmission_line_plot.jpg", plot = transplot,
       width = 7, height = 7, units = "in")


#### Symptoms plot
symptomPlot <- ggplot(data=transSummary, aes(x=week, y=meanPD)) +
  # geom_bar(position=position_dodge(), stat="identity", 
  #          aes(fill=trt)) +
  # geom_hline(aes(yintercept=50), linetype="dashed") +
  geom_line(aes(linetype=trt), size=1.25) +
  geom_point(aes(shape=trt), size=2.5) +
  geom_errorbar(aes(ymax=meanPD+sePD, ymin=meanPD-sePD), width=0.2) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = c(3,8,12)) + 
  scale_y_continuous(name = "Mean PD symptoms index") +
  # ylab("% insects on source plant") + 
  # ylim(c(0,100)) +
  # xlab("Weeks post inoculation") +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

symptomPlot
ggsave("results/figures/pd_symptom_line_plot.jpg", plot = symptomPlot,
       width = 7, height = 7, units = "in")


#### Xf pops in source plant plot
sourcexfplot <- ggplot(data=transSummary, aes(x=week, y=meancfu)) +
  # geom_bar(position=position_dodge(), stat="identity", 
  #          aes(fill=trt)) +
  # geom_hline(aes(yintercept=50), linetype="dashed") +
  geom_line(aes(linetype=trt), size=1.25) +
  geom_point(aes(shape=trt), size=2.5) +
  geom_errorbar(aes(ymax=meancfu+secfu, ymin=meancfu-secfu), width=0.2) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = c(3,8,12)) + 
  scale_y_continuous(name = "Xylella CFU/g plant tissue (log10)") +
  # ylab("% insects on source plant") + 
  # ylim(c(0,100)) +
  # xlab("Weeks post inoculation") +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

sourcexfplot
ggsave("results/figures/source_xf_line_plot_log.jpg", plot = sourcexfplot,
       width = 7, height = 7, units = "in")



#### Xf pops in source plant plot scatter plot
xfscatterplot <- ggplot(data=transdata, aes(x=week, y=(source.cfu.per.g/1000000))) +
  # geom_bar(position=position_dodge(), stat="identity", 
  #          aes(fill=trt)) +
  # geom_hline(aes(yintercept=50), linetype="dashed") +
  #geom_line(aes(linetype=trt), size=1.25) +
  geom_point(aes(shape=trt), size=2.5) +
  #geom_errorbar(aes(ymax=meancfu+secfu, ymin=meancfu-secfu), width=0.2) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = c(3,8,12)) + 
  scale_y_continuous(name = "Xylella CFU/g plant tissue (x 1,000,000)") +
  # ylab("% insects on source plant") + 
  # ylim(c(0,100)) +
  # xlab("Weeks post inoculation") +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

xfscatterplot
ggsave("results/figures/source_xf_pop_scatterplot.jpg", plot = xfscatterplot,
       width = 7, height = 7, units = "in")



##############################################################################################################
##############################################################################################################
#### 2017 data
##############################################################################################################

#### Importing and munging data
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



#####################################################################################################################################
#### Analysis of PD symptoms using ANCOVA
# pdMod1 includes week:genotype interaction, which tests for different slopes and intercepts
boxcox((PD_symptoms_index+1) ~ block + week*genotype, data = transdata, lambda = seq(-2, 2, by=0.5))
# Best transmformation is inverse sqrt; residuals don't look great, but better that with a quasipoisson GLM
pdMod1 <- lm(1/sqrt(PD_symptoms_index + 1) ~ block + week*genotype, data = transdata)
plot(pdMod1)
anova(pdMod1)
summary(pdMod1)


# Summarising and plotting
pdSummary <- transdata %>% group_by(week, genotype, trt) %>% 
  summarise(meanPD = mean(PD_symptoms_index, na.rm = TRUE), 
            sePD = sd(PD_symptoms_index, na.rm = TRUE)/sqrt(sum(!is.na(PD_symptoms_index))))


# PD symptoms plot
PDplot <- ggplot(data=pdSummary, aes(x=week, y=meanPD)) +
  geom_line(aes(linetype=genotype, colour = trt), size=1.25) +
  geom_point(aes(shape=genotype, colour = trt), size=3.5) +
  geom_errorbar(aes(ymax=meanPD+sePD, ymin=meanPD-sePD), width=0.2) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = c(2,5,8,14)) + 
  scale_y_continuous(name = "Mean PD symptom index",
                     limits = c(0,5)) +
  # ylab("% insects on source plant") + 
  # ylim(c(0,100)) +
  # xlab("Weeks post inoculation") +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

PDplot
ggsave("results/figures/2017_figures/pd_line_plot_2017.jpg", plot = PDplot,
       width = 7, height = 7, units = "in")




##############################################################################################################
#### Analysis of transmission data with binomial GLM

# How many trials have positive test plant and positive pre-screen plant?
transdata %>% dplyr::filter((grepl("PS plant positive", notes) & test_plant_infection == 1))
# Remove trials where test plant is positive and pre-screen plant is positive
transdata2 <- transdata %>% dplyr::filter(!(grepl("PS plant positive", notes) & test_plant_infection == 1))

transMod1 <- glm(test_plant_infection ~ block + week*genotype, family = "binomial", data = transdata2)
summary(transMod1)
transMod2 <- glm(test_plant_infection ~ block + week*trt, family = "binomial", data = transdata2)
transMod3 <- glm(test_plant_infection ~ block + week, family = "binomial", data = transdata2)
transMod4 <- glm(test_plant_infection ~ week*trt, family = "binomial", data = transdata2)
transMod5 <- glm(test_plant_infection ~ week + trt, family = "binomial", data = transdata2)
AICctab(transMod1, transMod2, transMod3, transMod4, transMod5, base = TRUE)
summary(transMod3)
# Note: best model doesn't include treatment or genotype, no difference there. And no trend with week; however, this is testing for a linear trend with week.

#### Summary by block
blockTrans <- transdata2 %>% group_by(block) %>% 
  summarise(percInfected = 100*(sum(test_plant_infection)/length(test_plant_infection)))

#### Summarize and plot data
transSummary <- transdata2 %>% group_by(week, genotype, trt) %>% 
  summarise(n = length(test_plant_infection),
            nInfected = sum(test_plant_infection),
            propInfected = nInfected/n)


# Transmission plot
transplot <- ggplot(data=transSummary, aes(x=week, y=propInfected)) +
  geom_line(aes(linetype=genotype, colour = trt), size=1.25) +
  geom_point(aes(shape=genotype, colour = trt), size=3.5) +
  #geom_errorbar(aes(ymax=meancfu+secfu, ymin=meancfu-secfu), width=0.2) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = c(2,5,8,14)) + 
  scale_y_continuous(name = "Proportion test plants positive for Xylella",
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

transplot

ggsave("results/figures/2017_figures/transmission_line_plot_2017.jpg", plot = transplot,
       width = 7, height = 7, units = "in")



##############################################################################################################
#### Fitting non-linear models to transmission data

# Load NLL functions for all models
source("R_functions/nll_functions_pdr1.R")


#### Fit a Holling Type IV model
#### Use model selection to assess if separate parameters are justified for S and R genotypes, and compare to linear model

# Setting up data
transSummarynl <- transdata2 %>% group_by(week, trt) %>% 
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

# #### Tests for Holling IV model  
# transR <- transSummarynl[transSummarynl$trt == "R",]
# 
# # Test NLL function
# with(transR, holling4NLL(a0, b0, c0, week, nInfected))
# lapply(1:nrow(transR), function(z) holling4NLL(a0, b0, c0, transR$week[z], transR$nInfected[z]))
# 
# holling4Rmod <- mle2(holling4NLL, data = list(week = transR$week, nInfected = transR$nInfected),
#                   start = list(a = a0, b = b0, c = c0),
#                   #optimizer = "optimx",
#                   control = list(maxit = 10000))
# 
# 
# #### Compare models that split and average over genotypes
# with(transSummarynl, holling4NLL(a0, b0, c0, week, nInfected))
# with(transSummarynl, holling4NLL.trt(a0, a0, b0, b0, c0, c0, week, nInfected))
# 
# holling4null <- mle2(holling4NLL, data = list(week = transSummarynl$week, nInfected = transSummarynl$nInfected),
#                      start = list(a = a0, b = b0, c = c0),
#                      optimizer = "optim", method = "Nelder-Mead",
#                      control = list(maxit = 10000))
# 
# holling4.trt <- mle2(holling4NLL.trt, data = list(week = transSummarynl$week, nInfected = transSummarynl$nInfected),
#                      start = list(a.r = a0, a.s = a0, b.r = b0, b.s = b0, c.r = c0, c.s = c0),
#                      # start = list(a = a0, b.r = b0, b.s = b0, c.r = c0, c.s = c0), # Initial parameters for when a is not split by trt
#                      optimizer = "optim", method = "Nelder-Mead",
#                      parameters = list(a~trt, b~trt, c~trt),
#                      control = list(maxit = 10000))
# 
# linearnull <- mle2(linearNLL, data = list(week = transSummarynl$week, nInfected = transSummarynl$nInfected),
#                    start = list(a = a.l0, b = b.l0),
#                    optimizer = "optim", method = "Nelder-Mead",
#                    control = list(maxit = 10000))
# 
# rickernull <- mle2(rickerNLL, data = list(week = transSummarynl$week, nInfected = transSummarynl$nInfected),
#                    start = list(a = a.rick0, b = b.rick0),
#                    optimizer = "optim", method = "Nelder-Mead",
#                    control = list(maxit = 10000))
# 
# ricker.trt <- mle2(rickerNLL.trt, data = list(week = transSummarynl$week, nInfected = transSummarynl$nInfected),
#                    start = list(a.r = a.rick0, a.s = a.rick0, b.r = b.rick0, b.s = b.rick0),
#                    optimizer = "optim", method = "Nelder-Mead",
#                    parameters = list(a~trt, b~trt, c~trt),
#                    control = list(maxit = 10000))
# 
# 
# # Compare linear parameters to those from a GLM
# glmMod <- glm(test_plant_infection ~ week, data = transdata2, family = "binomial")
# # The parameter estimates are way off. Not sure why. Not sure I've written the linearNLL correctly
# 
# 
# ICtab(holling4null, holling4.trt, 
#       linearnull, 
#       rickernull, ricker.trt,
#       type = "AICc", sort = TRUE, base = TRUE, nobs = 16)



##############################################################################################################
#### Fitting multiple non-linear models to Resistant and Susceptible trts separately

# Get data sets
transR <- transSummarynl[transSummarynl$trt == "R",]
transS <- transSummarynl[transSummarynl$trt == "S",]


#### Fitting resistant line data
transRresults <- optimizeTransModels(transR)
(modelSelectR <- transRresults[[2]])
opListR <- transRresults[[1]]

# Get model predictions for plotting
bestModR <- opListR$holling4Op # Holling4 model is the best
bestcoefR <- as.list(coef(bestModR))

newDatR <- data.frame(week = seq(2, 14, length = 101))
newDatR$propInfected <- with(newDatR, (bestcoefR$a*week^2)/(bestcoefR$b + bestcoefR$c*week + week^2))

#### Fitting susceptible line data
transSresults <- optimizeTransModels(transS)
(modelSelectS <- transSresults[[2]])
opListS <- transSresults[[1]]

# Get model predictions for plotting
bestModS <- opListS$rickerOp # Ricker model is the best
bestcoefS <- as.list(coef(bestModS))

newDatS <- data.frame(week = seq(2, 14, length = 101))
newDatS$propInfected <- with(newDatS, bestcoefS$a*week*exp(-bestcoefS$b*week))


#### Plotting results
# Plot S and R genotypes summarised together
# Use Ricker model for S and Holling 4 model for R
# Transmission plot
transplotNL <- ggplot(data=transSummarynl, aes(x=week, y=propInfected)) +
  #geom_line(aes(linetype=genotype, colour = trt), size=1.25) +
  geom_point(aes(shape = trt, colour = trt), size=3.5) +
  #geom_errorbar(aes(ymax=meancfu+secfu, ymin=meancfu-secfu), width=0.2) +
  # Defining the colors for the lines based on the default ggplot2 colors and ggplot_build()$data
  geom_smooth(data = newDatR, method = "loess", colour = "#F8766D", se = FALSE) +
  geom_smooth(data = newDatS, method = "loess", colour = "#00BFC4", se = FALSE) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = c(2,5,8,14)) + 
  scale_y_continuous(name = "Proportion test plants positive for Xylella",
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

transplotNL

ggsave("results/figures/2017_figures/transmission_non-linear_plot_2017.jpg", plot = transplotNL,
       width = 7, height = 7, units = "in")

