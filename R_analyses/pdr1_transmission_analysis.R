#### ANALYSIS OF PDR1 TRANSMISSION DATA

rm(list = ls())
# Load packages
my.packages <- c("tidyr", "dplyr", "data.table", "openxlsx", "ggplot2",
                 "MASS", "logistf", "multcomp", "bbmle", "lme4", "googlesheets")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/factor2numeric.R")


##############################################################################################################
##############################################################################################################
#### 2016 data
##############################################################################################################

#### Munging data

# Preference data
prefdata <- read.xlsx("data/2016_data/pdr1_preference_data.xlsx", sheet = "data")
str(prefdata)


# Separate leaf and symptom data from preference count data
leafdata <- prefdata[,c("week", "trt", "rep", "n_leaves_test", "n_leaves_source_start", 
                        "n_leaves_source_end", "n_pd_leaves", "n_ms_petioles", "pd_index")] %>% dplyr::filter(., !is.na(n_pd_leaves))
leafdata$n_leaves_test <- factor2numeric(leafdata$n_leaves_test)
leafdata$n_leaves_source_start <- factor2numeric(leafdata$n_leaves_source_start)
leafdata$n_leaves_source_end <- factor2numeric(leafdata$n_leaves_source_end)
# Some trials I only counted leaves at beginning, other trials only at end; combine total leaf data
leafdata$n_leaves_source <- leafdata %>% with(., cbind(n_leaves_source_start, n_leaves_source_end)) %>% rowMeans(., na.rm = TRUE)
leafdata$n_leaves_source[is.nan(leafdata$n_leaves_source)] <- NA
str(leafdata)


################################################################################################################
#### Constructing PD indices

#### PD index of Rashed et al. 2013:
# 0 = asymptomatic
# 1 = 1-2 scorched leaves
# 2 = 3-4 scorched leaves
# 3 = all leaves scorched; a few matchstick petioles
# 4 = all leaves heavily scorched and many matchstick petioles
# 5 = only a few leaves at the end of the cane

table(leafdata$n_pd_leaves, leafdata$n_ms_petioles)

#### Alternative PD index: 
prop_pd_leaves <- leafdata %>% with(., ifelse(n_leaves_source == 0, 1, n_pd_leaves/n_leaves_source))
leafdata$pd_index2 <- leafdata %>% with(., prop_pd_leaves + n_ms_petioles)
plot(x = leafdata$pd_index, y = leafdata$pd_index2)


####################################################################################################################
#### Import and combine culturing data, leaf data, and preference data
# Import culturing data
culturedata <- read.xlsx("data/2016_data/pdr1_culturing_data.xlsx", sheet = "Infection data")
# Remove repeat culturing data
culturedata <- culturedata %>% dplyr::filter(., notes != "repeat culture" | is.na(notes))
str(culturedata)
culturedata

# Import preference data: estimates of rate parameters from CM model for each cage
# choice 1 = source plants, choice 2 = test plants
# Ignoring variance around parameters.
# TO DO: check on convergence issues
paramDataCage <- readRDS("output/CMM_rate_parameters_per_cage.rds")
# Reshape paramDataCage
paramDataCage <- paramDataCage[,c("parameter", "estimate", "week.cage")] %>% spread(., key = parameter, value = estimate)

# Merge culturing and leaf data on week-trt-rep combination
transdata <- inner_join(culturedata, leafdata, by = c("week", "trt", "rep")) 

# Merge with preference data
transdata$week.cage <- transdata %>% with(., paste(week, trt, rep, sep=""))
transdata <- left_join(transdata, paramDataCage, by = "week.cage")
str(transdata)

saveRDS(transdata, file = "output/pdr1_transmission_preference_dataset.rds")


#############################################################################################################
#### Import vector acquisition data and combine with transmission-preference dataset
acqData <- readRDS("output/pdr1_2016_vector_cfu_from_qpcr.rds")
# Remove qPCR serial dilution and NTC rows
acqData <- acqData %>% dplyr::filter(!is.na(insect_code))
str(acqData)
# Remove one crazy outlier!
acqData <- acqData %>% dplyr::filter(!(tube_code == 117 & wellNumber == 60))
# Simplify data set
acqDataVector <- acqData %>% dplyr::select(week, trt, rep, tube_code, cfu) %>% group_by(week, trt, rep, tube_code) %>% summarise(vectorcfu = mean(cfu))
acqDataVector$week <- as.numeric(acqDataVector$week)
acqDataVector$rep <- as.numeric(acqDataVector$rep)
str(acqDataVector)


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

#### Constructing vector infection index
# In each trial, I have Xf pops for each of the vectors. I could include in transmission model:
# total Xf population among all vectors
# proportion of vectors infected (b/c some cages have fewer than 8 vectors)
# an evenness index of Xf pops among vectors
# Matt's transmission parameters paper might have something to say about this
# multiple possibilities can be evaluated using AIC

## Average duplicates for each sample
acqDataCage <- acqDataVector %>% group_by(week, trt, rep) %>% summarise(cagecfu = mean(vectorcfu),
                                                                        sdcfu = sd(vectorcfu),
                                                                        logCagecfu = mean(log10(vectorcfu + 1)),
                                                                        propVectorInfected = sum(vectorcfu > 0, na.rm = TRUE)/length(vectorcfu[!is.na(vectorcfu)]))
printTibble(acqDataCage)

## Merge with acquisition data at cage level with transmission-preference data set
acqDataCage$week.cage <- with(acqDataCage, paste(week, trt, rep, sep=""))

transdata <- readRDS("output/pdr1_transmission_preference_dataset.rds")

transdata <- left_join(transdata, acqDataCage, by = c("week.cage", "week", "trt", "rep"))
transdata
saveRDS(transdata, file = "output/pdr1_transmission_preference_dataset.rds")


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


#### symptom data
sympMod <- glm(pd_index ~ week*trt*log10(source.cfu.per.g+1), data = transdata, family = "poisson")
sympTrtMod <- glm(pd_index ~ week*trt, data = transdata, family = "poisson")
sympNoInterxnMod <- glm(pd_index ~ week + trt, data = transdata, family = "poisson")
AIC(sympMod, sympTrtMod, sympNoInterxnMod)
# sympTrtMod (with interaction) is best and quasipoisson distribution is best, AIC won't work with quasipoisson
sympTrtMod <- glm(pd_index ~ week*trt, data = transdata, family = "quasipoisson")
plot(sympTrtMod)
summary(sympTrtMod)


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
str(transdata)
summary(transdata)

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



#####################################################################################################################################
#### Analysis of PD symptoms
pdMod1 <- glm(PD_symptoms_index ~ block + week*genotype, family = "poisson", data = transdata)
summary(pdMod1)
pdMod2 <- glm(PD_symptoms_index ~ block + week*trt, family = "poisson", data = transdata)
AICctab(pdMod1, pdMod2, base = TRUE)
summary(pdMod2)

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
# How do I test for a non-linear trend with week?

#### Summary by block
blockTrans <- transdata2 %>% group_by(block) %>% 
  summarise(percInfected = 100*(sum(test_plant_infection)/length(test_plant_infection)))

#### Summarize and plot data
transSummary <- transdata2 %>% group_by(week, genotype, trt) %>% 
  summarise(n = length(test_plant_infection),
            nInfected = sum(test_plant_infection),
            percInfected = 100*(nInfected/n))


# Transmission plot
transplot <- ggplot(data=transSummary, aes(x=week, y=percInfected)) +
  geom_line(aes(linetype=genotype, colour = trt), size=1.25) +
  geom_point(aes(shape=genotype, colour = trt), size=3.5) +
  #geom_errorbar(aes(ymax=meancfu+secfu, ymin=meancfu-secfu), width=0.2) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = c(2,5,8,14)) + 
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


#### Tests for Holling IV model  
transR <- transSummarynl[transSummarynl$trt == "R",]

# Test NLL function
with(transR, holling4NLL(a0, b0, c0, week, nInfected))
lapply(1:nrow(transR), function(z) holling4NLL(a0, b0, c0, transR$week[z], transR$nInfected[z]))

holling4Rmod <- mle2(holling4NLL, data = list(week = transR$week, nInfected = transR$nInfected),
                  start = list(a = a0, b = b0, c = c0),
                  #optimizer = "optimx",
                  control = list(maxit = 10000))


#### Compare models that split and average over genotypes
with(transSummarynl, holling4NLL(a0, b0, c0, week, nInfected))
with(transSummarynl, holling4NLL.trt(a0, a0, b0, b0, c0, c0, week, nInfected))

holling4null <- mle2(holling4NLL, data = list(week = transSummarynl$week, nInfected = transSummarynl$nInfected),
                     start = list(a = a0, b = b0, c = c0),
                     optimizer = "optim", method = "Nelder-Mead",
                     control = list(maxit = 10000))

holling4.trt <- mle2(holling4NLL.trt, data = list(week = transSummarynl$week, nInfected = transSummarynl$nInfected),
                     start = list(a.r = a0, a.s = a0, b.r = b0, b.s = b0, c.r = c0, c.s = c0),
                     # start = list(a = a0, b.r = b0, b.s = b0, c.r = c0, c.s = c0), # Initial parameters for when a is not split by trt
                     optimizer = "optim", method = "Nelder-Mead",
                     parameters = list(a~trt, b~trt, c~trt),
                     control = list(maxit = 10000))

linearnull <- mle2(linearNLL, data = list(week = transSummarynl$week, nInfected = transSummarynl$nInfected),
                   start = list(a = a.l0, b = b.l0),
                   optimizer = "optim", method = "Nelder-Mead",
                   control = list(maxit = 10000))

rickernull <- mle2(rickerNLL, data = list(week = transSummarynl$week, nInfected = transSummarynl$nInfected),
                   start = list(a = a.rick0, b = b.rick0),
                   optimizer = "optim", method = "Nelder-Mead",
                   control = list(maxit = 10000))

ricker.trt <- mle2(rickerNLL.trt, data = list(week = transSummarynl$week, nInfected = transSummarynl$nInfected),
                   start = list(a.r = a.rick0, a.s = a.rick0, b.r = b.rick0, b.s = b.rick0),
                   optimizer = "optim", method = "Nelder-Mead",
                   parameters = list(a~trt, b~trt, c~trt),
                   control = list(maxit = 10000))


# Compare linear parameters to those from a GLM
glmMod <- glm(test_plant_infection ~ week, data = transdata2, family = "binomial")
# The parameter estimates are way off. Not sure why. Not sure I've written the linearNLL correctly


ICtab(holling4null, holling4.trt, 
      linearnull, 
      rickernull, ricker.trt,
      type = "AICc", sort = TRUE, base = TRUE, nobs = 16)



##############################################################################################################
#### Fitting multiple models to Resistant and Susceiptle trts separately

# Get data sets
transR <- transSummarynl[transSummarynl$trt == "R",]
transS <- transSummarynl[transSummarynl$trt == "S",]


#### Fitting resistant line data
transRresults <- optimizeTransModels(transR)
modelSelectR <- transRresults[[2]]
opListR <- transRresults[[1]]

# Get model predictions for plotting
bestModR <- opListR$holling4Op # Holling4 model is the best
bestcoefR <- as.list(coef(bestModR))

newDatR <- data.frame(week = seq(2, 14, length = 101))
newDatR$propInfected <- with(newDatR, (bestcoefR$a*week^2)/(bestcoefR$b + bestcoefR$c*week + week^2))

#### Fitting susceptible line data
transSresults <- optimizeTransModels(transS)
modelSelectS <- transSresults[[2]]
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
  scale_y_continuous(name = "Proportion test plants positive for X. fastidiosa",
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

