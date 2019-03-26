#### ANALYSIS OF PDR1 TRANSMISSION DATA

rm(list = ls())
# Load packages
my.packages <- c("tidyr", "dplyr", "data.table", "openxlsx", "MASS", "googlesheets", "caret",
                 "lme4", "ggplot2", "bbmle", "multcomp", "DHARMa", "glmnet", "glmnetUtils")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/factor2numeric.R")


######################################################################################################################
######################################################################################################################
#### Analysis of 2016 data
######################################################################################################################

#### Load combined and filtered transmission/preference data set
# Remove second 12-week set of trials
transdata <- readRDS("output/pdr1_transmission_preference_dataset.rds") %>% dplyr::filter(., week != 12.2)
transdata$source.cfu.per.g <- as.numeric(transdata$source.cfu.per.g)
transdata$test.plant.infection <- as.integer(transdata$test.plant.infection)
transdata$trt <- factor(transdata$trt)
str(transdata)

######################################################################################################################
#### Analysis of PD symptom data
# # Linear model with transformed symptom data
# boxcox(pd_index+1 ~ week*trt*log10(source.cfu.per.g+1), data = transdata)
# sympLM <- lm(log(pd_index+1) ~ week*trt*log10(source.cfu.per.g+1), data = transdata)
# plot(sympLM)
# summary(sympLM)
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


#############################################################################################################
#### source xf pop
## Data look overdispesed: try negative binomial GLM (MASS package)
sourcexfNB <- glm.nb(source.cfu.per.g ~ week*trt, data = transdata)
summary(sourcexfNB)
## Results are way off, don't make sense. Looks like model fitting didn't work
## Try quasi-poisson
sourcexfQP <- glm(source.cfu.per.g ~ week*trt, data = transdata, family = "quasipoisson")
plot(sourcexfQP)
summary(sourcexfQP)
## Even though the error variance is highly skewed, the results make sense


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



######################################################################################################################
#### Analyzing acquisition data at the per-vector level, with cage as a random effect
acqDataVector <- readRDS("output/pdr1_2016_vector_acquisition_dataset.rds")
str(acqDataVector)

# Analysis of CFU data
acqMod <- glmer(vectorcfu ~ week*trt + (1|week.cage), data = acqDataVector, family = "poisson")
summary(acqMod)

# Analysis of infection status
infectedMod <- glmer(vectorInfectious ~ week*trt + (1|week.cage), data = acqDataVector, family = "binomial")
plot(simulateResiduals(infectedMod))
summary(infectedMod)


#### Plotting
## Use transdata data set for plotting, b/c data are at the cage level
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


##############################################################################################################
#### Non-linear models of transmission



##############################################################################################################
##############################################################################################################
#### 2017 data
##############################################################################################################

transVCPdata <- readRDS("output/complete_2017_transmission-preference_dataset.rds")
str(transVCPdata)
with(transVCPdata, table(week, genotype))


##############################################################################################################
#### Analysis of PD symptoms using ANCOVA
# pdMod1 includes week:genotype interaction, which tests for different slopes and intercepts
boxcox((PD_symptoms_index+1) ~ block + week*genotype, data = transVCPdata, lambda = seq(-2, 2, by=0.5))
# Best transmformation is inverse sqrt; residuals don't look great, but better that with a quasipoisson GLM
pdMod1 <- lm(1/sqrt(PD_symptoms_index + 1) ~ block + week*genotype, data = transVCPdata)
plot(pdMod1)
anova(pdMod1)
summary(pdMod1)
## Analysis of PD symptoms using Partial Odds Logistic Regression
transVCPdata$PD_symptoms_index <- factor(transVCPdata$PD_symptoms_index, ordered = TRUE, levels = c("0", "1", "2", "3", "4", "5"))
olrMod <- polr(PD_symptoms_index ~ block + week*genotype, data = transVCPdata, Hess = TRUE, method = "logistic")
summary(olrMod)
# Calculate p values from t statistic
olrCoefs <- coef(summary(olrMod))
p <- pnorm(abs(olrCoefs[, "t value"]), lower.tail = FALSE)*2
(olrCoefs <- cbind(olrCoefs, "p-value" = p))
# Results: quasi-Poisson model, transformed LM model, and POLR model all give same result -> only week is significant positive 


##############################################################################################################
#### Analysis of Xf pops in source plants
## Linear model with transformation
boxcox(xfpop + 1 ~ block + week*genotype, data = transVCPdata, lambda = seq(-2, 2, by=0.5))
# Best transformation is either sqrt or log; go with square root because the residuals look better
sourcepopMod1 <- lm(sqrt(xfpop) ~ block + week*genotype, data = transVCPdata)
plot(simulateResiduals(sourcepopMod1))
## Residuals don't look very good
summary(sourcepopMod1)

## Generalized linear model with quasipoisson distribution
poispopMod1 <- glm(xfpop ~ block + week*genotype, data = transVCPdata, family = "quasipoisson")
plot(poispopMod1)
# Residuals look about the same as lm() with sqrt()
summary(poispopMod1)

## Negative binomial GLM
sourcexfNB <- glm.nb(xfpop ~ block + week*genotype, data = transVCPdata)
summary(sourcexfNB)
## Results are qualitatively similar to quasipoisson but returns a warning. Quasipoisson also works better for 2016 data; go with that.

#### NOTE ON ALTERNATIVE MODEL: Removing false negatives and all negatives had negligible effects on results.
#### These models had reduced significance of week main effect but were otherwise unchanged.
#### Log10 transformation works better if all negatives are removed


##############################################################################################################
#### Plotting PD symptoms, Xf pops, transmission, and acquisition from 2017

#### Plotting PD symptoms
## Summary
pdSummary <- transVCPdata %>% group_by(week, genotype, trt) %>% 
  summarise(meanPD = mean(PD_symptoms_index, na.rm = TRUE), 
            sePD = sd(PD_symptoms_index, na.rm = TRUE)/sqrt(sum(!is.na(PD_symptoms_index))))

## PD symptoms plot
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


#### Plotting source plant xf populations
## Summary
sourceSummary <- sourcedata2 %>% mutate(logxfpop = log10(xfpop+1), sqrtxfpop = sqrt(xfpop)) %>% 
  group_by(week, genotype, trt) %>% 
  summarise_at(c("logxfpop", "sqrtxfpop"), funs(mean = mean(.), n = sum(!is.na(.)), se = sd(.)/sqrt(sum(!is.na(.)))))

## Xf pops in source plant plot
sourcexfplot <- ggplot(data=sourceSummary, aes(x=week, y=sqrtxfpop_mean)) +
  geom_line(aes(linetype=genotype, colour = trt), size=1.25) +
  geom_point(aes(shape=genotype, colour = trt), size=3.5) +
  geom_errorbar(aes(ymax=sqrtxfpop_mean+sqrtxfpop_se, ymin=sqrtxfpop_mean-sqrtxfpop_se), width=0.2) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = c(2,5,8,14)) + 
  scale_y_continuous(name = "Xylella populations in source plant (square root)") +
  # limits = c(0,10)) +
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

ggsave("results/figures/2017_figures/source_xf_sqrt_line_plot_2017.jpg", plot = sourcexfplot,
       width = 7, height = 7, units = "in")




#### Plotting raw transmission
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

