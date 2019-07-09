#### SYNTHETIC ANALYSES OF COMBINED TRANSMISSION-RELATED DATA AND CHEMICAL DATA

rm(list = ls())
# Load packages
my.packages <- c("tidyr", "dplyr", "data.table", "ggplot2", "DHARMa",
                 "bbmle", "glmnetUtils", "caret",
                 "modeest", "cowplot", "broom", "car", "ggbiplot")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/factor2numeric.R")

#### For CMM parameters
## choice1 = source (Xylella infected) plant
## choice2 = test plant


#############################################################################################################
#############################################################################################################
#### 2016 Combined analysis of transmission data
#############################################################################################################

#### Load combined and filtered transmission/preference data set
# Remove second 12-week set of trials
transdata <- readRDS("output/pdr1_transmission_preference_dataset.rds") %>% dplyr::filter(., week != 12.2)
transdata$source.cfu.per.g <- as.numeric(transdata$source.cfu.per.g)
transdata$test.plant.infection <- as.integer(transdata$test.plant.infection)
transdata$trt <- factor(transdata$trt)
str(transdata)


#### Full linear model
FullModel <- glm(test.plant.infection ~ week*trt + propVectorInfected + source.cfu.per.g + p1 + p2 + mu1 + mu2 + pd_index, data = transdata, family = "binomial")
simResFullModel <- simulateResiduals(FullModel, n = 1000)
plot(simResFullModel)
diagnosticsFullModel <- testResiduals(simResFullModel)
## Model diagnostics suggest a good model
summary(FullModel)

#### Look at relationships among variables
transdata %>% dplyr::select(test.plant.infection, propVectorInfected, source.cfu.per.g, pd_index, p1, p2, mu1, mu2) %>% pairs()

transdata %>% group_by(trt) %>% summarise(propInfected = sum(test.plant.infection, na.rm = TRUE)/length(test.plant.infection[!is.na(test.plant.infection)]))

#### Transmission analysis using elastic net 
## Define lambda values
lambdas <- 2^seq(-1, -8, length = 20)
#### Using the glmnetUtils package
enetTransdata <- transdata %>% mutate(log.source.cfu = log10(source.cfu.per.g+1),
                                      # Trt code: R = 0, S = 1; so a positive relationship between trt and infection status indicates greater trans from Susceptible
                                      trtNumeric = as.numeric(trt) - 1,
                                      test.plant.infection = factor(ifelse(test.plant.infection == 1, "infected", "non_infected"))) %>%
  dplyr::select(test.plant.infection, trtNumeric, propVectorInfected, log.source.cfu, pd_index, mu1, mu2, p1, p2) %>%
  dplyr::filter(complete.cases(.))
## Switch factor level order to get correct sign on coefficient estimates
enetTransdata$test.plant.infection <- factor(enetTransdata$test.plant.infection, levels(enetTransdata$test.plant.infection)[c(2,1)])
## "non-infected" should be the first level
levels(enetTransdata$test.plant.infection)

str(enetTransdata)





#### Transmission analysis using elastic net 
## Define lambda values
lambdas <- 2^seq(0, -8, length = 20)

#### Cross-validation of alpha and lambda
enetTransCV <- cva.glmnet(test.plant.infection ~ ., data = enetTransdata, family = "binomial",
                          nfolds = 5, lambda = lambdas)
enetTransCV
#plot(enetTransCV)
minlossplot(enetTransCV)
## best alpha results are all over the place
## need to check with caret package is using elastic net
# ## Plot just the cv.glmnet output for alpha = 0, which is the first in the list of cv.glmnet objects
# plot(enetTransCV$modlist[[1]])
# enetTransCV$modlist[[1]]
# bestLambda <- enetTransCV$modlist[[1]]$lambda.min
# ridgeTrans <- glmnetUtils::glmnet(test.plant.infection ~ ., data = enetTransdata, family = "binomial", alpha = 0)
# enetResults1 <- coef(ridgeTrans, s = bestLambda)


#### Transmission analysis using LASSO
lassoTransCV <- cv.glmnet(test.plant.infection ~ ., data = enetTransdata, family = "binomial",
                          alpha = 1, nfolds = 5, lambda = lambdas)
plot(lassoTransCV)
coef(lassoTransCV, s = lassoTransCV$lambda.min)


#### Transmission analysis using elastic net and the caret package
## Glmnet wants the data to be matrices, not data frames.
x_train <- enetTransdata %>% dplyr::select(-test.plant.infection) %>% as.matrix()
## For caret, the y variable must be a factor
## and I need to re-order so that "non-infected" is first, to get the coefficient signs correct
y_train <- enetTransdata[,"test.plant.infection"] %>% as.matrix() %>% factor()
y_train <- factor(y_train, levels(y_train)[c(2,1)])

## Set up trainControl
train_control = trainControl(method = "cv",
                             number = 5, 
                             #selectionFunction = "oneSE",
                             returnResamp = "all",
                             classProbs = TRUE, 
                             summaryFunction = twoClassSummary,
                             savePredictions = "final")

#### Create a custom tuning grid.
## Find max lambda for cross validation
## Formula found here: https://stats.stackexchange.com/questions/144994/range-of-lambda-in-elastic-net-regression
## Need to center and scale predictors first
xenet <- enetTransdata %>% dplyr::select(-test.plant.infection) %>% scale(., center = TRUE, scale = TRUE)
yenet <- ifelse(enetTransdata$test.plant.infection == "infected", 1, 0)
lambdaMax <- apply(xenet, 2, function(x) sum(yenet*x)) %>% max()
## set upper limit for lamba range based on lambdaMax and alpha
## from this website: https://stats.stackexchange.com/questions/144994/range-of-lambda-in-elastic-net-regression
alphas = seq(0, 1, by = 0.1)
lambdaMaxAdj <- 1/(1-alphas)*lambdaMax
lambdaMaxAdj[lambdaMaxAdj == Inf] <- max(lambdaMaxAdj[lambdaMaxAdj != Inf])
enetTuneList <- vector("list", length(alphas))
for(i in 1:length(alphas)){
  enetTuneList[[i]] <- expand.grid(alpha = alphas[i], 
                                   lambda = seq(0, lambdaMaxAdj[i], length.out = 20))
}
enet_grid <- enetTuneList %>% rbindlist() %>% as.data.frame()

#### Run cross-validation once
enet = train(x_train, y_train, method = "glmnet",
             metric = "ROC",
             preProcess = c("center", "scale"),
             tuneGrid = enet_grid,
             trControl = train_control)
print(enet)
plot(enet)
enet$bestTune
(enetResults2 <- coef(enet$finalModel, s = enet$bestTune$lambda))

enet3 <- glmnetUtils::glmnet(test.plant.infection ~ ., data = enetTransdata, family = "binomial",
                             alpha = enet$bestTune$alpha, lambda = lambdas)
(enetResults3 <- coef(enet3, s = enet$bestTune$lambda))
#### The cva.glmnet and caret results largely match up
## Results are highly variable across CV runs. Need to run CV multiple times and average coefficient results

#### Run caret::train() through for loop
## Could do this with "repeatedcv" method in trainControl, but don't know how to extract variance estimate from repeats
## Run it in a for loop instead

## Number of runs
nrun <- 500

## Empty vectors for coefficients and best tune parameters
enetResultsList <- enetBestTuneList <- vector("list", nrun)

for(i in 1:nrun){
  enetTrain = train(x_train, y_train, method = "glmnet",
                    metric = "ROC",
                    preProcess = c("center", "scale"),
                    tuneGrid = enet_grid,
                    trControl = train_control)
  enetBestTuneList[[i]] <- enetTrain$bestTune
  enet <- glmnetUtils::glmnet(test.plant.infection ~ ., data = enetTransdata, family = "binomial",
                              alpha = enetTrain$bestTune$alpha, lambda = lambdas)
  enetResultsList[[i]] <- as.data.frame(t(as.matrix(coef(enet, s = enetTrain$bestTune$lambda))))
}

## Combine and summarize best tuning parameters
enetBestTune <- enetBestTuneList %>% rbindlist() %>% as.data.frame()
hist(enetBestTune$alpha)
hist(enetBestTune$lambda)
## Get mode of alpha and lambda
(modeTune16 <- apply(enetBestTune, 2, mfv))

#### Combine runs and summarize coefficient estimates
enetResultsData <- enetResultsList %>% rbindlist() %>% as.data.frame()
enetSummary16 <- data.frame(meancoef = apply(enetResultsData, 2, mean),
                          mediancoef = apply(enetResultsData, 2, median),
                          sdcoef = apply(enetResultsData, 2, sd))

saveRDS(list(modeTune16, enetSummary16), file = "results/multi-run_CV_elastic_net_transmission_results_2016.rds")


#### Read in 2016 elastic net results
enet16list <- readRDS("results/multi-run_CV_elastic_net_transmission_results_2016.rds")
enetSummary16 <- enet16list[[2]]

#### Plot Elastic Net results
## Vector of clear names for covariates
enetSummary16$niceNames <- factor(c("Intercept", "Resistance trait", "Proportion vectors infectious",
                                    "Xylella population size in infected plant", "PD disese severity", 
                                    "Leaving rate from infected plant", "Leaving rate from test plant",
                                    "Attraction rate from infected plant", "Attraction rate from test plant"))
levels(enetSummary16$niceNames)
## Reorder factor levels according to coefficient estimate
#enetSummary16$niceNames <- factor(enetSummary16$niceNames, levels = enetSummary16$niceNames[order(enetSummary16$meancoef, decreasing = FALSE)])
## Reorder factor levels according to original order
enetSummary16$niceNames <- with(enetSummary16, factor(niceNames, levels(niceNames)[c(3,8,7,9,6,1,2,4,5)]))
levels(enetSummary16$niceNames)

trans16enetPlot <- ggplot(enetSummary16, aes(y = niceNames, x = meancoef)) +
  geom_errorbarh(aes(xmin = meancoef-sdcoef, xmax = meancoef+sdcoef), colour = "black", height = 0.2) +
  geom_point(size = 3) +
  geom_vline(linetype = "longdash", xintercept = 0) +
  xlab("Coefficient estimate") + ylab("") + 
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

trans16enetPlot

ggsave(filename = "results/figures/2016_figures/transmission_2016_elastic_net_coefficients_plot.tiff", 
       plot = trans16enetPlot,
       width = 12, height = 7, units = "in")


######################################################################################################
######################################################################################################
#### 2017 Combined transmission analysis
######################################################################################################

#### Load transmission, acquisition, culturing, preference (and phenolic, when it's ready) data set
transVCPdata <- readRDS("output/complete_2017_transmission-preference_dataset.rds")
str(transVCPdata)
summary(transVCPdata)

#### For CMM parameters
## choice1 = source (Xylella infected) plant
## choice2 = test plant

#### Remove mu2 outlier
transVCPdata <- transVCPdata %>% dplyr::filter(mu2 < 7)

#### Full linear model
FullModel <- glm(test_plant_infection ~ week*genotype + PD_symptoms_index + propInfectious + log10(xfpop+1) + mu1 + mu2 + p1 + p2,
                 data = transVCPdata, family = "binomial")
plot(simulateResiduals(FullModel))
summary(FullModel)

#### Look at relationships among variables
transVCPdata %>% dplyr::select(-week, -block, -genotype, -trt, -rep, -notes, -nbugs, -totalInfectious, -plantID, -Rep2) %>% pairs()
## Look at just p1 vs p2
tiff("output/figures/attraction_rates_2017_scatterplot.tiff")
  with(transVCPdata, plot(x=p1, y=p2, xlab = "Attraction rate to source plant", ylab = "Attraction rate to test plant"))
dev.off()


transVCPdata %>% group_by(trt) %>% summarise(propInfected = sum(test_plant_infection)/length(test_plant_infection))

######################################################################################################
#### Transmission analysis using elastic net 
## Define lambda values
lambdas <- 2^seq(-1, -8, length = 20)

#### Using the glmnetUtils package
enetTransdata <- transVCPdata %>% mutate(log.xfpop = log10(xfpop+1),
                                         # Trt code: R = 0, S = 1; so a positive relationship between trt and infection status indicates greater trans from Susceptible
                                         trtNumeric = as.numeric(trt) - 1, 
                                         test_plant_infection = factor(ifelse(test_plant_infection == 1, "infected", "non_infected"))) %>%
  dplyr::select(test_plant_infection, trtNumeric, PD_symptoms_index, propInfectious, log.xfpop, mu1, mu2, p1, p2) %>%
  dplyr::filter(complete.cases(.))
# Switch factor level order to get correct sign on coefficient estimates
enetTransdata$test_plant_infection <- factor(enetTransdata$test_plant_infection, levels(enetTransdata$test_plant_infection)[c(2,1)])
str(enetTransdata)

## Set up genotype factor levels
glev <- list(genotype = levels(enetTransdata$genotype))

#### Cross-validation of alpha and lambda
enetTransCV <- cva.glmnet(test_plant_infection ~ ., data = enetTransdata, family = "binomial",
                          nfolds = 5, lambda = lambdas, xlev = glev)
enetTransCV
#plot(enetTransCV)
minlossplot(enetTransCV)
## best alpha results are all over the place
## need to check with caret package is using elastic net
# ## Plot just the cv.glmnet output for alpha = 0, which is the first in the list of cv.glmnet objects
# plot(enetTransCV$modlist[[1]])
# enetTransCV$modlist[[1]]
# bestLambda <- enetTransCV$modlist[[1]]$lambda.min
# ridgeTrans <- glmnetUtils::glmnet(test.plant.infection ~ ., data = enetTransdata, family = "binomial", alpha = 0)
# enetResults1 <- coef(ridgeTrans, s = bestLambda)


#### Transmission analysis using LASSO
## Define lambda values
lambdas <- 2^seq(-1, -8, length = 20)

lassoTransCV <- cv.glmnet(test_plant_infection ~ ., data = enetTransdata, family = "binomial",
                          alpha = 1, nfolds = 5, lambda = lambdas, xlev = glev)
plot(lassoTransCV)
coef(lassoTransCV, s = lassoTransCV$lambda.min)


#####################################################################################################################
#### Transmission analysis using elastic net and the caret package
## Glmnet wants the data to be matrices, not data frames.
x_train <- enetTransdata %>% dplyr::select(-test_plant_infection) %>% as.matrix()
## For caret, the y variable must be a factor
## and I need to re-order so that "non-infected" is first, to get the coefficient signs correct
y_train <- enetTransdata[,"test_plant_infection"] %>% as.matrix() %>% factor()
y_train <- factor(y_train, levels(y_train)[c(2,1)])

## Set up trainControl
train_control = trainControl(method = "cv",
                             number = 5, returnResamp = "all",
                             classProbs = TRUE, summaryFunction = twoClassSummary,
                             savePredictions = "final")

#### Create a custom tuning grid.
## Find max lambda for cross validation
## Formula found here: https://stats.stackexchange.com/questions/144994/range-of-lambda-in-elastic-net-regression
## Need to center and scale predictors first
xenet <- enetTransdata %>% dplyr::select(-test_plant_infection) %>% scale(., center = TRUE, scale = TRUE)
yenet <- ifelse(enetTransdata$test_plant_infection == "infected", 1, 0)
lambdaMax <- apply(xenet, 2, function(x) sum(yenet*x)) %>% max()
## set upper limit for lamba range based on lambdaMax and alpha
## from this website: https://stats.stackexchange.com/questions/144994/range-of-lambda-in-elastic-net-regression
alphas = seq(0, 1, by = 0.1)
lambdaMaxAdj <- 1/(1-alphas)*lambdaMax
lambdaMaxAdj[lambdaMaxAdj == Inf] <- max(lambdaMaxAdj[lambdaMaxAdj != Inf])
enetTuneList <- vector("list", length(alphas))
for(i in 1:length(alphas)){
  enetTuneList[[i]] <- expand.grid(alpha = alphas[i], 
                                   lambda = seq(0, lambdaMaxAdj[i], length.out = 20))
}
enet_grid <- enetTuneList %>% rbindlist() %>% as.data.frame()

## Run cross-validation once
enet = train(x_train, y_train, method = "glmnet",
             metric = "ROC",
             preProcess = c("center", "scale"),
             tuneGrid = enet_grid,
             trControl = train_control)
print(enet)
plot(enet)
enet$bestTune
(enetResults2 <- coef(enet$finalModel, s = enet$bestTune$lambda))

enet3 <- glmnetUtils::glmnet(test_plant_infection ~ ., data = enetTransdata, family = "binomial",
                             alpha = enet$bestTune$alpha, lambda = lambdas)
(enetResults3 <- coef(enet3, s = enet$bestTune$lambda))


#### The cva.glmnet and caret results largely match up
#### Results are highly variable across CV runs. Need to run CV multiple times and average coefficient results
## Run caret::train() through for loop

## Number of runs
nrun <- 500

## Empty vectors for coefficients and best tune parameters
enetResultsList <- enetBestTuneList <- vector("list", nrun)

for(i in 1:nrun){
  enetTrain = train(x_train, y_train, method = "glmnet",
                    metric = "ROC",
                    preProcess = c("center", "scale"),
                    tuneGrid = enet_grid,
                    trControl = train_control)
  enetBestTuneList[[i]] <- enetTrain$bestTune
  enet <- glmnetUtils::glmnet(test_plant_infection ~ ., data = enetTransdata, family = "binomial",
                              alpha = enetTrain$bestTune$alpha, lambda = lambdas)
  enetResultsList[[i]] <- as.data.frame(t(as.matrix(coef(enet, s = enetTrain$bestTune$lambda))))
}

## Combine and summarize best tuning parameters
enetBestTune17 <- enetBestTuneList %>% rbindlist() %>% as.data.frame()
hist(enetBestTune17$alpha)
hist(enetBestTune17$lambda)
## Get mode of alpha and lambda
modeTune17 <- apply(enetBestTune17, 2, mfv)

#### Combine runs and summarize coefficient estimates
enetResultsData <- enetResultsList %>% rbindlist() %>% as.data.frame()
enetSummary17 <- data.frame(meancoef = apply(enetResultsData, 2, mean),
                          mediancoef = apply(enetResultsData, 2, median),
                          sdcoef = apply(enetResultsData, 2, sd))

saveRDS(list(modeTune17, enetSummary17), file = "results/multi-run_CV_elastic_net_transmission_results_2017.rds")


#### Read in 2017 elastic net results
enet17list <- readRDS("results/multi-run_CV_elastic_net_transmission_results_2017.rds")
enetSummary17 <- enet17list[[2]]

#### Plot Elastic Net results
## Vector of clear names for covariates
enetSummary17$niceNames <- factor(c("Intercept", "Resistance trait", "PD disese severity", 
                                    "Proportion vectors infectious", "Xylella population size in infected plant", 
                                    "Leaving rate from infected plant", "Leaving rate from test plant",
                                    "Attraction rate from infected plant", "Attraction rate from test plant"))
levels(enetSummary17$niceNames)
## Reorder factor levels according to coefficient estimate
#enetSummary17$niceNames <- factor(enetSummary17$niceNames, levels = enetSummary17$niceNames[order(enetSummary17$meancoef, decreasing = FALSE)])
## Reorder factor levels according to original order
enetSummary17$niceNames <- with(enetSummary17, factor(niceNames, levels(niceNames)[c(3,8,7,9,6,1,2,4,5)]))
levels(enetSummary17$niceNames)


trans17enetPlot <- ggplot(enetSummary17, aes(y = niceNames, x = meancoef)) +
  geom_errorbarh(aes(xmin = meancoef-sdcoef, xmax = meancoef+sdcoef), colour = "black", height = 0.2) +
  geom_point(size = 3) +
  geom_vline(linetype = "longdash", xintercept = 0) +
  xlab("Coefficient estimate") +
  scale_y_discrete(labels = NULL, name = NULL) +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

trans17enetPlot

ggsave(filename = "results/figures/2017_figures/transmission_2017_elastic_net_coefficients_plot.tiff", 
       plot = trans17enetPlot,
       width = 12, height = 7, units = "in")


#### Plotting transmission vs. acquisition
transAcq17plot <- ggplot(transVCPdata, aes(y = jitter(test_plant_infection, amount = 0.07), x = propInfectious)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, colour = "black") +
  scale_x_continuous(name = "Proportion of vectors infectious",
                     limits = c(0,1),
                     breaks = seq(0,1,0.25)) +
  scale_y_continuous(name = "Probability of transmission",
                     limits = c(-0.1,1.1),
                     breaks = seq(0,1,0.25)) +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

transAcq17plot
ggsave(filename = "results/figures/2017_figures/transmission_acquisition_2017_scatter_plot.tiff", 
       plot = transAcq17plot,
       width = 7, height = 7, units = "in")


#### Transmission vs. mu2 plot
transleave17plot <- ggplot(transVCPdata, aes(y = jitter(test_plant_infection, amount = 0.07), x = mu2)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, colour = "black") +
  xlab("Leaving rate from test plant") + 
  scale_y_continuous(name = "Probability of transmission",
                     limits = c(-0.1,1.1),
                     breaks = seq(0,1,0.25)) +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

transleave17plot
ggsave(filename = "results/figures/2017_figures/transmission_mu2_2017_scatter_plot.tiff", 
       plot = transleave17plot,
       width = 7, height = 7, units = "in")



#####################################################################################################
#### Combining Elastic Net results from 2016 and 2017 into multi-panel figure

enfigure <- plot_grid(trans16enetPlot, trans17enetPlot,
                      align = "h", ncol = 2, nrow = 1, rel_widths = c(2,1),
                      labels = c("(a)", "(b)"), label_x = 0.85, label_y = 0.98,
                      label_size = 10)

enfigure

ggsave(filename = "results/figures/elastic_net_coefficients_both_years_plots.tiff",
       plot = enfigure,
       width = 20, height = 10, units = "cm", dpi = 300, compression = "lzw")



#####################################################################################################
#####################################################################################################
#### Analysis of phytochemsitry data -- phenolics only
#####################################################################################################

#### Cleaning dataset
phenPrefTransData <- readRDS("output/chemistry_preference_transmission_dataset.rds")
summary(phenPrefTransData)
with(phenPrefTransData, table(trt, week))


#### Code for choices
# p1, mu1 = source plant
# p2, mu2 = test plant

#### Log transform response attraction and leaving rates and select only variables of interest
## compounds are in columns 25:ncol
pptData <- phenPrefTransData %>% mutate(logp1 = log(p1), logmu1 = log(mu1)) %>% # Log transform p1 and mu1 for linear models
  dplyr::select(week, PD_symptoms_index, trt, Rep2, xfpop, mu1, p1, logp1, logmu1, 25:ncol(phenPrefTransData))


######################################################################################################
#### MANOVA of total phenolics between resistant and susceptible vines
## Set up formula
manovaDF <- phenPrefTransData %>% dplyr::select(week, trt, 25:ncol(phenPrefTransData)) %>% 
  dplyr::filter(., complete.cases(.)) 
manovaYtotals <- manovaDF %>% dplyr::select(contains("total"))
manovaformulatotals <- as.formula(c("cbind(", paste(names(manovaYtotals), collapse = ", "), ") ~ week*trt"))
manovaformulatotals

manovaMod <- manova(manovaformulatotals, data = manovaDF)
manovaTests <- c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")
for(i in 1:length(manovaTests)){
  print(summary(manovaMod, test = manovaTests[i]))
}
## RESULTS: trt main effect, week main effect, and week:trt interaction are significant
summary.aov(manovaMod) ## Univariate ANCOVAs


## Make a data set for plotting
## Select only totals, gather into a "long" version, and calculate means and SE for each week-trt combination
totalchemMeans <- manovaDF %>% dplyr::select(week, trt, contains("total")) %>%
  gather(key = compound, value = concentration, total.phenolics.stem, total.phenolics.leaf) %>%
  group_by(week, trt, compound) %>% 
  summarise_at("concentration",
               list(~mean(.), se = ~sd(.)/sqrt(length(.)))) %>%
  dplyr::rename(resistance_status = trt)
print.data.frame(totalchemMeans)

#### Plot total phytochemicals over time
totalChemPlot <- ggplot(data=totalchemMeans, aes(x=week, y=mean, shape = compound, linetype = resistance_status,
                                                 group = interaction(compound, resistance_status))) +
  geom_line(colour = "black", size=1.25) +
  geom_point(colour = "black", size=2.5) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), width=0.2, linetype = 1) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = c(2,5,8,14), limits = c(2,14)) + 
  scale_y_continuous(name = "Concentration (ppm)") +
                     #limits = c(0,4000)) +
  # 007 = open circles, dashed line
  # 092 = open triangles, dashed line
  # 094 = closed circles, solid line
  # 102 = closed triangles, solid line
  # scale_shape_manual(values = c(1,2,3)) +
  # scale_linetype_manual(values = c(2,1)) +
  theme_bw(base_size=16) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank())
totalChemPlot

ggsave("results/figures/2017_figures/total_phytochemicals_timeseries_plot.jpg", plot = totalChemPlot,
       width = 10, height = 7, units = "in")


##############################################################################################################
#### Principal Component Analysis of chemistry data and transmission data

## Remove variables of totals and grouped compounds
# List for inclusion of volatiles
# pptData <- pptData %>% dplyr::select(-contains("total"), -monoterpenoids, -green.leafy.volatiles,
#                                      -isoprenoids.and.related, -sequiterpenoids, -early.unknowns, -late.unknowns)

pptData <- pptData %>% dplyr::select(-contains("total"))

#### PCA ordination plots for each time point
figList <- pcaWeekList <- pptWeekData <- pcaANOVAList <- goodPCscores <- pvaluesList <- vector("list", length(unique(pptData$week)))

for(i in 1:length(unique(pptData$week))){
  pptData.i <- pptData %>% 
    dplyr::select(-PD_symptoms_index, -xfpop) %>%
    dplyr::filter(week == unique(week)[i]) %>% 
    dplyr::filter(., complete.cases(.))
  pptWeekData[[i]] <- pptData.i
  chemVars.i <- pptData.i %>% 
    dplyr::select(which(names(pptData.i) == "protocatechuic.acid.hexoside"):ncol(pptData.i))
  pcaWeekList[[i]] <- prcomp(chemVars.i, scale = TRUE)
  pcaVar.i <- pcaWeekList[[i]]$sdev^2
  # Proportion variance explained
  pve.i <- pcaVar.i/sum(pcaVar.i)
  print(paste("week == ", unique(pptData$week)[i]))
  print(round(pve.i, 3))
  # Cummulative variance explained vs. PC
  cumpve.i <- cumsum(pve.i)
  ## Select PCs that explain a cumulative 95% of variance
  goodPCindex.i <- which(round(cumpve.i, digits = 2) <= 0.99)
  goodPCscores.i <- as.data.frame(pcaWeekList[[i]]$x[,goodPCindex.i])
  goodPCscores[[i]] <- names(goodPCscores.i)
  ## cbind good PCs to dat
  pcregData.i <- cbind(pptData.i, goodPCscores.i)
  #print(names(pcregData.i))
  #### Run ANOVAs in a loop
  pcaANOVAList[[i]] <- vector("list", length(goodPCindex.i))
  for(j in 1:length(goodPCindex.i)){
    pc.j <- names(goodPCscores.i)[j]
    anova.j <- aov(pcregData.i[,which(names(pcregData.i) == pc.j)] ~ trt, data = pcregData.i)
    pcaANOVAList[[i]][[j]] <- summary(anova.j)
    #print(pc.j)
    #print(pcaANOVAList[[i]][[j]])
  }
}  
  
#aovweek14list <- pcaANOVAList[[4]]

for(i in 1:length(unique(pptData$week))){
  pcaANOVAList.i <- pcaANOVAList[[i]]
  pvalues <- sapply(1:length(pcaANOVAList.i), function(x) pcaANOVAList.i[[x]][[1]][1,"Pr(>F)"], simplify = TRUE)  
  pvaluesList[[i]] <- pvalues
  nsig <- sum(pvalues < 0.07)
  pvaluesIndexSorted <- order(pvalues)
  sigpvals.i <- c(NA, NA)
  if(nsig == 0){
    sigpvals.i <- 1:2
  }
  if(nsig == 1 & pvaluesIndexSorted[1] != 1){
    sigpvals.i <- c(pvaluesIndexSorted[1], 1)
  }
  if(nsig == 1 & pvaluesIndexSorted[1] == 1){
    sigpvals.i <- pvaluesIndexSorted[1:2]
  }
  if(nsig >= 2){
    sigpvals.i <- pvaluesIndexSorted[1:2]
  }
  figList[[i]] <- ggbiplot(pcaWeekList[[i]], choices = sigpvals.i, obs.scale = 1, var.scale = 1, groups = pptWeekData[[i]]$trt, 
                           ellipse = TRUE, circle = FALSE, ellipse.prob = 0.95, var.axes = FALSE,
                           theme(axis.line = element_line(colour = "black"),
                                 text = element_text(size = 12),
                                 legend.position = "none"))
}


pcaWeekFigure <- plot_grid(figList[[1]], figList[[2]], figList[[3]], figList[[4]],
                           align = "h", ncol = 2, nrow = 2, 
                           labels = c("(a)", "(b)", "(c)", "(d)"), #label_x = 0.85, label_y = 0.98,
                           label_size = 16)

pcaWeekFigure

ggsave(filename = "results/figures/2017_figures/chemistry_pca_weeks_plots.tiff",
       plot = pcaWeekFigure,
       width = 30, height = 30, units = "cm", dpi = 300, compression = "lzw")



##############################################################################################
#### Look at the PCA loadings for week 8
PCloadings8 <- pcaWeekList[[3]]$rotation
pptDataWeek8 <- pptWeekData[[3]]

## PC4 loadings
loadings4 <- data.frame(compounds = attr(PCloadings8[,4], "names"),
                        loadings = as.numeric(PCloadings8[,4]))
loadings4 %>% arrange(-loadings)


## Make sure I'm interpreting the loadings correctly
(chem8Means <- pptDataWeek8 %>% 
    group_by(trt) %>% 
    ## First two compounds should be greater in R (negative loadings); last two compounds should be greater in S (positive loadings)
    summarise_at(c("procyanidin.B.gallate.2", "epicatechin.gallate", "quercitrin", "malvidin.glucoside.oenin"), 
                 list(~mean(.), se = ~sd(.)/sqrt(length(.)))))

##############################################################################################
#### Look at the PCA loadings for week 14
PCloadings14 <- pcaWeekList[[4]]$rotation
pptDataWeek14 <- pptWeekData[[4]]

## PC1 loadings
loadings1 <- data.frame(compounds = attr(PCloadings14[,1], "names"),
                        loadings = as.numeric(PCloadings14[,1]))
loadings1 %>% arrange(-loadings)

## PC2 loadings
loadings2 <- data.frame(compounds = attr(PCloadings14[,2], "names"),
                        loadings = as.numeric(PCloadings14[,2]))
loadings2 %>% arrange(-loadings)


## Make sure I'm interpreting the loadings correctly
(chem14Means <- pptDataWeek14 %>% 
    group_by(trt) %>% 
    ## First three compounds should be greater in R (negative loadings); last compound should be greater in S (positive loadings)
    summarise_at(c("coutaric.acid.1", "deltaviniferin", "piceid", "geraniol"), 
                 list(~mean(.), se = ~sd(.)/sqrt(length(.)))))

###########################################################################################################
#### PC REGRESSION
## PCs vs. genotype, only week 14
## As a post-hoc test, should reduce critical level alpha when evaluating p-values
pptDataWeek14 <- pptWeekData[[4]]
pcaWeek14 <- pcaWeekList[[4]]
summary(pcaWeek14)
str(pcaWeek14)
table(pptDataWeek14$trt)

pcaPlotWeek14 <- ggbiplot(pcaWeek14, choices = 1:2, obs.scale = 1, var.scale = 1, groups = pptDataWeek14$trt, 
                          ellipse = TRUE, circle = FALSE, ellipse.prob = 0.95, var.axes = TRUE,
                          theme(axis.line = element_line(colour = "black"),
                                text = element_text(size = 12)))

pcaPlotWeek14

ggsave(filename = "results/figures/2017_figures/chemistry_pca_week14_plot.tiff",
       plot = pcaPlotWeek14,
       width = 30, height = 30, units = "cm", dpi = 300, compression = "lzw")

## Scree plot of variance captured by PCs
# Variance explained by each component
pcaVar <- pcaWeek14$sdev^2
# Proportion variance explained
pve <- pcaVar/sum(pcaVar)
pve
# Plot variance captured vs. PC
plot(pve, xlab = "Principal Component", ylab = "Proportion of Variance Explained", ylim = c(0,1), type = "b")
## Results: First PC explains 50% of all variance; 2nd PC explains 13%; next 3 are relatively similar, but each explains < 10% of variance
# Cummulative variance explained vs. PC
cumpve <- cumsum(pve)
plot(cumsum(pve), xlab="Principal Component", ylab = "Cumulative Proportion of Variance Explained", 
     ylim=c(0,1), type='b')


#### PCA regression with p1 and mu1
## Select PCs that explain a cumulative 95% of variance
goodPCindex <- which(round(cumpve, digits = 2) <= 0.99)
goodPCscores <- as.data.frame(pcaWeek14$x[,goodPCindex])
## cbind good PCs to dat
chemPC14 <- cbind(pptDataWeek14, goodPCscores)
head(chemPC14)

#### Run ANOVAs in a loop
pcaANOVAList <- vector("list", length(goodPCindex))
for(i in 1:length(goodPCindex)){
  pc.i <- names(goodPCscores)[i]
  anova.i <- aov(chemPC14[,which(names(chemPC14) == pc.i)] ~ trt, data = chemPC14)
  pcaANOVAList[[i]] <- tidy(anova.i, conf.level = 0.95)
  print(pc.i)
  print(pcaANOVAList[[i]])
}
## PC1 and PC2 are significantly different between genotypes
## How are they different?
pca14Summary <- chemPC14 %>% group_by(trt) %>% summarise_at(c("PC1", "PC2"), list(~mean(.), se = ~sd(.)/sqrt(length(.))))

#### Look at the PCA loadings 
PCloadings <- pcaWeek14$rotation

## PC1 loadings
loadings1 <- data.frame(compounds = attr(PCloadings[,1], "names"),
                        loadings = as.numeric(PCloadings[,1]))
loadings1 %>% arrange(-loadings)

## PC2 loadings
loadings2 <- data.frame(compounds = attr(PCloadings[,2], "names"),
                        loadings = as.numeric(PCloadings[,2]))
loadings2 %>% arrange(-loadings)


## Make sure I'm interpreting the loadings correctly
(chem14Means <- pptDataWeek14 %>% 
  group_by(trt) %>% 
  summarise_at(c("coutaric.acid.1", "deltaviniferin", "piceid", "geraniol"), 
               list(~mean(.), se = ~sd(.)/sqrt(length(.)))))

###########################################################################################################
#### PC REGRESSION
## PCs vs. genotype, only week 8
## As a post-hoc test, should reduce critical level alpha when evaluating p-values
pptDataWeek8 <- pptWeekData[[3]]
pcaWeek8 <- pcaWeekList[[3]]
summary(pcaWeek8)
str(pcaWeek8)
table(pptDataWeek8$trt)

pcaPlotWeek8 <- ggbiplot(pcaWeek8, choices = c(1,4), obs.scale = 1, var.scale = 1, groups = pptDataWeek8$trt, 
                          ellipse = TRUE, circle = FALSE, ellipse.prob = 0.95, var.axes = TRUE,
                          theme(axis.line = element_line(colour = "black"),
                                text = element_text(size = 12)))

pcaPlotWeek8

ggsave(filename = "results/figures/2017_figures/chemistry_pca_week8_plot.tiff",
       plot = pcaPlotWeek8,
       width = 30, height = 30, units = "cm", dpi = 300, compression = "lzw")

## Scree plot of variance captured by PCs
# Variance explained by each component
pcaVar <- pcaWeek8$sdev^2
# Proportion variance explained
pve <- pcaVar/sum(pcaVar)
pve
# Plot variance captured vs. PC
plot(pve, xlab = "Principal Component", ylab = "Proportion of Variance Explained", ylim = c(0,1), type = "b")
## Results: First PC explains 30% of all variance; 2nd PC explains 13%; next 2 are relatively similar, but each explains ~ 10% of variance
# Cummulative variance explained vs. PC
cumpve <- cumsum(pve)
cumpve
plot(cumsum(pve), xlab="Principal Component", ylab = "Cumulative Proportion of Variance Explained", 
     ylim=c(0,1), type='b')


#### PCA regression with p1 and mu1
## Select PCs that explain a cumulative 95% of variance
goodPCindex <- which(round(cumpve, digits = 2) <= 0.95)
goodPCscores <- as.data.frame(pcaWeek8$x[,goodPCindex])
## cbind good PCs to dat
chemPC8 <- cbind(pptDataWeek8, goodPCscores)
head(chemPC8)

#### Run ANOVAs in a loop
pcaANOVAList <- vector("list", length(goodPCindex))
for(i in 1:length(goodPCindex)){
  pc.i <- names(goodPCscores)[i]
  anova.i <- aov(chemPC8[,which(names(chemPC8) == pc.i)] ~ trt, data = chemPC8)
  pcaANOVAList[[i]] <- tidy(anova.i, conf.level = 0.95)
  print(pc.i)
  print(pcaANOVAList[[i]])
}
## PC4 is strongly significantly different between genotypes
## How are they different?
pca8Summary <- chemPC8 %>% group_by(trt) %>% summarise_at("PC4", list(~mean(.), se = ~sd(.)/sqrt(length(.))))
pca8Summary

#### Look at the PCA loadings 
PCloadings <- pcaWeek8$rotation

## PC1 loadings
loadings4 <- data.frame(compounds = attr(PCloadings[,4], "names"),
                        loadings = as.numeric(PCloadings[,4]))
loadings4 %>% arrange(-loadings)




#################################################################################################
#### PC REGRESSION with vector attraction and leaving rates
#### Run PCA over all weeks, for PC regression
## Use the untransformed data
## Create a data set of just the chemistry data without totals
chemVars <- pptData %>% dplyr::select(10:ncol(pptData)) %>% dplyr::select(-contains("total"))
names(chemVars)

#### PCA of chemistry variables across all time points
pcaOut <- prcomp(chemVars, scale = TRUE)
summary(pcaOut)
str(pcaOut)
#### Biplot using ggbiplot
## Biplot of PC1 vs PC2
pc12biplot <- ggbiplot(pcaOut, choices = 1:2, obs.scale = 1, var.scale = 1, groups = pptData$trt, 
                       ellipse = TRUE, circle = FALSE, ellipse.prob = 0.95)
pc12biplot
ggsave(pc12biplot, filename = "results/figures/2017_figures/PCA_1-2_biplot_phytochemistry.jpg")


## Biplot of PC1 vs PC3
pc13biplot <- ggbiplot(pcaOut, choices = c(1,3), obs.scale = 1, var.scale = 1, groups = pptData$trt, 
                       ellipse = TRUE, circle = FALSE, ellipse.prob = 0.95)
pc13biplot
## Biplot of PC2 vs PC3
pc23biplot <- ggbiplot(pcaOut, choices = c(2,3), obs.scale = 1, var.scale = 1, groups = pptData$trt, 
                       ellipse = TRUE, circle = FALSE, ellipse.prob = 0.95)
pc23biplot


## Scree plot of variance captured by PCs
# Variance explained by each component
pcaVar <- pcaOut$sdev^2
# Proportion variance explained
pve <- pcaVar/sum(pcaVar)
pve
# Plot variance captured vs. PC
plot(pve, xlab = "Principal Component", ylab = "Proportion of Variance Explained", ylim = c(0,1), type = "b")
## Results: First PC explains 43% of all variance; next 3 PCs are relatively similar, but each explains < 10% of variance
# Cummulative variance explained vs. PC
cumpve <- cumsum(pve)
plot(cumsum(pve), xlab="Principal Component", ylab = "Cumulative Proportion of Variance Explained", 
     ylim=c(0,1), type='b')


#### PCA regression with p1 and mu1
## Select PCs that explain a cumulative 95% of variance
goodPCindex <- which(round(cumpve, digits = 2) <= 0.90)
goodPCscores <- as.data.frame(pcaOut$x[,goodPCindex])
## cbind good PCs to dat
chemPC <- cbind(pptData, goodPCscores)

#### PC Regression of attraction rate
PCattrformulaDF <- chemPC %>% dplyr::select(contains("PC"))
PCattrformula <- as.formula(c("logp1 ~", paste(names(PCattrformulaDF), collapse = "+")))

PCattrMod <- lm(PCattrformula, data = chemPC)
plot(simulateResiduals(PCattrMod))
summary(PCattrMod)
## Results: PC1 and 2 are significantly related to p1; PCs 4,5, and 10 are marginially significant
PCloadings <- pcaOut$rotation

## PC1 loadings
loadings1 <- data.frame(compounds = attr(PCloadings[,1], "names"),
                                     loadings = as.numeric(PCloadings[,1]))
loadings1 %>% arrange(abs(loadings))

## PC2 loadings
loadings2 <- data.frame(compounds = attr(PCloadings[,2], "names"),
                        loadings = as.numeric(PCloadings[,2]))
loadings2 %>% arrange(abs(loadings))


#### PC Regression of leaving rate
PCleaveformula <- as.formula(c("logmu1 ~", paste(names(PCattrformulaDF), collapse = "+")))

PCleaveMod <- lm(PCleaveformula, data = chemPC)
plot(simulateResiduals(PCleaveMod))
summary(PCleaveMod)
## Results: only PC9 is significantly related to mu1
PCloadings <- pcaOut$rotation

## PC9 loadings
loadings9 <- data.frame(compounds = attr(PCloadings[,9], "names"),
                        loadings = as.numeric(PCloadings[,9]))
loadings9 %>% arrange(abs(loadings))





#####################################################################################################
#### Elastic Net analysis of per-cage preference and phenolics

## Set up trainControl
train_control = trainControl(method = "cv",
                             number = 5, returnResamp = "all",
                             savePredictions = "final")

## Number of runs
nrun <- 500


#### Elastic Net of attraction rates 
attrEnetData <- pptData %>% dplyr::filter(., complete.cases(.)) %>%
  dplyr::select(logp1, PD_symptoms_index, 10:ncol(pptData)) 

#### Create a custom tuning grid.
## Find max lambda for cross validation
## Formula found here: https://stats.stackexchange.com/questions/144994/range-of-lambda-in-elastic-net-regression
## Need to center and scale predictors first
xenet <- attrEnetData %>% dplyr::select(-logp1) %>%
  scale(., center = TRUE, scale = TRUE)
yenet <- attrEnetData$logp1 %>% as.numeric()
lambdaMax <- apply(xenet, 2, function(x) sum(yenet*x)) %>% max()
## set upper limit for lamba range based on lambdaMax and alpha
## from this website: https://stats.stackexchange.com/questions/144994/range-of-lambda-in-elastic-net-regression
alphas = seq(0, 1, by = 0.1)
lambdaMaxAdj <- 1/(1-alphas)*lambdaMax
lambdaMaxAdj[lambdaMaxAdj == Inf] <- max(lambdaMaxAdj[lambdaMaxAdj != Inf])
enetTuneList <- vector("list", length(alphas))
for(i in 1:length(alphas)){
  enetTuneList[[i]] <- expand.grid(alpha = alphas[i], 
                                   lambda = seq(0, lambdaMaxAdj[i], length.out = 20))
}
enet_grid <- enetTuneList %>% rbindlist() %>% as.data.frame()

## Run cross-validation once
enet = train(xenet, yenet, method = "glmnet",
             tuneGrid = enet_grid,
             trControl = train_control)
print(enet)
plot(enet)
enet$bestTune
(enetResults2 <- coef(enet$finalModel, s = enet$bestTune$lambda))



#### The cva.glmnet and caret results largely match up
#### Results are highly variable across CV runs. Need to run CV multiple times and average coefficient results
## Run caret::train() through for loop


## Empty vectors for coefficients and best tune parameters
enetResultsList <- enetBestTuneList <- vector("list", nrun)

for(i in 1:nrun){
  enetTrain = train(xenet, yenet, method = "glmnet",
                    metric = "RMSE",
                    tuneGrid = enet_grid,
                    trControl = train_control)
  enetBestTuneList[[i]] <- enetTrain$bestTune
  enetResultsList[[i]] <- as.data.frame(t(as.matrix(coef(enetTrain$finalModel, s = enetTrain$bestTune$lambda))))
}

saveRDS(list(enetBestTuneList, enetResultsList), file = "output/elastic_net_attraction_chemistry_multiple_runs_results.rds")

enetList <- readRDS("output/elastic_net_attraction_chemistry_multiple_runs_results.rds")
enetBestTuneList <- enetList[[1]]
enetResultsList <- enetList[[2]]

## Combine and summarize best tuning parameters
enetBestTuneAttr <- enetBestTuneList %>% rbindlist() %>% as.data.frame()
hist(enetBestTuneAttr$alpha)
hist(enetBestTuneAttr$lambda)
## Get mode of alpha and lambda
(modeTuneAttr <- apply(enetBestTuneAttr, 2, mfv) %>% t() %>% as.data.frame())

#### Combine runs and summarize coefficient estimates
enetResultsData <- enetResultsList %>% rbindlist() %>% as.data.frame()
enetSummaryAttr <- data.frame(compounds = c("Intercept", attr(xenet, "dimnames")[[2]]),
                              meancoef = apply(enetResultsData, 2, mean),
                              mediancoef = apply(enetResultsData, 2, median),
                              sdcoef = apply(enetResultsData, 2, sd))
enetSummaryAttr %>% arrange(abs(meancoef))

#### Plot Elastic Net results
## Reduce the number of compounds for plotting, and remove Intercept because it inflates the scale
plotAttrData <- enetSummaryAttr %>% dplyr::filter(compounds != "Intercept" & abs(mediancoef) > 0)
## Vector of clear names for covariates
levels(plotAttrData$compounds)
## Reorder factor levels according to coefficient estimate
plotAttrData$compounds <- with(plotAttrData, factor(compounds, levels = compounds[order(meancoef, decreasing = FALSE)]))
levels(plotAttrData$compounds)


attrEnetPlot <- ggplot(plotAttrData, aes(y = compounds, x = meancoef)) +
  geom_errorbarh(aes(xmin = meancoef-sdcoef, xmax = meancoef+sdcoef), colour = "black", height = 0.2) +
  geom_point(size = 3) +
  geom_vline(linetype = "longdash", xintercept = 0) +
  scale_x_continuous("Coefficient estimate", limits = c(-0.2,0.1),
                     labels = c(-0.2,-0.1,0,0.1),
                     breaks = c(-0.2,-0.1,0,0.1)) +
  #scale_y_discrete(labels = NULL, name = NULL) +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        text = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

attrEnetPlot

ggsave(filename = "results/figures/2017_figures/attraction_chemistry_elastic_net_coefficients_plot.tiff", 
       plot = attrEnetPlot,
       width = 7, height = 7, units = "in")



#### Which compounds were selected by both PCR and Elastic Net for attraction?
loadings2selected <- loadings2 %>% dplyr::filter(loadings >= 0.1)
plotAttrData[plotAttrData$compounds %in% loadings2selected$compounds,]
loadings2selected[loadings2selected$compounds %in% plotAttrData$compounds,]




#### Elastic Net of leaving rates 
leaveEnetData <- pptData %>% dplyr::filter(., complete.cases(.)) %>%
  dplyr::select(logmu1, PD_symptoms_index, 10:ncol(pptData))

#### Create a custom tuning grid.
## Find max lambda for cross validation
## Formula found here: https://stats.stackexchange.com/questions/144994/range-of-lambda-in-elastic-net-regression
## Need to center and scale predictors first
xenet <- leaveEnetData %>% dplyr::select(-logmu1) %>%
  scale(., center = TRUE, scale = TRUE)
yenet <- leaveEnetData$logmu1 %>% as.numeric()
lambdaMax <- apply(xenet, 2, function(x) sum(yenet*x)) %>% max()
## set upper limit for lamba range based on lambdaMax and alpha
## from this website: https://stats.stackexchange.com/questions/144994/range-of-lambda-in-elastic-net-regression
alphas = seq(0, 1, by = 0.1)
lambdaMaxAdj <- 1/(1-alphas)*lambdaMax
lambdaMaxAdj[lambdaMaxAdj == Inf] <- max(lambdaMaxAdj[lambdaMaxAdj != Inf])
enetTuneList <- vector("list", length(alphas))
for(i in 1:length(alphas)){
  enetTuneList[[i]] <- expand.grid(alpha = alphas[i], 
                                   lambda = seq(0, lambdaMaxAdj[i], length.out = 20))
}
enet_grid <- enetTuneList %>% rbindlist() %>% as.data.frame()

## Run cross-validation once
enet = train(xenet, yenet, method = "glmnet",
             tuneGrid = enet_grid,
             trControl = train_control)
print(enet)
plot(enet)
enet$bestTune
(enetResults2 <- coef(enet$finalModel, s = enet$bestTune$lambda))


#### Results are highly variable across CV runs. Need to run CV multiple times and average coefficient results
## Run caret::train() through for loop


## Empty vectors for coefficients and best tune parameters
enetResultsList2 <- enetBestTuneList2 <- vector("list", nrun)

for(i in 1:nrun){
  enetTrain = train(xenet, yenet, method = "glmnet",
                    metric = "RMSE",
                    tuneGrid = enet_grid,
                    trControl = train_control)
  enetBestTuneList2[[i]] <- enetTrain$bestTune
  enetResultsList2[[i]] <- as.data.frame(t(as.matrix(coef(enetTrain$finalModel, s = enetTrain$bestTune$lambda))))
}

saveRDS(list(enetBestTuneList2, enetResultsList2), file = "output/elastic_net_leaving_chemistry_multiple_runs_results.rds")

enetList <- readRDS("output/elastic_net_leaving_chemistry_multiple_runs_results.rds")
enetBestTuneList2 <- enetList[[1]]
enetResultsList2 <- enetList[[2]]

## Combine and summarize best tuning parameters
enetBestTuneleave <- enetBestTuneList2 %>% rbindlist() %>% as.data.frame()
hist(enetBestTuneleave$alpha)
hist(enetBestTuneleave$lambda)
## Get mode of alpha and lambda
(modeTuneleave <- apply(enetBestTuneleave, 2, mfv) %>% t() %>% as.data.frame())

#### Combine runs and summarize coefficient estimates
enetResultsData <- enetResultsList2 %>% rbindlist() %>% as.data.frame()
enetSummaryleave <- data.frame(compounds = c("Intercept", attr(xenet, "dimnames")[[2]]),
                              meancoef = apply(enetResultsData, 2, mean),
                              mediancoef = apply(enetResultsData, 2, median),
                              sdcoef = apply(enetResultsData, 2, sd))
enetSummaryleave %>% arrange(abs(meancoef))


#### Plot Elastic Net results
## Reduce the number of compounds for plotting, and remove Intercept because it inflates the scale
plotleaveData <- enetSummaryleave %>% dplyr::filter(compounds != "Intercept" & abs(mediancoef) > 0)
## Vector of clear names for covariates
levels(plotleaveData$compounds)
## Reorder factor levels according to coefficient estimate
plotleaveData$compounds <- with(plotleaveData, factor(compounds, levels = compounds[order(meancoef, decreasing = FALSE)]))
levels(plotleaveData$compounds)


leaveEnetPlot <- ggplot(plotleaveData, aes(y = compounds, x = meancoef)) +
  geom_errorbarh(aes(xmin = meancoef-sdcoef, xmax = meancoef+sdcoef), colour = "black", height = 0.2) +
  geom_point(size = 3) +
  geom_vline(linetype = "longdash", xintercept = 0) +
  xlab("Coefficient estimate") +
  #scale_y_discrete(labels = NULL, name = NULL) +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        text = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

leaveEnetPlot

ggsave(filename = "results/figures/2017_figures/leaving_chemistry_elastic_net_coefficients_plot.tiff", 
       plot = leaveEnetPlot,
       width = 7, height = 14, units = "in")




# # Standardize all continuous explanatory variables
# xlasso <- datalasso %>% dplyr::select(-trt, -Rep2, -p1, -mu1) %>% scale(., center = TRUE, scale = TRUE) %>% as.data.frame()
# groupdata <- datalasso %>% dplyr::select(trt, Rep2)
# datalasso2 <- cbind(xlasso, groupdata)
# # Define data classes
# datalasso2$Rep2 <- factor(datalasso2$Rep2)
# datalasso2$trt <- factor(datalasso2$trt)
# datalasso2$trt <- as.numeric(datalasso2$trt) - 1
# # Define trt as a binary numeric variable: R = 0, S = 1
# # Create separate attraction and leaving datasets, and put p1 and mu1 back into datasets
# attrlasso <- leavelasso <- datalasso2
# attrlasso$p1 <- datalasso$p1
# leavelasso$mu1 <- datalasso$mu1
# 
# str(attrlasso)
# str(leavelasso)
# 
# #### Lambda values for cross-validation
# lambdaValues <- seq(500, 0, by = -5)
# 
# # Create formula from column names, rather than writing them all out
# attrformulaDF <- attrlasso %>% dplyr::select(-Rep2, -p1)
# attrformula <- as.formula(c("p1~", paste(names(lassoformulaDF), collapse = "+")))
# attrformula
# 
# # cv.glmmLasso fails with all of the covariates. Need to look into this. Just using the max that will run for now
# cv.attr <- cv.glmmLasso(dat = attrlasso,
#                         form.fixed = p1 ~ week + PD_symptoms_index + ferulic.acid + caftaric.acid.stem + 
#                           procyanidin.B1.stem + catechin.stem + procyanidin.B2.stem + 
#                           epicatechin.stem + pc1.stem + pc3.stem + pc5.stem + pc6.stem + 
#                           pc7.stem + pc8.stem + pc9.stem + pc10.stem + pc11.stem + 
#                           pc12.stem + rutin.stem + res1.stem + res2.stem + res3.stem +
#                           res4.stem + caftaric.acid.leaf + procyanidin.B1.leaf + catechin.leaf + 
#                           procyanidin.B2.leaf + epicatechin.leaf + pc1.leaf + pc3.leaf + 
#                           pc5.leaf + pc6.leaf + pc7.leaf + pc8.leaf + pc9.leaf + pc10.leaf +
#                           pc11.leaf + pc12.leaf,
#                         form.rnd = list(Rep2 = ~1),
#                         lambda = lambdaValues,
#                         family = poisson(link = log))
# 
# bicPath <- cv.attr$BIC_path
# cvRes <- cbind(lambdaValues, bicPath) %>% as.data.frame
# bestLambda <- dplyr::filter(cvRes, bicPath == min(bicPath, na.rm = TRUE))
# 
# 
# attrmod <- glmmLasso(attrformula,
#                      rnd = list(Rep2 = ~1),
#                      lambda = 20,
#                      family = poisson(link = log),
#                      data = attrlasso,
#                      final.re = TRUE,
#                      control = list(print.iter = TRUE))
# 
# 
# summary(attrmod)
# 
# saveRDS(cv.attr, file = "output/cv_attraction_glmmLasso.rds")
# 
# 
# #### Plotting attraction results
# attrResults <- data.frame(params = names(attrmod$coefficients),
#                           estimates = attrmod$coefficients,
#                           SE = attrmod$fixerror) %>% arrange(estimates)
# attrResults$params <- factor(attrResults$params, levels = attrResults$params[order(attrResults$estimates, decreasing = FALSE)])
# 
# attrPlot <- ggplot(attrResults, aes(y = params, x = estimates)) +
#   geom_errorbarh(aes(xmin = estimates-SE, xmax = estimates+SE), colour = "black", height = 0.2) +
#   geom_point(size = 3) +
#   geom_vline(linetype = "longdash", xintercept = 0) +
#   xlab("Coefficient estimate") + ylab("Covariate") + 
#   theme_bw() + 
#   theme(axis.line = element_line(colour = "black"),
#         text = element_text(size = 14),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(colour = "black"),
#         panel.background = element_blank()) 
# attrPlot
# 
# ggsave("results/figures/2017_figures/attraction_phenolics_lasso_plot_2017.jpg", plot = attrPlot,
#        width = 7, height = 7, units = "in")
# 
# 
# #### Analysis of leaving rates
# # Create formula from column names, rather than writing them all out
# leaveformulaDF <- leavelasso %>% dplyr::select(-Rep2, -mu1)
# leaveformula <- as.formula(c("mu1~", paste(names(lassoformulaDF), collapse = "+")))
# leaveformula
# 
# # cv.glmmLasso fails with all of the covariates. Need to look into this. Just using the max that will run for now
# cv.leave <- cv.glmmLasso(dat = leavelasso,
#                          form.fixed = mu1 ~ week + PD_symptoms_index + ferulic.acid + caftaric.acid.stem + 
#                            procyanidin.B1.stem + catechin.stem + procyanidin.B2.stem + 
#                            epicatechin.stem + pc1.stem + pc3.stem + pc5.stem + pc6.stem + 
#                            pc7.stem + pc8.stem + pc9.stem + pc10.stem + pc11.stem + 
#                            pc12.stem + rutin.stem + res1.stem + res2.stem + res3.stem + 
#                            res4.stem + caftaric.acid.leaf + procyanidin.B1.leaf + catechin.leaf + 
#                            procyanidin.B2.leaf + epicatechin.leaf + pc1.leaf + pc3.leaf + 
#                            pc5.leaf + pc6.leaf + pc7.leaf + pc8.leaf + pc9.leaf + pc10.leaf + 
#                            pc11.leaf + pc12.leaf,
#                          form.rnd = list(Rep2 = ~1),
#                          lambda = seq(500, 0, by = -5),
#                          family = poisson(link = log))
# 
# bicPath <- cv.leave$BIC_path
# cvRes <- cbind(lambdaValues, bicPath) %>% as.data.frame
# bestLambda <- dplyr::filter(cvRes, bicPath == min(bicPath, na.rm = TRUE))
# bestLambda
# 
# leavemod <- glmmLasso(leaveformula,
#                       rnd = list(Rep2 = ~1),
#                       lambda = 10,
#                       family = poisson(link = log),
#                       data = leavelasso,
#                       final.re = TRUE,
#                       control = list(print.iter = TRUE))
# 
# 
# summary(leavemod)
# 
# saveRDS(cv.leave, file = "output/cv_leaving_glmmLasso.rds")
# 
# 
# #### Plotting leaveaction results
# leaveResults <- data.frame(params = names(leavemod$coefficients),
#                            estimates = leavemod$coefficients,
#                            SE = leavemod$fixerror) %>% arrange(estimates)
# leaveResults$params <- factor(leaveResults$params, levels = leaveResults$params[order(leaveResults$estimates, decreasing = FALSE)])
# 
# leavePlot <- ggplot(leaveResults, aes(y = params, x = estimates)) +
#   geom_errorbarh(aes(xmin = estimates-SE, xmax = estimates+SE), colour = "black", height = 0.2) +
#   geom_point(size = 3) +
#   geom_vline(linetype = "longdash", xintercept = 0) +
#   xlab("Coefficient estimate") + ylab("Covariate") + 
#   theme_bw() + 
#   theme(axis.line = element_line(colour = "black"),
#         text = element_text(size = 14),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(colour = "black"),
#         panel.background = element_blank()) 
# leavePlot
# 
# ggsave("results/figures/2017_figures/leaving_phenolics_lasso_plot_2017.jpg", plot = leavePlot,
#        width = 7, height = 7, units = "in")
# 
# 
# 
# 
# 
# 
# 
# 
