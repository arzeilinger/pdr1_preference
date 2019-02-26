#### SYNTHETIC ANALYSES OF COMBINED TRANSMISSION-RELATED DATA AND CHEMICAL DATA

rm(list = ls())
# Load packages
my.packages <- c("tidyr", "dplyr", "data.table", "ggplot2", "DHARMa",
                 "bbmle", "glmmLasso", "lmmen", "glmnetUtils", "caret",
                 "modeest")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/factor2numeric.R")


#### Load transmission, acquisition, culturing, preference (and phenolic, when it's ready) data set
transVCPdata <- readRDS("output/complete_2017_transmission-preference_dataset.rds")
str(transVCPdata)
summary(transVCPdata)

#### For CMM parameters
## choice1 = source (Xylella infected) plant
## choice2 = test plant


######################################################################################################
#### Combined transmission analysis 2017
######################################################################################################
FullModel <- glm(test_plant_infection ~ week*genotype + PD_symptoms_index + propInfectious + log10(xfpop+1) + mu1 + mu2 + p1 + p2,
                 data = transVCPdata, family = "binomial")
plot(simulateResiduals(FullModel))
summary(FullModel)

#### Look at relationships among variables
transVCPdata %>% dplyr::select(-week, -block, -genotype, -trt, -rep, -nbugs, -totalInfectious, -plantID, -Rep2) %>% pairs()



#### Transmission analysis using elastic net 
## Define lambda values
lambdas <- 2^seq(-1, -8, length = 20)
#### Using the glmnetUtils package
enetTransdata <- transVCPdata %>% mutate(log.xfpop = log10(xfpop+1),
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
## Create a custom tuning grid.
enet_grid = expand.grid(alpha = seq(0, 1, by = 0.1),
                        lambda = lambdas) # doing an exponential scale allows to explore more space of lambda values
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
nrun <- 100

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
enetBestTune <- enetBestTuneList %>% rbindlist() %>% as.data.frame()
hist(enetBestTune$alpha)
hist(enetBestTune$lambda)
## Get mode of alpha and lambda
modeTune <- apply(enetBestTune, 2, mfv)

#### Combine runs and summarize coefficient estimates
enetResultsData <- enetResultsList %>% rbindlist() %>% as.data.frame()
enetSummary <- data.frame(meancoef = apply(enetResultsData, 2, mean),
                          mediancoef = apply(enetResultsData, 2, median),
                          sdcoef = apply(enetResultsData, 2, sd))







saveRDS(enetResultsList, file = "output/elastic_net_multi_cv_coefficients.rds")




#####################################################################################################
#### GLMM LASSO analysis of per-cage preference and phenolics
#####################################################################################################

#### Cleaning dataset for glmmLASSO
phenPrefTransData <- readRDS("output/full_phenolics_preference_transmission_dataset.rds")
## Select only variables for analysis
datalasso <- phenPrefTransData %>% dplyr::select(week, trt, Rep2, p1, mu1, PD_symptoms_index, 
                                                 ferulic.acid, contains(".stem"), contains(".leaf"), 
                                                 contains("pc"), contains("res")) %>%
  dplyr::filter(., complete.cases(.))
# Standardize all continuous explanatory variables
xlasso <- datalasso %>% dplyr::select(-trt, -Rep2, -p1, -mu1) %>% scale(., center = TRUE, scale = TRUE) %>% as.data.frame()
groupdata <- datalasso %>% dplyr::select(trt, Rep2)
datalasso2 <- cbind(xlasso, groupdata)
# Define data classes
datalasso2$Rep2 <- factor(datalasso2$Rep2)
datalasso2$trt <- factor(datalasso2$trt)
datalasso2$trt <- as.numeric(datalasso2$trt) - 1
# Define trt as a binary numeric variable: R = 0, S = 1
# Create separate attraction and leaving datasets, and put p1 and mu1 back into datasets
attrlasso <- leavelasso <- datalasso2
attrlasso$p1 <- datalasso$p1
leavelasso$mu1 <- datalasso$mu1

str(attrlasso)
str(leavelasso)

#### Lambda values for cross-validation
lambdaValues <- seq(500, 0, by = -5)

# Create formula from column names, rather than writing them all out
attrformulaDF <- attrlasso %>% dplyr::select(-Rep2, -p1)
attrformula <- as.formula(c("p1~", paste(names(lassoformulaDF), collapse = "+")))
attrformula

# cv.glmmLasso fails with all of the covariates. Need to look into this. Just using the max that will run for now
cv.attr <- cv.glmmLasso(dat = attrlasso,
                        form.fixed = p1 ~ week + PD_symptoms_index + ferulic.acid + caftaric.acid.stem + 
                          procyanidin.B1.stem + catechin.stem + procyanidin.B2.stem + 
                          epicatechin.stem + pc1.stem + pc3.stem + pc5.stem + pc6.stem + 
                          pc7.stem + pc8.stem + pc9.stem + pc10.stem + pc11.stem + 
                          pc12.stem + rutin.stem + res1.stem + res2.stem + res3.stem +
                          res4.stem + caftaric.acid.leaf + procyanidin.B1.leaf + catechin.leaf + 
                          procyanidin.B2.leaf + epicatechin.leaf + pc1.leaf + pc3.leaf + 
                          pc5.leaf + pc6.leaf + pc7.leaf + pc8.leaf + pc9.leaf + pc10.leaf +
                          pc11.leaf + pc12.leaf,
                        form.rnd = list(Rep2 = ~1),
                        lambda = lambdaValues,
                        family = poisson(link = log))

bicPath <- cv.attr$BIC_path
cvRes <- cbind(lambdaValues, bicPath) %>% as.data.frame
bestLambda <- dplyr::filter(cvRes, bicPath == min(bicPath, na.rm = TRUE))


attrmod <- glmmLasso(attrformula,
                     rnd = list(Rep2 = ~1),
                     lambda = 20,
                     family = poisson(link = log),
                     data = attrlasso,
                     final.re = TRUE,
                     control = list(print.iter = TRUE))


summary(attrmod)

saveRDS(cv.attr, file = "output/cv_attraction_glmmLasso.rds")


#### Plotting attraction results
attrResults <- data.frame(params = names(attrmod$coefficients),
                          estimates = attrmod$coefficients,
                          SE = attrmod$fixerror) %>% arrange(estimates)
attrResults$params <- factor(attrResults$params, levels = attrResults$params[order(attrResults$estimates, decreasing = FALSE)])

attrPlot <- ggplot(attrResults, aes(y = params, x = estimates)) +
  geom_errorbarh(aes(xmin = estimates-SE, xmax = estimates+SE), colour = "black", height = 0.2) +
  geom_point(size = 3) +
  geom_vline(linetype = "longdash", xintercept = 0) +
  xlab("Coefficient estimate") + ylab("Covariate") + 
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 
attrPlot

ggsave("results/figures/2017_figures/attraction_phenolics_lasso_plot_2017.jpg", plot = attrPlot,
       width = 7, height = 7, units = "in")


#### Analysis of leaving rates
# Create formula from column names, rather than writing them all out
leaveformulaDF <- leavelasso %>% dplyr::select(-Rep2, -mu1)
leaveformula <- as.formula(c("mu1~", paste(names(lassoformulaDF), collapse = "+")))
leaveformula

# cv.glmmLasso fails with all of the covariates. Need to look into this. Just using the max that will run for now
cv.leave <- cv.glmmLasso(dat = leavelasso,
                         form.fixed = mu1 ~ week + PD_symptoms_index + ferulic.acid + caftaric.acid.stem + 
                           procyanidin.B1.stem + catechin.stem + procyanidin.B2.stem + 
                           epicatechin.stem + pc1.stem + pc3.stem + pc5.stem + pc6.stem + 
                           pc7.stem + pc8.stem + pc9.stem + pc10.stem + pc11.stem + 
                           pc12.stem + rutin.stem + res1.stem + res2.stem + res3.stem + 
                           res4.stem + caftaric.acid.leaf + procyanidin.B1.leaf + catechin.leaf + 
                           procyanidin.B2.leaf + epicatechin.leaf + pc1.leaf + pc3.leaf + 
                           pc5.leaf + pc6.leaf + pc7.leaf + pc8.leaf + pc9.leaf + pc10.leaf + 
                           pc11.leaf + pc12.leaf,
                         form.rnd = list(Rep2 = ~1),
                         lambda = seq(500, 0, by = -5),
                         family = poisson(link = log))

bicPath <- cv.leave$BIC_path
cvRes <- cbind(lambdaValues, bicPath) %>% as.data.frame
bestLambda <- dplyr::filter(cvRes, bicPath == min(bicPath, na.rm = TRUE))
bestLambda

leavemod <- glmmLasso(leaveformula,
                      rnd = list(Rep2 = ~1),
                      lambda = 10,
                      family = poisson(link = log),
                      data = leavelasso,
                      final.re = TRUE,
                      control = list(print.iter = TRUE))


summary(leavemod)

saveRDS(cv.leave, file = "output/cv_leaving_glmmLasso.rds")


#### Plotting leaveaction results
leaveResults <- data.frame(params = names(leavemod$coefficients),
                           estimates = leavemod$coefficients,
                           SE = leavemod$fixerror) %>% arrange(estimates)
leaveResults$params <- factor(leaveResults$params, levels = leaveResults$params[order(leaveResults$estimates, decreasing = FALSE)])

leavePlot <- ggplot(leaveResults, aes(y = params, x = estimates)) +
  geom_errorbarh(aes(xmin = estimates-SE, xmax = estimates+SE), colour = "black", height = 0.2) +
  geom_point(size = 3) +
  geom_vline(linetype = "longdash", xintercept = 0) +
  xlab("Coefficient estimate") + ylab("Covariate") + 
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 
leavePlot

ggsave("results/figures/2017_figures/leaving_phenolics_lasso_plot_2017.jpg", plot = leavePlot,
       width = 7, height = 7, units = "in")








