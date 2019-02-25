#### SYNTHETIC ANALYSES OF COMBINED TRANSMISSION-RELATED DATA AND CHEMICAL DATA

rm(list = ls())
# Load packages
my.packages <- c("tidyr", "dplyr", "data.table", "ggplot2",
                 "bbmle", "glmmLasso", "lmmen")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/factor2numeric.R")


##############################################################################################################
#### Combine preference, phenolic, and transmission data sets
##############################################################################################################

#### Load transmission, acquisition, culturing, preference data set
transVCPdata <- readRDS("output/complete_2017_transmission-preference_dataset.rds")


## Filter phenolic data set to only Treatment == "Both" as these correspond to the xf_plants or source plants in the trials
phenData <- readRDS("output/full_phenolic_data_pdr1_2017.rds")
phenData <- phenData %>% dplyr::filter(Treatment == "Both")
phenData$Rep <- as.numeric(phenData$Rep)
# Make "Res" groups capitalized
phenData$Res <- with(phenData, ifelse(Res == "r", "R", "S"))
# Add a genotype column for merging with other data sets
phenData$genotype <- with(phenData, ifelse(Res == "R", "094", "092"))


#### Read in phenolic dataset
phenData <- readRDS("output/full_phenolic_data_pdr1_2017.rds")


#### Merge CMM preference parameter data set, phenolic data set, and transmission data set
phenPrefTransData <- transVCPdata %>% 
  left_join(., phenData, by = c("week" = "Week", "genotype", "trt" = "Res", "Rep2" = "Rep")) %>% 
  arrange(week)
str(phenPrefTransData)
summary(phenPrefTransData)

# Save final data set
saveRDS(phenPrefTransData, file = "output/full_phenolics_preference_transmission_dataset.rds")



#####################################################################################################
#####################################################################################################
#### GLMM LASSO analysis of per-cage preference and phenolics

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


######################################################################################################
#### Combined transmission analysis






