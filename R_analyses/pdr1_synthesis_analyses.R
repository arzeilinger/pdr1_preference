#### SYNTHETIC ANALYSES OF COMBINED TRANSMISSION-RELATED DATA AND CHEMICAL DATA

rm(list = ls())
# Load packages
my.packages <- c("tidyr", "dplyr", "data.table", "ggplot2", "DHARMa",
                 "bbmle", "glmnetUtils", "caret",
                 "modeest", "cowplot", "broom", "car", "ggbiplot")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/factor2numeric.R")


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



#####################################################################################################
#####################################################################################################
#### Elastic Net analysis of per-cage preference and phenolics
#####################################################################################################

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



####################################################################################################################
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
## Remove Intercept because it inflates the scale
plotleaveData <- enetSummaryleave %>% dplyr::filter(compounds != "Intercept")
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

