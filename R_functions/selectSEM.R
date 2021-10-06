#### Function to run model selection for SEM  

selectSEM <- function(dat){
  #### Function to run model selection for SEM analysis of transmission data
  ## returns Fisher's C statistics and AICc values
  ## models are ordered by no. parameters in descending order
  require(piecewiseSEM); require(dplyr)
  assign("dat", dat, envir = .GlobalEnv)
  
  #############################################################################
  #### Specifying models
  
  #### transModel = best hypothesis on processes driving transmission
  transModel <- psem(
    glm(transmission ~ acquisition + p2 + mu2, "binomial", dat), # Transmission model
    lm(acquisition ~ p1 + mu1 + XylellaPopulation, dat), # Vector acquisition model
    lm(p1 ~ pd_index, dat), # Vector attraction model
    lm(mu1 ~ pd_index, dat), # Vector leaving model
    lm(XylellaPopulation ~ resistanceTrait + week, dat), # Bacterial pop growth model
    lm(pd_index ~ resistanceTrait + XylellaPopulation, dat), # Symptom development model
    ## Specify correlated errors
    p1 %~~% mu2,
    p2 %~~% mu1,
    p1 %~~% p2,
    mu1 %~~% mu2
  )
  outputTrans <- summary(transModel, .progressBar = FALSE)
  
  #### noPDModel = transModel without PD index
  noPDModel <- psem(
    glm(transmission ~ acquisition + p2 + mu2, "binomial", dat),
    lm(acquisition ~ p1 + mu1 + XylellaPopulation, dat),
    lm(XylellaPopulation ~ resistanceTrait + week, dat),
    ## Correlated errors
    p1 %~~% mu2,
    p2 %~~% mu1,
    p1 %~~% p2,
    mu1 %~~% mu2,
    pd_index ~ 1
  )
  outputNoPD <- summary(noPDModel, .progressBar = FALSE)
  
  #### attractModel = only attraction rates, no leaving rates included
  attractModel <- psem(
    glm(transmission ~ acquisition + p2, "binomial", dat),
    lm(acquisition ~ p1 + XylellaPopulation, dat),
    lm(p1 ~ pd_index, dat),
    lm(XylellaPopulation ~ resistanceTrait + week, dat),
    lm(pd_index ~ resistanceTrait + XylellaPopulation, dat),
    p1 %~~% p2,
    mu1 ~ 1,
    mu2 ~ 1
  )
  outputAttract <- summary(attractModel, .progressBar = FALSE)
  
  #### attractModel2 = only attraction rates without pd_index included
  attractModel2 <- psem(
    glm(transmission ~ acquisition + p2, "binomial", dat),
    lm(acquisition ~ p1 + XylellaPopulation, dat),
    lm(XylellaPopulation ~ resistanceTrait + week, dat),
    p1 %~~% p2,
    pd_index ~ 1,
    mu1 ~ 1,
    mu2 ~ 1
  )
  outputAttract2 <- summary(attractModel2, .progressBar = FALSE)
  
  #### leaveModel = only leaving rates and no PD symptoms index
  leaveModel <- psem(
    glm(transmission ~ acquisition + mu2, "binomial", dat),
    lm(acquisition ~ mu1 + XylellaPopulation, dat),
    lm(XylellaPopulation ~ resistanceTrait + week, dat),
    mu1 %~~% mu2,
    pd_index ~ 1,
    p1 ~ 1,
    p2 ~ 1
  )
  outputLeave <- summary(leaveModel, .progressBar = FALSE)
  
  #### noPrefModel = no preference components nor PD disease index included
  noPrefModel <- psem(
    glm(transmission ~ acquisition, "binomial", dat),
    lm(acquisition ~ XylellaPopulation, dat),
    lm(XylellaPopulation ~ resistanceTrait + week, dat),
    p1 ~ 1,
    p2 ~ 1,
    mu1 ~ 1,
    mu2 ~ 1,
    pd_index ~ 1
  )
  outputNoPref <- summary(noPrefModel, .progressBar = FALSE)
  
  #### Compile and summarize results
  outputList <- list(outputTrans, outputNoPD, outputAttract, outputAttract2, outputLeave, outputNoPref)
  modelNames <- c("transModel", "noPDModel", "attractModel", "attractModel2", "leaveModel", "noPrefModel")
  names(outputList) <- modelNames
  aiccScores <- sapply(1:length(outputList), function(x) outputList[[x]]$IC$AICc, simplify = TRUE)
  daicc <- aiccScores - min(aiccScores)
  FCstat <- sapply(1:length(outputList), function(x) outputList[[x]]$Cstat$Fisher.C, simplify = TRUE)
  FCn <- sapply(1:length(outputList), function(x) outputList[[x]]$Cstat$df, simplify = TRUE)
  FCp <- sapply(1:length(outputList), function(x) outputList[[x]]$Cstat$P.Value, simplify = TRUE)
  results <- data.frame(models = modelNames,
                        aiccScores, # AICc scores for each model
                        daicc, # delta AICc
                        FCstat, # Fisher's C statistic
                        FCn, # sample size for Fisher's C
                        FCp) # P-value for Fisher's C
  results <- arrange(results, daicc)
  return(list(results, outputList))
}