#### Function to run model selection for SEM  

selectSEM <- function(dat){
  #### Function to run model selection for SEM analysis of transmission data
  ## returns Fisher's C statistics and AICc values
  ## models are ordered by no. parameters in descending order
  require(piecewiseSEM)
  ## directPrefModel = additional direct links between feeding preference and transmission
  directPrefModel <- psem(
    glm(transmission ~ acquisition + p2 + mu2 + p1 + mu1, "binomial", dat),
    lm(acquisition ~ p1 + mu1 + XylellaPopulation, dat),
    lm(p1 ~ pd_index, dat),
    lm(mu1 ~ pd_index, dat),
    lm(XylellaPopulation ~ resistanceTrait + week, dat),
    lm(pd_index ~ resistanceTrait + XylellaPopulation, dat),
    p1 %~~% mu2,
    p2 %~~% mu1,
    p1 %~~% p2,
    mu1 %~~% mu2
  )
  outputDirectPref <- summary(directPrefModel, .progressBar = FALSE)
  ## transModel = best hypothesis on processes driving transmission
  transModel <- psem(
    glm(transmission ~ acquisition + p2 + mu2, "binomial", dat),
    lm(acquisition ~ p1 + mu1 + XylellaPopulation, dat),
    lm(p1 ~ pd_index, dat),
    lm(mu1 ~ pd_index, dat),
    lm(XylellaPopulation ~ resistanceTrait + week, dat),
    lm(pd_index ~ resistanceTrait + XylellaPopulation, dat),
    p1 %~~% mu2,
    p2 %~~% mu1,
    p1 %~~% p2,
    mu1 %~~% mu2
  )
  outputTrans <- summary(transModel, .progressBar = FALSE)
  ## attractModel = only attraction rates, no leaving rates included
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
  ## attractModel2 = only attraction rates without pd_index included
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
  ## leaveModel = only leaving rates and no PD symptoms index
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
  ## noPrefModel = no preference components nor PD disease index included
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
  outputList <- list(outputDirectPref, outputTrans, outputAttract, outputAttract2, outputLeave, outputNoPref)
  names(outputList) <- c("directPrefModel", "transModel", "attractModel", "attractModel2", "leaveModel", "noPrefModel")
  aiccScores <- sapply(1:length(outputList), function(x) outputList[[x]]$IC$AICc, simplify = TRUE)
  daicc <- aiccScores - min(aiccScores)
  FCstat <- sapply(1:length(outputList), function(x) outputList[[x]]$Cstat$Fisher.C, simplify = TRUE)
  FCn <- sapply(1:length(outputList), function(x) outputList[[x]]$Cstat$df, simplify = TRUE)
  FCp <- sapply(1:length(outputList), function(x) outputList[[x]]$Cstat$P.Value, simplify = TRUE)
  results <- data.frame(models = c("directPrefModel", "transModel", "attractModel", "attractModel2", "leaveModel", "noPrefModel"),
                        aiccScores,
                        daicc,
                        FCstat,
                        FCn,
                        FCp)
  return(list(results, outputList))
}