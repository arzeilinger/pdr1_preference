#### Functions for Parameter Estimation using Maximum Likelihood Estimation

# Source external scripts
#source("R_functions/factor2numeric.R")
#source("R_functions/simulateData.R")

#### AICc function for optimx objects
aicc <- function(k, L, n){
  kl <- length(k)
  aic <- 2*kl + 2*L
  aicc <- aic + ((2*kl*(kl + 1))/(n - kl - 1))
  aicc
}


#### Probability models as functions
P1.func <- function(p1, p2, mu1, mu2, tau, N, n1m1, n2m1){
  # P1 equation as developed in Appendix A
  P1f <- ((1/2)*p1*((2*mu2)/(mu2*p1 + mu1*(mu2 + p2)) + 
                      (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                         (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                       n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                        4*(mu2*p1 + mu1*(mu2 + p2))) + 
                      (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*
                         (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))
  ifelse(P1f < 0, 0, P1f) # ifelse function to prevent underflow
}

P2.func <- function(p1, p2, mu1, mu2, tau, N, n1m1, n2m1){
  # P2 equation as developed in Appendix A
  P2f <- ((mu1*p2)/(mu2*p1 + mu1*(mu2 + p2)) - 
            ((mu1 - mu2 + p1 - p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
               (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                        4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                  (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
            exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                  tau)/(4*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))) - 
            (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                   tau)*(mu1 - mu2 + p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
               (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                       4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                  (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
            (4*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))
  ifelse(P2f < 0, 0, P2f) # ifelse function to prevent underflow
}


# P1.func(p1 = 0.2, p2 = 0.2, mu1 = 0.2, mu2 = 0.2, tau = 20, N = 20, n1m1 = 0, n2m1 = 0)
# P2.func(p1 = 0.2, p2 = 0.2, mu1 = 0.2, mu2 = 0.2, tau = 20, N = 20, n1m1 = 0, n2m1 = 0)
# # Check functions: output from P1.func and P2.func should be equal


###########################################################################
#### Negative log likelihood functions for optimx()

## Fixed Model Variant 
NLL.fixed <- function(par = pars, y = dat){
  # Specify parameters from 'par' vector
  p1f <- par[1]
  p2f <- p1f
  mu1f <- par[2]
  mu2f <- mu1f
  tf <- y$t # Specify vector of time points from data set 'y'
  tm1f <- c(0, tf[1:length(tf)-1]) # specify tm1 ("t minus 1") vector from the tf vector
  tauf <- tf - tm1f # calculate tau = t - tm1
  # specify n1, n2, and n3 from data set 'y'
  n1f <- y$n1 
  n2f <- y$n2
  n3f <- y$n3
  n.all <- cbind(n1f, n2f, n3f)
  # calculate n1m1 ("n1 minus 1") and n2m1 from n1 vector 
  n1m1f <- c(0, n1f[1:length(n1f)-1])
  n2m1f <- c(0, n2f[1:length(n2f)-1])
  Nf <- y$N # specify total sample size, N
  # Calculate P1, P2, P3
  P1 <- P1.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P2 <- P2.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P3 <- 1 - P1 - P2
  p.all <- cbind(P1, P2, P3)
  # Calculate negative log likelihood and correct for possible undeflow or overflow of calculations
  if(all(p.all > 0 & p.all < 1)) 
    -sum(sapply(1:nrow(dat), function(z) dmultinom(x = n.all[z,], prob = p.all[z,], log = TRUE), simplify = TRUE))
  else 1e07
}

## Free Attraction Model variant
NLL.p.choice <- function(par = pars, y = dat){
  # Specify parameters from 'par' vector
  p1f <- par[1]
  p2f <- par[2]
  mu1f <- par[3]
  mu2f <- mu1f
  tf <- y$t # Specify vector of time points from data set 'y'
  tm1f <- c(0, tf[1:length(tf)-1]) # specify tm1 ("t minus 1") vector from the tf vector
  tauf <- tf - tm1f # calculate tau = t - tm1
  # specify n1, n2, and n3 from data set 'y'
  n1f <- y$n1 
  n2f <- y$n2
  n3f <- y$n3
  n.all <- cbind(n1f, n2f, n3f)
  # calculate n1m1 ("n1 minus 1") and n2m1 from n1 vector 
  n1m1f <- c(0, n1f[1:length(n1f)-1])
  n2m1f <- c(0, n2f[1:length(n2f)-1])
  Nf <- y$N # specify total sample size, N
  # Calculate P1, P2, P3
  P1 <- P1.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P2 <- P2.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P3 <- 1 - P1 - P2
  p.all <- cbind(P1, P2, P3)
  # Calculate negative log likelihood and correct for possible undeflow or overflow of calculations
  if(all(p.all > 0 & p.all < 1)) 
    -sum(sapply(1:nrow(dat), function(x) dmultinom(x = n.all[x,], prob = p.all[x,], log = TRUE), simplify = TRUE))
  else 1e07
}

## Free Leaving Model variant
NLL.mu.choice <- function(par = pars, y = dat){
  # Specify parameters from 'par' vector
  p1f <- par[1]
  p2f <- p1f
  mu1f <- par[2]
  mu2f <- par[3]
  tf <- y$t # Specify vector of time points from data set 'y'
  tm1f <- c(0, tf[1:length(tf)-1]) # specify tm1 ("t minus 1") vector from the tf vector
  tauf <- tf - tm1f # calculate tau = t - tm1
  # specify n1, n2, and n3 from data set 'y'
  n1f <- y$n1 
  n2f <- y$n2
  n3f <- y$n3
  n.all <- cbind(n1f, n2f, n3f)
  # calculate n1m1 ("n1 minus 1") and n2m1 from n1 vector 
  n1m1f <- c(0, n1f[1:length(n1f)-1])
  n2m1f <- c(0, n2f[1:length(n2f)-1])
  Nf <- y$N # specify total sample size, N
  # Calculate P1, P2, P3
  P1 <- P1.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P2 <- P2.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P3 <- 1 - P1 - P2
  p.all <- cbind(P1, P2, P3)
  # Calculate negative log likelihood and correct for possible undeflow or overflow of calculations
  if(all(p.all > 0 & p.all < 1)) 
    -sum(sapply(1:nrow(dat), function(x) dmultinom(x = n.all[x,], prob = p.all[x,], log = TRUE), simplify = TRUE))
  else 1e07
}

## Free Choice Model variant
NLL.choice <- function(par = pars, y = dat){
  # Specify parameters from 'par' vector
  p1f <- par[1]
  p2f <- par[2]
  mu1f <- par[3]
  mu2f <- par[4]
  tf <- y$t # Specify vector of time points from data set 'y'
  tm1f <- c(0, tf[1:length(tf)-1]) # specify tm1 ("t minus 1") vector from the tf vector
  tauf <- tf - tm1f # calculate tau = t - tm1
  # specify n1, n2, and n3 from data set 'y'
  n1f <- y$n1 
  n2f <- y$n2
  n3f <- y$n3
  n.all <- cbind(n1f, n2f, n3f)
  # calculate n1m1 ("n1 minus 1") and n2m1 from n1 vector 
  n1m1f <- c(0, n1f[1:length(n1f)-1])
  n2m1f <- c(0, n2f[1:length(n2f)-1])
  Nf <- y$N # specify total sample size, N
  # Calculate P1, P2, P3
  P1 <- P1.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P2 <- P2.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P3 <- 1 - P1 - P2
  p.all <- cbind(P1, P2, P3)
  # Calculate negative log likelihood and correct for possible undeflow or overflow of calculations
  if(all(p.all > 0 & p.all < 1)) 
    -sum(sapply(1:nrow(n.all), function(v) dmultinom(x = n.all[v,], prob = p.all[v,], log = TRUE), simplify = TRUE))
  else 1e07
}


#### List of the NLL functions
NLLlist <- list(NLL.fixed, NLL.p.choice, NLL.mu.choice, NLL.choice)
names(NLLlist) <- c("fixed", "p.choice", "mu.choice", "choice")


#################################################################################################
#### Function to fit data to each model variant and perform model selection with AIC
optimizeCMM <- function(dat = dat, lowerConstraint = 0.0001, upperConstraint = 100, aiccN = dat$N[1]){
  # lowerConstraint = lower inequality constraint for optimizer
  # upperConstraint = upper inequality constraint for optimizer
  # aiccN = N for AICc (corrected for small sample size)
  assign("dat", dat, envir = .GlobalEnv)
  ## Fixed Model Optimization
  op.fixed <- optimx(par = runif(2, lowerConstraint, upperConstraint), # define starting parameter values, randomly chosen from uniform distribution
                     fn = NLL.fixed, y = dat, # define NLL function and data set
                     gr = gr.fixed, # define gradient function
                     method = c("spg"), # define optimization method
                     lower = rep(lowerConstraint, 2), # set lower inequality constraint
                     upper = rep(upperConstraint, 2), # set upper inequality constraint
                     hessian = TRUE,
                     control = list(all.methods=FALSE, 
                                    ftol = 1e-20, maxit = 100000)) # set convergence tolerance and maximum number of iterations
  ## Free Attraction Model Optimization
  op.p.choice <- optimx(par = runif(3, lowerConstraint, upperConstraint), 
                        fn = NLL.p.choice, y = dat, 
                        gr = gr.p.choice,
                        method = c("spg"), 
                        lower = rep(lowerConstraint, 3),
                        upper = rep(upperConstraint, 3),
                        hessian = TRUE,
                        control = list(all.methods=FALSE,
                                       ftol = 1e-20, maxit = 100000))
  ## Free Leaving Model Optimization
  op.mu.choice <- optimx(par = runif(3, lowerConstraint, upperConstraint), 
                         fn = NLL.mu.choice, y = dat, 
                         gr = gr.mu.choice,
                         method = c("spg"), 
                         lower = rep(lowerConstraint, 3),
                         upper = rep(upperConstraint, 3),
                         hessian = TRUE,
                         control = list(all.methods=FALSE, 
                                        ftol = 1e-20, maxit = 100000))
  ## Free Model Optimization
  op.choice <- optimx(par = runif(4, lowerConstraint, upperConstraint), 
                      fn = NLL.choice, y = dat, 
                      gr = gr.choice,
                      method = c("spg"),
                      lower = rep(lowerConstraint, 4),
                      upper = rep(upperConstraint, 4),
                      hessian = TRUE,
                      control = list(all.methods=FALSE, 
                                     ftol = 1e-20, maxit = 100000))
  op.list <- list(op.fixed, op.p.choice, op.mu.choice, op.choice)
  names(op.list) <- c("fixed", "p.choice", "mu.choice", "choice")
  #### Model Selection
  ## AICc calculations for model variants
  ## Using aicc() function from 'source_functions.R' file
  aic.fixed <- aicc(k = coef(op.fixed), L = op.fixed$value, n = aiccN)
  aic.p.choice <- aicc(k = coef(op.p.choice), L = op.p.choice$value, n = aiccN)
  aic.mu.choice <- aicc(k = coef(op.mu.choice), L = op.mu.choice$value, n = aiccN)
  aic.choice <- aicc(k = coef(op.choice), L = op.choice$value, n = aiccN)
  aic.comp <- data.frame("model" = c("fixed", "p.choice", "mu.choice", "choice"),
                         "AICc" = c(aic.fixed, aic.p.choice, aic.mu.choice, aic.choice))
  aic.comp$dAICc <- abs(min(aic.comp$AICc) - aic.comp$AICc)
  aic.comp$df <- c(2, 3, 3, 4)
  # Sort AIC results from lowest to highest dAIC
  modelSelect <- aic.comp[order(aic.comp$dAICc),]
  resultsList <- list(op.list, modelSelect)
  names(resultsList) <- c("op.list", "modelSelect")
  print(dat[1,])
  print(resultsList$modelSelect)
  return(resultsList)
}



#### Function to extract Hessian matrix of best model and look at inter-parameter correlations
getParCorrelations <- function(resultsList = resultsList){
  require(dplyr)
  # resultsList is a list of optimx model fits and model selection, returned from optimizeCMM()
  op.list <- resultsList$op.list
  modelSelect <- resultsList$modelSelect
  # Extract the best model and fitted model object
  bestModelName <- modelSelect[modelSelect$dAICc == 0, "model"]
  bestModelFit <- op.list[[which(names(op.list) == bestModelName)]]
  # Extract parameter estimates into a vector
  bestParams <- bestModelFit[,grep("p", names(bestModelFit))] %>% as.numeric()
  # Extract NLL function from NLLlist, which has the same names as modelSelect and op.list
  #nllFunction <- NLLlist[[which(names(NLLlist) == bestModelName)]]
  Hessian <- attr(bestModelFit, "details")["spg", "nhatend"][[1]]
  #Hessian <- hessian(func = nllFunction, x = bestParams)
  # Calculate variance-covariance matrix
  vcovMatrix <- Hessian %>% solve()
  # Calculate correlation matrix; should have 1's along diagonal
  corMatrix <- vcovMatrix %>% cov2cor()
  output <- list(bestModelName, vcovMatrix, corMatrix)
  names(output) <- c("bestModelName", "vcovMatrix", "corMatrix")
  return(output)
}


#### Function to extract correlation matrix and make it a table with parameters identified
extractCorrelationMatrix <- function(corrMatrixList){
  # corrMatrixList should be the output of getParrCorrelations() function
  # namely, this is a list containing elements:
  # bestModelName = the best model variant, for which vcov and corr matrices are calculated
  # corrMatrix = the correlation matrix from the best model variant, calculated from the variance-covariance matrix
  correlationTable <- corrMatrixList$corMatrix %>% round(., digits = 3) %>% as.data.frame()
  bestModelName <- as.character(corrMatrixList$bestModelName)
  if(bestModelName == "choice") {
    paramNames <- c("p1", "p2", "mu1", "mu2")
    I} else {
      if(bestModelName == "p.choice") {
        paramNames <- c("p1", "p2", "mu")
      } else {
        if(bestModelName == "mu.choice") {
          paramNames <- c("p", "mu1", "mu2")  
        } else {
          paramNames <- c("p", "mu")
        }}}
  names(correlationTable) <- paramNames
  correlationTable <- cbind(paramNames, correlationTable)
  return(correlationTable)
}

##########################################################################################
#### Function to extract parameter estimates and NLL values, 
#### calculate variance based on quadratic approximation method,
#### and return data.frame
# Quadratic approximation based on description in Bolker (2008) pgs. 196 - 201.
# Method can only be used if MLE is close to global maximum.

mleTable <- function(resultsList = resultsList){
  require(optimx); require(tidyr); require(dplyr)
  # resultsList is a list of optimx model fits and model selection, returned from optimizeCMM()
  opList <- resultsList$op.list
  modelSelect <- resultsList$modelSelect
  # Select all models included in selection procedure
  # models <- modelSelect[, "model"]
  # opList <- lapply(models, function(x) op.list[[which(names(op.list) == x)]])
  # names(opList) <- models
  extractFromList <- function(x){
    # p is a vector of estimates for each parameter of the full model
    # hess is the Hessian matrix extracted from the optimx object
    # var is the variance calculated from the Hessian/obvserved information matrix
    # varfull is a vector of the variances for each parameter in the full model
    model <- opList[[x]]
    model.name <- names(opList)[x]
    if(model.name == "fixed") {
        p <- as.numeric(model[,c("p1","p1","p2","p2")])
        hess <- attr(model, "details")["spg", "nhatend"][[1]]
        var <- tryCatch(diag(solve(hess)), 
                        error = function(e) c(-99, -99))
        varfull <- c(var[1], var[1], var[2], var[2])
        if(any(var < 0)){
          warning("negative variance calculated, starting browser")
          #browser()
        }
    } else {
        if(model.name == "p.choice") {
            p <- as.numeric(model[,c("p1","p2","p3","p3")]) 
            hess <- attr(model, "details")["spg", "nhatend"][[1]]
            var <- tryCatch(diag(solve(hess)), 
                            error = function(e) c(-99, -99, -99))
            varfull <- c(var[1], var[2], var[3], var[3])
            if(any(var < 0)){
              warning("negative variance calculated, starting browser")
              #browser()
            }
        } else {
            if(model.name == "mu.choice") {
                p <- as.numeric(model[,c("p1","p1","p2","p3")]) 
                hess <- attr(model, "details")["spg", "nhatend"][[1]]
                var <- tryCatch(diag(solve(hess)), 
                                error = function(e) c(-99, -99, -99))
                varfull <- c(var[1], var[1], var[2], var[3])
                if(any(var < 0)){
                  warning("negative variance calculated, starting browser")
                  #browser()
                }
            } else {
                p <- as.numeric(model[,c("p1","p2","p3","p4")]) # Estimates for Free Choice Model
                hess <- attr(model, "details")["spg", "nhatend"][[1]]
                var <- tryCatch(diag(solve(hess)), 
                                error = function(e) c(-99, -99, -99, -99))
                varfull <- c(var[1], var[2], var[3], var[4])
                if(any(var < 0)){
                  warning("negative variance calculated, starting browser")
                  #browser()
                }
            }}}
    NLL <- model$value # Negative Log Likelihood estimate
    return(c(model.name, p, varfull, NLL))
  }
  moddat <- as.data.frame(t(sapply(1:length(opList), extractFromList, simplify = TRUE)))
  names(moddat) <- c("model", "p1", "p2", "mu1", "mu2", "p1var", "p2var", "mu1var", "mu2var", "NLL")
  moddat <- moddat %>% left_join(., modelSelect, by = "model") %>% arrange(dAICc)
  # Convert factor variables to numeric variables
  moddat$model <- as.character(moddat$model)
  factors <- which(sapply(1:ncol(moddat), function(x) is.factor(moddat[,x]), simplify = TRUE))
  for(i in factors){
    moddat[,i] <- factor2numeric(moddat[,i])
  }
  return(moddat)
}



############################################################################################################
#### Function to average results of good models
averageModels <- function(modelResults = modelResults, dAIC.threshold = 7){
  # modelResults should be a data.frame with a column for each parameter and dAIC
  # calculations appended to original data.frame
  # dAIC.threshold is the threshold for what is considered a good model,
  # convention is dAIC.threshold <= 2, Burnham et al. 2011 Behav Ecol Sociobiol say dAIC.threshold <= 7
  tml <- sum(exp(-modelResults$dAICc/2)) # Total marginal likelihood; requires all models in selection procedure
  # Select only good models
  goodModels <- modelResults[modelResults$dAICc <= dAIC.threshold, ]
  goodModels$gml <- exp(-goodModels$dAICc/2) # Relative or marginal likelihoods for each model
  goodModels$wts <- goodModels$gml/tml # Calculate relative weights for each model
  parNames <- c("p1", "p2", "mu1", "mu2")
  modelAverageFunction <- function(x){
    par.name <- x
    var.name <- paste(par.name, "var", sep = "")
    par.estimates <- as.numeric(goodModels[,par.name])
    parAv <- sum(goodModels$wts*par.estimates)
    var <- as.numeric(goodModels[,c(var.name)])
    varAv <- sum(goodModels$wts*(var + (par.estimates - parAv)^2))
    return(c(par.name, parAv, varAv))
  }
  modelAverage <- parNames %>% sapply(., modelAverageFunction, simplify = TRUE) %>% t() %>% as.data.frame()
  names(modelAverage) <- c("parameter", "estimate", "variance")
  return(modelAverage)
}


