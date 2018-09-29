#### Negative log-likelihood functions for analysis of 2017 PdR1 transmission data
#### Functions inculde linear, Holling Type IV, and Ricker models, 
#### each with a NLL function that estimates across R and S trts and a function that splits estimates, using mle2()'s parameters statement

#### Linear NLL models
linearNLL <- function(a, b, week, nInfected){
  probTrans <- a + b*week
  nll <- -sum(dbinom(nInfected, prob = probTrans, size = 16, log = TRUE))
  ifelse(is.nan(nll), 1e07, nll)
  #return(probTrans)
}

linearNLL.trt <- function(a.r, a.s, b.r, b.s, week, nInfected){
  a <- c(a.r, a.s)[transSummarynl$trt]
  b <- c(b.r, b.s)[transSummarynl$trt]
  probTrans <- a + b*week
  nll <- -sum(dbinom(nInfected, prob = probTrans, size = 16, log = TRUE))
  ifelse(is.nan(nll), 1e07, nll)
  #return(probTrans)
}


#### Holling Type IV NLL models
holling4NLL <- function(a, b, c, week, nInfected){
  probTrans <- (a*week^2)/(b + c*week + week^2)
  nll <- -sum(dbinom(nInfected, prob = probTrans, size = 16, log = TRUE))
  ifelse(is.nan(nll), 1e07, nll)
  #return(probTrans)
}

holling4NLL.trt <- function(a.r, a.s, b.r, b.s, c.r, c.s, week, nInfected){
  a <- c(a.r, a.s)[transSummarynl$trt]
  b <- c(b.r, b.s)[transSummarynl$trt]
  c <- c(c.r, c.s)[transSummarynl$trt]
  probTrans <- (a*week^2)/(b + c*week + week^2)
  nll <- -sum(dbinom(nInfected, prob = probTrans, size = 16, log = TRUE))
  ifelse(is.nan(nll), 1e07, nll)
  #return(probTrans)
}


#### Ricker NLL models
rickerNLL <- function(a, b, week, nInfected){
  probTrans <- a*week*exp(-b*week)
  nll <- -sum(dbinom(nInfected, prob = probTrans, size = 16, log = TRUE))
  ifelse(is.nan(nll), 1e07, nll)
  #return(probTrans)
}

rickerNLL.trt <- function(a.r, a.s, b.r, b.s, week, nInfected){
  a <- c(a.r, a.s)[transSummarynl$trt]
  b <- c(b.r, b.s)[transSummarynl$trt]
  probTrans <- a*week*exp(-b*week)
  nll <- -sum(dbinom(nInfected, prob = probTrans, size = 16, log = TRUE))
  ifelse(is.nan(nll), 1e07, nll)
  #return(probTrans)
}


#### Logistic NLL model
## For now, I will only fit genotype data separately, so I don't need a .trt model too. I may add it later
logisticGrowthNLL <- function(a, b, week, nInfected){
  probTrans <- exp(a + b*week)/(1 + exp(a + b*week))
  nll <- -sum(dbinom(nInfected, prob = probTrans, size = 16, log = TRUE))
  ifelse(is.nan(nll), 1e07, nll)
}


#### Michaelis-Menten NLL model
## For now, I will only fit genotype data separately, so I don't need a .trt model too. I may add it later
MMNLL <- function(a, b, week, nInfected){
  probTrans <- a*week/(b + week)
  nll <- -sum(dbinom(nInfected, prob = probTrans, size = 16, log = TRUE))
  ifelse(is.nan(nll), 1e07, nll)
}


######################################################################################################
#### Function to optimize all three candidate models and return model selection

optimizeTransModels <- function(dat){
  # Linear optimization
  linearOp <-  mle2(linearNLL, data = list(week = dat$week, nInfected = dat$nInfected),
                    start = list(a = a.l0, b = b.l0),
                    optimizer = "optim", method = "Nelder-Mead",
                    control = list(maxit = 10000))
  # Holling Type IV optimization
  holling4Op <- mle2(holling4NLL, data = list(week = dat$week, nInfected = dat$nInfected),
                       start = list(a = a0, b = b0, c = c0),
                       #optimizer = "optimx",
                       control = list(maxit = 10000))
  # Ricker optimization
  rickerOp <- mle2(rickerNLL, data = list(week = dat$week, nInfected = dat$nInfected),
                   start = list(a = a.rick0, b = b.rick0),
                   optimizer = "optim", method = "Nelder-Mead",
                   control = list(maxit = 10000))
  # Logistic Growth optimization
  logisticOp <- mle2(logisticGrowthNLL, data = list(week = dat$week, nInfected = dat$nInfected),
                   start = list(a = a.logist0, b = b.logist0),
                   optimizer = "optim", method = "Nelder-Mead",
                   control = list(maxit = 10000))
  # Michaelis-Menten optimization
  MMOp <- mle2(MMNLL, data = list(week = dat$week, nInfected = dat$nInfected),
                   start = list(a = a.mm0, b = b.mm0),
                   optimizer = "optim", method = "Nelder-Mead",
                   control = list(maxit = 10000))
  # Model selection using AICc
  modelSelect <- ICtab(linearOp, holling4Op, rickerOp, logisticOp, MMOp,
                       type = "AICc", sort = TRUE, delta = TRUE, base = TRUE, 
                       nobs = 16)
  # Return a named list of the optimization results and the model selection
  opList <- list(linearOp, holling4Op, rickerOp, logisticOp, MMOp)
  names(opList) <- c("linearOp", "holling4Op", "rickerOp", "logisticOp", "MMOp")
  return(list(opList, modelSelect))
}

