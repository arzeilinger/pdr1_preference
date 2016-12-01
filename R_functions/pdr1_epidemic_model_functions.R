#### FUNCTIONS FOR SECI EPIDEMIC MODEL


#### Function to specify model equations
# Take the SECI model and then add two parameters determining vector movement and feeding
# phi[p,k] = attraction rate for k = C,I
# phi[mu,k] = leaving rate for k = C,I
# In Madden et al., phi is the number of plants visited per unit time -- is my interpretation of phi valid?
# T = time spent feeding per plant visit
# Host turnover changed to nu (used to be a)
# Vector turnover changed to lambda (used to be mu)


SECIMovementModel <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    # Vector-SECI model with preference decomposed to attraction and leaving
    dS <- nu*(E + C + I) - (phi_muc*(1-exp(-beta_c*T))*V_c*S)/(I + S + E + C) - (phi_mui*(1-exp(-beta_i*T))*V_i*S)/(I + S + E + C)
    dE <- (phi_muc*(1-exp(-beta_c*T))*V_c*S)/(I + S + E + C) + (phi_mui*(1-exp(-beta_i*T))*V_i*S)/(I + S + E + C) - (delta + nu)*E
    dC <- delta*E - nu*C - gamma*C
    dI <- gamma*C - nu*I
    dU <- lambda*(V_c + V_i) - (phi_pc*(1 - exp(-alpha_c*T))*C*U)/(I + S + E + C) - (phi_pi*(1 - exp(-alpha_i*T))*I*U)/(I + S + E + C)
    dV_c <- (phi_pc*(1 - exp(-alpha_c*T))*C*U)/(I + S + E + C) - lambda*V_c
    dV_i <- (phi_pi*(1 - exp(-alpha_i*T))*I*U)/(I + S + E + C) - lambda*V_i
    return(list(c(dS, dE, dC, dI, dU, dV_c, dV_i)))
  })
}



#### Function to return all State variables over all time steps for a single combination of parameters
SECIMDynamics <- function(x){
  x <- as.numeric(x)
  # Submit the parameter values as a row vector in the following order:
  Pars <- c(alpha_c = x[1], alpha_i = x[2], 
            beta_c = x[3], beta_i = x[4],
            phi_pc = x[5], phi_pi = x[6], phi_muc = x[7], phi_mui = x[8], 
            delta = x[9], gamma = x[10], 
            nu = x[11], lambda = x[12], 
            T = x[13])
  State <- c(S = S0, E = E0, C = C0, I = I0, U = U0, V_c = V0, V_i = V0) # Define starting values for states
  Time <- seq(0, 5000, by = 1) # time steps
  model.out <- as.data.frame(ode(func = SECIMovementModel,
                                 y = State,
                                 parms = Pars,
                                 times = Time))
  return(model.out)
}


#### Function to run simulations and output infective host (I) and infectious vector (V) densities
#### For multiple combinations of parameters
SECIMSimulations <- function(x){
  x <- as.numeric(x)
  # Submit the parameter values as a row vector in the following order:
  Pars <- c(alpha_c = x[1], alpha_i = x[2], 
            beta_c = x[3], beta_i = x[4],
            phi_pc = x[5], phi_pi = x[6], phi_muc = x[7], phi_mui = x[8], 
            delta = x[9], gamma = x[10], 
            nu = x[11], lambda = x[12], 
            T = x[13])
  State <- c(S = S0, E = E0, C = C0, I = I0, U = U0, V_c = V0, V_i = V0) # Define starting values for states
  Time <- seq(0, 5000, by = 1) # time steps
  model.out <- as.data.frame(ode(func = SECIMovementModel,
                                 y = State,
                                 parms = Pars,
                                 times = Time))
  # Need to remove some time steps because under/over flow results in NAs
  model.nona <- model.out[!is.na(model.out$I),]
  model.dat <- model.nona[nrow(model.nona),c("time", "C", "I", "V_c", "V_i")]
  #model.dat$trt <- x[14]
  return(model.dat)
}
