# load the necessary functions 
source('functions_immune_history.R')
source('functions_trial_sim.R')

# function to compute efficacy as a function of the various parameters
calc_efficacy <- function(age.dist, # age distribution for Indonesia
                          life.expectancy, # age-weighted life expectancy 
                          population_structure = c('exponential','rectangular'), # assumption about population structure for R0 calculation.
                          Ct, # coverage in treated clusters
                          Cc, # control in control clusters
                          FOI, # force of infection
                          epsilon, # factor reduction in R0 from intevention
                          R0, # Can either be a positive real number, or the string "estimate"
                          rho.tt, # proportion of time that treated individuals spend in treated clusters
                          rho.cc = rho.tt, # proportion of time that control individuals spend in control clusters (assume equal for checkerboard)
                          inv.cross.im = 0.5, # assumed that the period of heterologous protection is by default 2 yrs
                          Nt = 1, # population in treated clusters
                          Nc = 1) # population in control clusters
{

  # Get the value of population_structure
  population_structure <- match.arg(population_structure)
    
  # mean time to first infection from one serotype 
  mean.first.inf <- 1 / FOI 
  
  # calculate the fraction initially susceptible
  S0 <- calc.prop.susc.ci(age.dist = age.dist,
                          FOI = FOI,
                          inv.cross.im = inv.cross.im)

  # calculate the R0
  if (R0 == 'estimate_average_age'){
    if(population_structure == 'exponential')
    {
      R0 <- 1 + (life.expectancy / mean.first.inf)
    } else if (population_structure == 'rectangular')
    {
      R0 <- life.expectancy / mean.first.inf
    }
  } else if (R0 == 'estimate_final_size'){
    R0  <-  FOI / (S0 * (1 - exp(-FOI)))
  } else if (!is.numeric(R0) | R0 < 0) {
    stop("R0 should either be a positive real number, or one of the strings 'estimate_final_size' or 'estimate_average_age'")
  }
      
    
  # get the movement parameters
  rho.tc <- 1 - rho.tt
  rho.ct <- 1 - rho.cc
  
  # compute the IARs 
  IAR.t.bestcase <- ifelse(R0 * (1 - epsilon) * S0 > 1, optim(par = S0, fn = function(par){loss.one(pi = par, S0 = S0, R0 = R0 * (1 - epsilon))}, lower = log(R0 * (1 - epsilon) * S0) / R0 / (1 - epsilon), upper = S0, method = 'Brent')$par, 0)
  IAR.c.bestcase <- ifelse(R0 * S0 > 1, optim(par = S0, fn = function(par){loss.one(pi = par, S0 = S0, R0 = R0)}, lower = log(R0 * S0) / R0, upper = S0, method = 'Brent')$par, 0)
  ## Sf.t.bestcase <- S0 - IAR.t.bestcase
  ## Sf.c.bestcase <- S0 - IAR.c.bestcase

  IAR.t.mosquito <- ifelse(R0 * (1 - Ct * epsilon) * S0 > 1, optim(par = S0, fn = function(par){loss.one(pi = par, S0 = S0, R0 = R0 * (1 - Ct * epsilon))}, lower = log(R0 * (1 - epsilon * Ct) * S0) / R0 / (1 - epsilon * Ct), upper = S0, method = 'Brent')$par, 0)
  IAR.c.mosquito <- ifelse(R0 * (1 - Cc * epsilon) * S0 > 1, optim(par = S0, fn = function(par){loss.one(pi = par, S0 = S0, R0 = R0 * (1 - Cc * epsilon))}, lower = log(R0 * (1 - epsilon * Cc) * S0) / R0 / (1 - epsilon * Cc), upper = S0, method = 'Brent')$par, 0)
  ## Sf.t.mosquito <- S0 - IAR.t.mosquito
  ## Sf.c.mosquito <- S0 - IAR.c.mosquito
  
  IAR.t.human <- IAR.t.mosquito * rho.tt + IAR.c.mosquito * rho.tc
  IAR.c.human <- IAR.c.mosquito * rho.cc + IAR.t.mosquito * rho.ct
  ## Sf.t.human <- S0 - IAR.t.human
  ## Sf.c.human <- S0 - IAR.c.human
  
  IAR.suppression <- optim(c(S0,S0),function(par)
    loss.two(par[1], par[2], rho.tt = rho.tt, rho.tc = rho.tc, rho.cc = rho.cc, rho.ct = rho.ct, Cc = Cc, Ct = Ct, epsilon = epsilon, S0 = S0, R0 = R0))$par
  IAR.t.suppression <- IAR.suppression[1]
  IAR.c.suppression <- IAR.suppression[2]
  ## Sf.t.suppression <- S0 - IAR.t.suppression
  ## Sf.c.suppression <- S0 - IAR.c.suppression
  
  # calculate efficacy
  ## eff.bestcase <- 1 - (IAR.t.bestcase / IAR.c.bestcase) * (Sf.c.bestcase / Sf.t.bestcase)
  ## eff.mosquito <- 1 - (IAR.t.mosquito / IAR.c.mosquito) * (Sf.c.mosquito / Sf.t.mosquito)
  ## eff.human <- 1 - (IAR.t.human / IAR.c.human) * (Sf.c.human / Sf.t.human)
  ## eff.suppression <- 1 - (IAR.t.suppression / IAR.c.suppression) * (Sf.c.suppression / Sf.t.suppression)
  eff.bestcase <- 1 - (IAR.t.bestcase / IAR.c.bestcase) * (1 - IAR.c.bestcase) / (1 - IAR.t.bestcase)
  eff.mosquito <- 1 - (IAR.t.mosquito / IAR.c.mosquito) * (1 - IAR.c.mosquito) / (1 - IAR.t.mosquito)
  eff.human <- 1 - (IAR.t.human / IAR.c.human) * (1 - IAR.c.human) / (1 - IAR.t.human)
  eff.suppression <- 1 - (IAR.t.suppression / IAR.c.suppression) * (1 - IAR.c.suppression) / (1 - IAR.t.suppression)

  
  # return output 
  return(c(eff.bestcase = eff.bestcase,
           eff.mosquito = eff.mosquito,
           eff.human = eff.human,
           eff.suppression = eff.suppression))
}


