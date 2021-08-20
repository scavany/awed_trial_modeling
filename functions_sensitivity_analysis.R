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
  
  # calculate the R0
  if(population_structure == 'exponential')
  {
    R0 <- 1 + (life.expectancy / mean.first.inf)
  } else if (population_structure == 'rectangular')
  {
    R0 <- life.expectancy / mean.first.inf
  }
  
  # calculate the fraction initially susceptible
  S0 <- calc.prop.susc.ci(age.dist = age.dist,
                          FOI = FOI,
                          inv.cross.im = inv.cross.im)
  
  # get the movement parameters
  rho.tc <- 1 - rho.tt
  rho.ct <- 1 - rho.cc
  
  # compute the IARs 
  IAR.t.bestcase <- optim(par = c(0.9), fn = function(par){loss.one(pi = par, S0 = S0, R0 = R0 * (1 - epsilon))}, lower = c(0), upper = c(1), method = 'Brent')$par - (1-S0)
  IAR.c.bestcase <- optim(par = c(0.9), fn = function(par){loss.one(pi = par, S0 = S0, R0 = R0)}, lower = c(0), upper = c(1), method = 'Brent')$par - (1-S0)
  
  IAR.t.mosquito <- optim(par = c(0.9), fn = function(par){loss.one(pi = par, S0 = S0, R0 = R0 * (1 - Ct * epsilon))}, lower = c(0), upper = c(1), method = 'Brent')$par - (1-S0)
  IAR.c.mosquito <- optim(par = c(0.9), fn = function(par){loss.one(pi = par, S0 = S0, R0 = R0 * (1 - Cc * epsilon))}, lower = c(0), upper = c(1), method = 'Brent')$par - (1-S0)
  
  IAR.t.human <- IAR.t.mosquito * rho.tt + IAR.c.mosquito * rho.tc
  IAR.c.human <- IAR.c.mosquito * rho.cc + IAR.t.mosquito * rho.ct
  
  IAR.suppression <- optim(c(0.9,0.9),function(par)
    loss.two(par[1],par[2], rho.tt = rho.tt, rho.tc = rho.tc, rho.cc = rho.cc, rho.ct = rho.ct, Cc = Cc, Ct = Ct, epsilon = epsilon, S0 = S0, R0 = R0),lower=c(0,0),upper=c(1,1),method='BFGS',
    control = list(reltol=1e-12))$par - (1 - S0)
  IAR.t.suppression <- IAR.suppression[1]
  IAR.c.suppression <- IAR.suppression[2]
  
  # calculate efficacy
  eff.bestcase <- 1 - (IAR.t.bestcase / IAR.c.bestcase)
  eff.mosquito <- 1 - (IAR.t.mosquito / IAR.c.mosquito)
  eff.human <- 1 - (IAR.t.human / IAR.c.human)
  eff.suppression <- 1 - (IAR.t.suppression / IAR.c.suppression)
  
  # return output 
  return(c(eff.bestcase = eff.bestcase,
         eff.mosquito = eff.mosquito,
         eff.human = eff.human,
         eff.suppression = eff.suppression))
}


