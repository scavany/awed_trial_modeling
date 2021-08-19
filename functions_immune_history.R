# install necessary packages
if(!require(deSolve)){install.packages('deSolve'); library(deSolve)}

# FOI model without cross-immunity
prop.exposed.model.noci <- function(t, x, params)
{
  # state variables 
  Zero <- x[1]
  One <- x[2]
  Two <- x[3]
  Three <- x[4]
  Four <- x[5]
  
  # parameters
  FOI <- params[1]
  
  # ODE 
  dZero <- - 4 * FOI * Zero
  dOne <- 4 * FOI * Zero - 3 * FOI * One
  dTwo <- 3 * FOI * One - 2 * FOI * Two
  dThree <- 2 * FOI * Two - 1 * FOI * Three 
  dFour <- 1 * FOI * Three
  
  list(c(dZero, dOne, dTwo, dThree, dFour))
}

# FOI model with cross-immmunity
prop.exposed.model.ci <- function(t, x, params)
{
  # state variables
  Zero <- x[1]
  OneCI <- x[2]
  One <- x[3]
  TwoCI <- x[4]
  Two <- x[5]
  ThreeCI <- x[6]
  Three <- x[7]
  FourCI <- x[8]
  Four <- x[9]
  
  # parameters
  FOI <- params[1]
  inv.cross.im <- params[2]
  
  # ODE 
  dZero <- - 4 * FOI * Zero
  dOneCI <- 4 * FOI * Zero - inv.cross.im * OneCI
  dOne <- inv.cross.im * OneCI - 3 * FOI * One
  dTwoCI <- 3 * FOI * One - inv.cross.im * TwoCI
  dTwo <- inv.cross.im * TwoCI - 2 * FOI * Two
  dThreeCI <- 2 * FOI * Two - inv.cross.im * ThreeCI
  dThree <- inv.cross.im * ThreeCI - 1 * FOI * Three
  dFourCI <- 1 * FOI * Three - inv.cross.im * FourCI
  dFour <- inv.cross.im * FourCI
  
  list(c(dZero, dOneCI, dOne, dTwoCI, dTwo, dThreeCI, dThree, dFourCI, dFour))
}

# function to calculate the proportion susceptible with cross-immunity
calc.prop.susc.ci <- function(age.dist, FOI, inv.cross.im)
{
  seroprob.ci <- ode(c(1,rep(0,8)),age.dist[,1], prop.exposed.model.ci,parms=c(FOI,inv.cross.im))
  prop.by.age <- age.dist[,2] / sum(age.dist[,2])
  prop.susc.ci <- sum(seroprob.ci[,2] * prop.by.age) + (3/4) * sum(seroprob.ci[,4] * prop.by.age) + (2/4) * sum(seroprob.ci[,6] * prop.by.age) + (1/4) * sum(seroprob.ci[,8] * prop.by.age)
  return(prop.susc.ci)
}
