# set working directory
setwd('~/Documents/awed_trial_modeling/')

# clear existing workspace
rm(list = ls())

# load functions
source('functions_trial_sim.R')
source('functions_immune_history.R')

# load the coverage data and calculate mean coverage 
dat_cov <- read.csv('AWED_Data.csv')
Ct <- mean(dat_cov$wMel[dat_cov$Treated == 1])
Cc <- mean(dat_cov$wMel[dat_cov$Treated == 0])

# load the age distribution 
pop.age = read.csv('./pop_by_age_Indonesia.csv')
age.dist = cbind(
  rep(pop.age$AgeGrpStart,each=pop.age$AgeGrpSpan) + (0:4),
  rep(pop.age$PopTotal/5,each=pop.age$AgeGrpSpan))
age.dist[,2] = age.dist[,2] * 1000

# calculate the life expectancy as the weighted average of the life expectancies across age groups
life.expectancy <- sum(pop.age[,ncol(pop.age)] * (age.dist[,2] / sum(age.dist[,2])))

# specify transmission parameters
FOI = 0.0317767
inv.cross.im <- 0.5

# calculate the initially susceptible 
S0 = calc.prop.susc.ci(age.dist = age.dist,
                       FOI = FOI,
                       inv.cross.im = inv.cross.im)
# S0 = 1 - 1e-7 #slightly less than 1 or optim solver breaks

# calculate R0 
#R0 <- 1 + life.expectancy * FOI
Sf = S0 * exp(-FOI)
R0 = 3.5#(log(Sf) - log(S0)) / (Sf - S0)

# look at effect on IAR 
## rho.vec = seq(0.8,1,by=0.001)
## b.98 = 10
## b.90 = which.min(abs(rho.vec-0.90))
## b.vec = c(b.98,b.90)
load('./fig_1.RData',verbose=TRUE)
epsilon.vec = seq(0,1,by=0.005)
epsilon <- approxfun(efficacy.epsilon.fullmodel,epsilon.vec)(efficacy.trial[2])
Nt = 1
delta.vec  = seq(0,1e4,by=50)
Nc = Nc_checker(delta.vec)
rho.tt.vec = rho_tt_checker(60,delta.vec)
rho.cc.vec = rho_cc_checker(60,delta.vec)
rho.tc.vec = 1 - rho.tt
rho.ct.vec = 1 - rho.cc
IAR.bestcase = IAR.mosquito = IAR.human = IAR.hummoz = IAR.humsupp = IAR.fullmodel = matrix(0,2,length(delta.vec))
## S.f.bestcase = S.f.mosquito = S.f.human = S.f.fullmodel = matrix(0,2,length(epsilon.vec))

for(jj in 1:length(delta.vec)){
  rho.tt <- rho.tt.vec[jj]
  rho.cc <- rho.cc.vec[jj]
  rho.tc <- rho.tc.vec[jj]
  rho.ct <- rho.ct.vec[jj]
    
  IAR.bestcase[1,jj] <- ifelse(R0 * (1 - epsilon) * S0 > 1, optim(par = S0, fn = function(par){loss.one(pi = par, S0 = S0, R0 = R0 * (1 - epsilon))}, lower = log(R0 * (1 - epsilon) * S0) / R0 / (1 - epsilon), upper = S0, method = 'Brent')$par, 0)
  IAR.bestcase[2,jj] <- ifelse(R0 * S0 > 1, optim(par = S0, fn = function(par){loss.one(pi = par, S0 = S0, R0 = R0)}, lower = log(R0 * S0) / R0, upper = S0, method = 'Brent')$par, 0)
  ## S.f.bestcase[1,jj] <- S0 - IAR.bestcase[1,jj]
  ## S.f.bestcase[2,jj] <- S0 - IAR.bestcase[2,jj]
  
  
  IAR.mosquito[1,jj] <- ifelse(R0 * (1 - Ct * epsilon) * S0 > 1, optim(par = S0, fn = function(par){loss.one(pi = par, S0 = S0, R0 = R0 * (1 - Ct * epsilon))}, lower = log(R0 * (1 - epsilon * Ct) * S0) / R0 / (1 - epsilon * Ct), upper = S0, method = 'Brent')$par, 0)
  IAR.mosquito[2,jj] <- ifelse(R0 * (1 - Cc * epsilon) * S0 > 1, optim(par = S0, fn = function(par){loss.one(pi = par, S0 = S0, R0 = R0 * (1 - Cc * epsilon))}, lower = log(R0 * (1 - epsilon * Cc) * S0) / R0 / (1 - epsilon * Cc), upper = S0, method = 'Brent')$par, 0)
  ## S.f.mosquito[1,jj] <- S0 - IAR.mosquito[1,jj]
  ## S.f.mosquito[2,jj] <- S0 - IAR.mosquito[2,jj]

  IAR.human[1,jj] <- IAR.bestcase[1,jj] * rho.tt + IAR.bestcase[2,jj] * rho.tc
  IAR.human[2,jj] <- IAR.bestcase[2,jj] * rho.cc + IAR.bestcase[1,jj] * rho.ct
  
  IAR.hummoz[1,jj] <- IAR.mosquito[1,jj] * rho.tt + IAR.mosquito[2,jj] * rho.tc
  IAR.hummoz[2,jj] <- IAR.mosquito[2,jj] * rho.cc + IAR.mosquito[1,jj] * rho.ct
  ## S.f.human[1,jj] <- S0 - IAR.human[1,jj]
  ## S.f.human[2,jj] <- S0 - IAR.human[2,jj]
  
  IAR.humsupp[,jj] = plogis(optim(qlogis(c(S0,S0)),function(par)
    loss.two(plogis(par[1]),plogis(par[2]), rho.tt = rho.tt, rho.tc = rho.tc, rho.cc = rho.cc, rho.ct = rho.ct, Cc = 0, Ct = 1, epsilon = epsilon, S0 = S0, R0 = R0))$par)
  
  IAR.fullmodel[,jj] = plogis(optim(qlogis(c(S0,S0)),function(par)
    loss.two(plogis(par[1]),plogis(par[2]), rho.tt = rho.tt, rho.tc = rho.tc, rho.cc = rho.cc, rho.ct = rho.ct, Cc = Cc, Ct = Ct, epsilon = epsilon, S0 = S0, R0 = R0))$par)
  ## S.f.fullmodel[1,jj] <- S0 - IAR.fullmodel[1,jj]
  ## S.f.fullmodel[2,jj] <- S0 - IAR.fullmodel[2,jj]
}

# compute efficacy
## ARR method
## efficacy.bestcase <- 1 - (IAR.bestcase[1,] / IAR.bestcase[2,])
## efficacy.mosquito <- 1 - (IAR.mosquito[1,] / IAR.mosquito[2,])
## efficacy.human <- 1 - (IAR.human[1,] / IAR.human[2,])
## efficacy.fullmodel <- 1 - (IAR.fullmodel[1,] / IAR.fullmodel[2,])
## OR1 method
## efficacy.bestcase <- 1 - ((IAR.bestcase[1,] / IAR.bestcase[2,]) * (S.f.bestcase[2,] / S.f.bestcase[1,]))
## efficacy.mosquito <- 1 - ((IAR.mosquito[1,] / IAR.mosquito[2,]) * (S.f.mosquito[2,] / S.f.mosquito[1,]))
## efficacy.human <- 1 - ((IAR.human[1,] / IAR.human[2,]) * (S.f.human[2,] / S.f.human[1,]))
## efficacy.fullmodel <- 1 - ((IAR.fullmodel[1,] / IAR.fullmodel[2,]) * (S.f.fullmodel[2,] / S.f.fullmodel[1,]))
## OR2 method
efficacy.bestcase <- 1 - (IAR.bestcase[1,] / IAR.bestcase[2,]) * (1 - IAR.bestcase[2,]) / (1 - IAR.bestcase[1,])
efficacy.mosquito <- 1 - (IAR.mosquito[1,] / IAR.mosquito[2,]) * (1 - IAR.mosquito[2,]) / (1 - IAR.mosquito[1,])
efficacy.human <- 1 - (IAR.human[1,] / IAR.human[2,]) * (1 - IAR.human[2,]) / (1 - IAR.human[1,])
efficacy.hummoz <- 1 - (IAR.hummoz[1,] / IAR.hummoz[2,]) * (1 - IAR.hummoz[2,]) / (1 - IAR.hummoz[1,])
efficacy.humsupp <- 1 - (IAR.humsupp[1,] / IAR.humsupp[2,]) * (1 - IAR.humsupp[2,]) / (1 - IAR.humsupp[1,])
efficacy.fullmodel <- 1 - (IAR.fullmodel[1,] / IAR.fullmodel[2,]) * (1 - IAR.fullmodel[2,]) / (1 - IAR.fullmodel[1,])

# compute the fraction of the bias attributable to each phenomenon
# assuming the following embedding, in this case: bestcase -> human -> hummoz -> fullmodel
frac.bias.human <-  (efficacy.bestcase - efficacy.human) / (efficacy.bestcase - efficacy.fullmodel)
frac.bias.hummoz <- (efficacy.human - efficacy.hummoz) / (efficacy.bestcase - efficacy.fullmodel)
frac.bias.fullmodel <- 1 - frac.bias.hummoz - frac.bias.human

# save output 
save(delta.vec,
     IAR.bestcase,
     IAR.mosquito,
     IAR.human,
     IAR.hummoz,
     IAR.humsupp,
     IAR.fullmodel,
     efficacy.bestcase,
     efficacy.mosquito,
     efficacy.human,
     efficacy.hummoz,
     efficacy.humsupp,
     efficacy.fullmodel,
     frac.bias.human,
     frac.bias.hummoz,
     frac.bias.fullmodel,
     file = './fig_2_delta.RData')
