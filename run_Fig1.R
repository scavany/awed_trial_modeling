# set working directory
setwd('~/Documents/awed_trial_modeling/')

# clear existing workspace
rm(list = ls())

# load necessary functions
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

# calculate R0 
#R0 <- 1 + life.expectancy * FOI
Sf = S0 * exp(-FOI)
R0 = (log(Sf) - log(S0)) / (Sf - S0)

## Analysis 1 - How does efficacy change with epsilon? ## 
epsilon.vec = seq(0,1,by=0.005)
Nt = Nc = 1
rho.tt = rho.cc = rho_tt_checker(b = 100, d = 1e3)
rho.tc = 1 - rho.tt
rho.ct = 1 - rho.cc
efficacy.epsilon.bestcase = numeric(length(epsilon.vec))
efficacy.epsilon.fullmodel = numeric(length(epsilon.vec))
for(jj in 1:length(epsilon.vec)){
  epsilon = epsilon.vec[jj]
  IAR.t.mosquito <- optim(par = c(0.9), fn = function(par){loss.one(pi = par, S0 = S0, R0 = R0 * (1 - 1 * epsilon))}, lower = c(0), upper = c(1), method = 'Brent')$par - (1-S0)
  IAR.c.mosquito <- optim(par = c(0.9), fn = function(par){loss.one(pi = par, S0 = S0, R0 = R0 * (1 - 0 * epsilon))}, lower = c(0), upper = c(1), method = 'Brent')$par - (1-S0)
  S.f.t.mosquito <- 1 - (IAR.t.mosquito + (1 - S0))
  S.f.c.mosquito <- 1 - (IAR.c.mosquito + (1 - S0))
  
  efficacy.epsilon.bestcase[jj] = 1 - IAR.t.mosquito / IAR.c.mosquito
  efficacy.epsilon.bestcase[jj] = 1 - ((IAR.t.mosquito / IAR.c.mosquito) * (S.f.c.mosquito / S.f.t.mosquito))
  
  IAR.fullmodel = optim(c(0.9,0.9),function(par)
    loss.two(par[1],par[2], rho.tt = rho.tt, rho.tc = rho.tc, rho.cc = rho.cc, rho.ct = rho.ct, Cc = Cc, Ct = Ct, epsilon = epsilon, R0 = R0, S0 = S0),lower=c(0,0),upper=c(1,1),method='BFGS',
    control = list(reltol=1e-12))$par - (1 - S0)
  S.f.t.fullmodel <- 1 - (IAR.fullmodel[1] + (1-S0))
  S.f.c.fullmodel <- 1 - (IAR.fullmodel[2] + (1-S0))
  
  #efficacy.epsilon.fullmodel[jj] = 1 - IAR.fullmodel[1] / IAR.fullmodel[2] 
  efficacy.epsilon.fullmodel[jj] = 1 - ((IAR.fullmodel[1] / IAR.fullmodel[2]) * (S.f.c.fullmodel / S.f.t.fullmodel))
}

plot(epsilon.vec, efficacy.epsilon.bestcase, type = 'l')
lines(epsilon.vec, efficacy.epsilon.fullmodel, lty = 2)
abline(h = c(0.653,0.771,0.849))
lines(0:1, 0:1, col = 'red')

epsilon.fn = approxfun(efficacy.epsilon.fullmodel,epsilon.vec)
efficacy.trial = c(0.653,0.771,0.849)
epsilon.trial = epsilon.fn(efficacy.trial)

epsilon.fn.bestcase = approxfun(efficacy.epsilon.bestcase, epsilon.vec)
epsilon.trial.bestcase = epsilon.fn.bestcase(efficacy.trial)

## Analysis 2 - How does efficacy change with rho 
rho.vec = seq(0.5,1,by=0.001)
efficacy.rho = matrix(0,length(rho.vec),3)
for(jj in 1:3){
  epsilon = epsilon.trial[jj]
  for(ii in 1:length(rho.vec)){
    rho.tt = rho.cc = rho.vec[ii]
    rho.tc = 1 - rho.tt
    rho.ct = 1 - rho.cc
    
    #IAR.t.mosquito <- optim(par = c(0.9), fn = function(par){loss.one(pi = par, S0 = S0, R0 = R0 * (1 - Ct * epsilon))}, lower = c(0), upper = c(1), method = 'Brent')$par - (1-S0)
    #IAR.c.mosquito <- optim(par = c(0.9), fn = function(par){loss.one(pi = par, S0 = S0, R0 = R0 * (1 - Cc * epsilon))}, lower = c(0), upper = c(1), method = 'Brent')$par - (1-S0)
    
    #IAR.t.human <- IAR.t.mosquito * rho.tt + IAR.c.mosquito * rho.tc
    #IAR.c.human <- IAR.c.mosquito * rho.cc + IAR.t.mosquito * rho.ct
    
    IAR = optim(c(0.9,0.9),function(par)
      loss.two(par[1],par[2], rho.tt = rho.tt, rho.tc = rho.tc, rho.cc = rho.cc, rho.ct = rho.ct, Cc = Cc, Ct = Ct, epsilon = epsilon, R0 = R0, S0 = S0),lower=c(0,0),upper=c(1,1),method='BFGS',
      control = list(reltol=1e-12))$par - (1 - S0)
    S.f.t <- 1 - (IAR[1] + (1-S0))
    S.f.c <- 1 - (IAR[2] + (1-S0))
    #efficacy.rho[ii,jj] = 1 - IAR.t.human / IAR.c.human
    efficacy.rho[ii,jj] = 1 - ((IAR[1] / IAR[2]) * (S.f.c / S.f.t))
  }
}

## Analysis 3 - What is the epsilon implied by inferred epsilon as a function of rho
epsilon.implied.vec = seq(0.4,1,by=0.005)
rho.implied.vec = seq(0.75,1,by=0.001)
efficacy.implied = matrix(0,length(rho.implied.vec),length(epsilon.implied.vec))
for(jj in 1:length(epsilon.implied.vec)){
  epsilon = epsilon.implied.vec[jj]
  for(ii in 1:length(rho.implied.vec)){
    
    rho.tt = rho.cc = rho.implied.vec[ii]
    rho.tc = 1 - rho.tt
    rho.ct = 1 - rho.cc
    
    #IAR.t.mosquito <- optim(par = c(0.9), fn = function(par){loss.one(pi = par, S0 = S0, R0 = R0 * (1 - Ct * epsilon))}, lower = c(0), upper = c(1), method = 'Brent')$par - (1-S0)
    #IAR.c.mosquito <- optim(par = c(0.9), fn = function(par){loss.one(pi = par, S0 = S0, R0 = R0 * (1 - Cc * epsilon))}, lower = c(0), upper = c(1), method = 'Brent')$par - (1-S0)
    
    #IAR.t.human <- IAR.t.mosquito * rho.tt + IAR.c.mosquito * rho.tc
    #IAR.c.human <- IAR.c.mosquito * rho.cc + IAR.t.mosquito * rho.ct
    
    IAR = optim(c(0.9,0.9),function(par)
      loss.two(par[1],par[2], rho.tt = rho.tt, rho.tc = rho.tc, rho.cc = rho.cc, rho.ct = rho.ct, Cc = Cc, Ct = Ct, epsilon = epsilon, R0 = R0, S0 = S0),lower=c(0,0),upper=c(1,1),method='BFGS',
      control = list(reltol=1e-12))$par - (1 - S0)
    S.f.t <- 1 - (IAR[1] + (1-S0))
    S.f.c <- 1 - (IAR[2] + (1-S0))
    #efficacy.implied[ii,jj] = 1 - IAR.t.human / IAR.c.human
    efficacy.implied[ii,jj] = 1 - ((IAR[1] / IAR[2]) * (S.f.c / S.f.t))
  }
}

epsilon.implied.fn = approxfun(efficacy.implied[nrow(efficacy.implied),],epsilon.implied.vec)
epsilon.implied.mat = matrix(epsilon.implied.fn(efficacy.implied),nrow(efficacy.implied),ncol(efficacy.implied))
epsilon.implied.trial = epsilon.implied.fn(c(0.653,0.771,0.849))

##consider scale and movement and time at risk 
b.vec = seq(1,1e3,1)
rho.vec = numeric(length(b.vec))
for(ii in 1:length(b.vec)){
  rho.vec[ii] = rho_tt_checker(b.vec[ii],1e3)
}
b.95 = b.vec[which.min(abs(rho.vec-0.95))]
b.90 = b.vec[which.min(abs(rho.vec-0.90))]
b.98 = b.vec[which.min(abs(rho.vec-0.98))]

## Analysis 4 - what would epsilon have been under these different values of b?
rho.b.vec = seq(0.8,1,by=0.001)
epsilon.b.98 = c(epsilon.implied.vec[which.min(abs(
  efficacy.implied[which.min(
    abs(rho.implied.vec-rho_tt_checker(100,1e3))),] -
    efficacy.trial[1]))],
  epsilon.implied.vec[which.min(abs(
    efficacy.implied[which.min(
      abs(rho.implied.vec-rho_tt_checker(100,1e3))),] -
      efficacy.trial[2]))],
  epsilon.implied.vec[which.min(abs(
    efficacy.implied[which.min(
      abs(rho.implied.vec-rho_tt_checker(100,1e3))),] -
      efficacy.trial[3]))])

epsilon.b.90 = c(epsilon.implied.vec[which.min(abs(
  efficacy.implied[which.min(
    abs(rho.implied.vec-rho_tt_checker(b.90,1e3))),] -
    efficacy.trial[1]))],
  epsilon.implied.vec[which.min(abs(
    efficacy.implied[which.min(
      abs(rho.implied.vec-rho_tt_checker(b.90,1e3))),] -
      efficacy.trial[2]))],
  epsilon.implied.vec[which.min(abs(
    efficacy.implied[which.min(
      abs(rho.implied.vec-rho_tt_checker(b.90,1e3))),] -
      efficacy.trial[3]))])

## Analysis 4 - what if the dimensions of the checkerboard were bigger or smaller?
epsilon.delta.vec = epsilon.trial
b.delta.vec = 100
delta.vec = seq(1e2,1e4,length.out=200)
efficacy.delta = matrix(0,3,length(delta.vec))
for(ii in 1:length(epsilon.delta.vec)){
  epsilon = epsilon.delta.vec[ii]
  for(jj in 1:length(delta.vec)){
    delta = delta.vec[jj]
    rho.tt = rho_tt_checker(b.delta.vec,delta)
    rho.cc = rho_cc_checker(b.delta.vec,delta)
    rho.tc = 1 - rho.tt
    rho.ct = 1 - rho.cc

    #IAR.t.mosquito <- optim(par = c(0.9), fn = function(par){loss.one(pi = par, S0 = S0, R0 = R0 * (1 - Ct * epsilon))}, lower = c(0), upper = c(1), method = 'Brent')$par - (1-S0)
    #IAR.c.mosquito <- optim(par = c(0.9), fn = function(par){loss.one(pi = par, S0 = S0, R0 = R0 * (1 - Cc * epsilon))}, lower = c(0), upper = c(1), method = 'Brent')$par - (1-S0)
    
    #IAR.t.human <- IAR.t.mosquito * rho.tt + IAR.c.mosquito * rho.tc
    #IAR.c.human <- IAR.c.mosquito * rho.cc + IAR.t.mosquito * rho.ct
    
    IAR = optim(c(0.9,0.9),function(par)
      loss.two(par[1],par[2], rho.tt = rho.tt, rho.tc = rho.tc, rho.cc = rho.cc, rho.ct = rho.ct, Cc = Cc, Ct = Ct, epsilon = epsilon, R0 = R0, S0 = S0),lower=c(0,0),upper=c(1,1),method='BFGS',
      control = list(reltol=1e-12))$par - (1 - S0)
    
    S.f.t <- 1 - (IAR[1] + (1-S0))
    S.f.c <- 1 - (IAR[2] + (1-S0))
    
    #efficacy.delta[ii,jj] = 1 - (IAR.t.human / IAR.c.human)
    efficacy.delta[ii,jj] = 1 - ((IAR[1] / IAR[2]) * (S.f.c / S.f.t))
  }
}

# plot things about the scale of movement and time at risk
b.vec = seq(1,1e3,1)
rho.vec.by.b = numeric(length(b.vec))
for(ii in 1:length(b.vec)){
  rho.vec.by.b[ii] = rho_tt_checker(b.vec[ii],1e3)
}

# save output for figure generating 
save(efficacy.epsilon.bestcase,
     efficacy.epsilon.fullmodel,
     epsilon.vec,
     efficacy.trial,
     epsilon.trial,
     epsilon.trial.bestcase,
     rho.vec,
     efficacy.rho,
     efficacy.implied,
     epsilon.implied.vec,
     epsilon.implied.mat,
     rho.implied.vec,
     epsilon.implied.mat,
     epsilon.implied.trial,
     delta.vec,
     efficacy.delta,
     b.vec,
     rho.vec.by.b,
     file = './fig_1.RData')

