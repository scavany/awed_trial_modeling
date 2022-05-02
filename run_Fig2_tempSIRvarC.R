# set working directory
setwd('~/Documents/awed_trial_modeling/')

# clear existing workspace
rm(list = ls())

# install necessary packages
if(!require(pacman)){install.packages('pacman'); library(pacman)}
p_load(deSolve,data.table)

# load functions
source('functions_trial_sim.R')
source('functions_immune_history.R')

# load the coverage data and calculate mean coverage 
dat_cov <- fread('AWED_Data.csv')
Ct <- dat_cov[Treated == 1,.(Cluster,wMel,Date)]
Cc <- dat_cov[Treated == 0,.(Cluster,wMel,Date)]
Ct.clusters <- unique(Ct$Cluster)
Cc.clusters <- unique(Cc$Cluster)
Ct.mean <- Ct[,.(wMel.mean=mean(wMel)),by=Date]
Cc.mean <- Cc[,.(wMel.mean=mean(wMel)),by=Date]
setorder(Ct.mean,Date)
setorder(Cc.mean,Date)
Ct.fn <- approxfun(as.numeric(Ct.mean$Date - min(Ct.mean$Date)),Ct.mean$wMel.mean)
Cc.fn <- approxfun(as.numeric(Cc.mean$Date - min(Cc.mean$Date)),Cc.mean$wMel.mean)
wMel.start.date <- min(Ct$Date)

# load the age distribution 
pop.age = read.csv('./pop_by_age_Indonesia.csv')
age.dist = cbind(
  rep(pop.age$AgeGrpStart,each=pop.age$AgeGrpSpan) + (0:4),
  rep(pop.age$PopTotal/5,each=pop.age$AgeGrpSpan))
age.dist[,2] = age.dist[,2] * 1000

# calculate the life expectancy as the weighted average of the life expectancies across age groups
## life.expectancy <- sum(pop.age[,ncol(pop.age)] * (age.dist[,2] / sum(age.dist[,2])))
life.expectancy <- pop.age[1,"LifeExpectancyAtBirth"]

## specify transmission parameters
load("foi_yogya.RData")
FOI = foi.yogya
inv.cross.im <- 0.5

# calculate the initially susceptible 
valency = calc.prop.susc.ci.byvalency(age.dist = age.dist,
                                     FOI = FOI,
                                     inv.cross.im = inv.cross.im)
S0 <- sum(valency[c(1,3,5,7)] * c(1,.75,.5,.25))

# calculate R0 
#R0 <- 1 + life.expectancy * FOI
Sf = S0 * exp(-4*FOI)
R0 = (log(Sf) - log(S0)) / (Sf - S0)
## R0 <- R0a.val
## R0 <- R0b.val

## Other parms
gamma = 1/7 #recovery rate
gamma.wan <- gamma * (inv.cross.im/365.25) / (gamma - inv.cross.im/365.25) #CI waning rate
mu <- 1/(life.expectancy*365.25) #mortality rate
tvec.full <- seq(-1,366*11,1)
day0 <- as.Date("2006-01-01")
dvec <- tvec.full + day0
start.date <- as.Date("2006-01-01"); end.date <- as.Date("2007-12-31")
tvec <- tvec.full[which(dvec >= start.date & dvec <= end.date)]
load("./tempSIR_optimout.RData",verbose=TRUE)
beta.amp <- model.fit$par[["beta.amp"]]
offset <- model.fit$par[["offset"]]
load("./yearly_initial_conditions.RData",verbose=TRUE)
state.init <- as.numeric(yearly.ics[year==year(start.date),-1])
names(state.init) <- colnames(yearly.ics)[-1]
state.init["cum.inc"] <- 0

# look at effect on IAR 
rho.vec = seq(0.8,1,by=0.001)
b.98 = 10
b.90 = which.min(abs(rho.vec-0.90))
b.vec = c(b.98,b.90)
epsilon.vec = seq(0,1,by=0.005)
Nt = 1
delta  = 1e3
Nc = Nc_checker(delta)
rho.tt = rho_tt_checker(44,delta)
rho.cc = rho_cc_checker(44,delta)
rho.tc = 1 - rho.tt
rho.ct = 1 - rho.cc
IAR.bestcase = IAR.mosquito = IAR.human = IAR.hummoz = IAR.humsupp = IAR.fullmodel = matrix(0,2,length(epsilon.vec))
## S.f.bestcase = S.f.mosquito = S.f.human = S.f.fullmodel = matrix(0,2,length(epsilon.vec))

for(jj in 1:length(epsilon.vec)){
  print(jj)
  epsilon = epsilon.vec[jj]
  
  state.t.bestcase <- state.init
  parms.t.bestcase <- c(R0=R0,epsilon=epsilon,C=1,N=1,gamma=gamma,gamma.wan=gamma.wan,
                        beta.amp=beta.amp,
                        mu=mu,offset=offset)
  IAR.bestcase[1,jj] <- max(as.data.frame(ode(state.t.bestcase,tvec,SIR.onepatch.births.sero,parms.t.bestcase))$cum.inc)
  state.c.bestcase <- state.init
  parms.c.bestcase <- c(R0=R0,epsilon=epsilon,C=0,N=1,gamma=gamma,gamma.wan=gamma.wan,
                        beta.amp=beta.amp,
                        mu=mu,offset=offset)
  IAR.bestcase[2,jj] <- max(as.data.frame(ode(state.c.bestcase,tvec,SIR.onepatch.births.sero,parms.c.bestcase))$cum.inc)
  ## S.f.bestcase[1,jj] <- S0 - IAR.bestcase[1,jj]
  ## S.f.bestcase[2,jj] <- S0 - IAR.bestcase[2,jj]
  
  state.t.mosquito <- state.init
  parms.t.mosquito <- c(R0=R0,epsilon=epsilon,N=1,gamma=gamma,gamma.wan=gamma.wan,
                        beta.amp=beta.amp,
                        mu=mu,offset=offset)
  IAR.mosquito[1,jj] <- max(as.data.frame(ode(state.t.mosquito,tvec,SIR.onepatch.births.varC.sero,parms.t.mosquito,C.fn=Ct.fn))$cum.inc)
  state.c.mosquito <- state.init
  parms.c.mosquito <- c(R0=R0,epsilon=epsilon,N=1,gamma=gamma,gamma.wan=gamma.wan,
                        beta.amp=beta.amp,
                        mu=mu,offset=offset)
  IAR.mosquito[2,jj] <- max(as.data.frame(ode(state.c.mosquito,tvec,SIR.onepatch.births.varC.sero,parms.c.mosquito,C.fn=Cc.fn))$cum.inc)
  ## S.f.mosquito[1,jj] <- S0 - IAR.mosquito[1,jj]
  ## S.f.mosquito[2,jj] <- S0 - IAR.mosquito[2,jj]

  IAR.human[1,jj] <- IAR.bestcase[1,jj] * rho.tt + IAR.bestcase[2,jj] * rho.tc
  IAR.human[2,jj] <- IAR.bestcase[2,jj] * rho.cc + IAR.bestcase[1,jj] * rho.ct
  
  IAR.hummoz[1,jj] <- IAR.mosquito[1,jj] * rho.tt + IAR.mosquito[2,jj] * rho.tc
  IAR.hummoz[2,jj] <- IAR.mosquito[2,jj] * rho.cc + IAR.mosquito[1,jj] * rho.ct
  ## S.f.human[1,jj] <- S0 - IAR.human[1,jj]
  ## S.f.human[2,jj] <- S0 - IAR.human[2,jj]
  
  state.humsupp <- rep(state.init,each=2)
  parms.humsupp <- list(R0=R0,epsilon=epsilon,C=c(1,0),N=1,gamma=gamma,gamma.wan=gamma.wan,
                        beta.amp=beta.amp,
                        mu=mu,offset=offset,
                        rho.ij=matrix(c(rho.tt,rho.tc,rho.ct,rho.cc),ncol = 2, byrow = T))
  out.temp <- as.data.frame(ode(state.humsupp,tvec,SIR.twopatch.births.sero,parms.humsupp))
  IAR.humsupp[,jj] <- apply(out.temp[names(out.temp)=="cum.inc"],2,max)

  ## S.f.fullmodel[1,jj] <- S0 - IAR.fullmodel[1,jj]
  ## S.f.fullmodel[2,jj] <- S0 - IAR.fullmodel[2,jj]
  state.fullmodel <- rep(state.init,each=2)
  parms.fullmodel <- list(R0=R0,epsilon=epsilon,N=1,gamma=gamma,gamma.wan=gamma.wan,
                        beta.amp=beta.amp,
                        mu=mu,offset=offset,
                        rho.ij=matrix(c(rho.tt,rho.tc,rho.ct,rho.cc),ncol = 2, byrow = T))
  out.temp <- as.data.frame(ode(state.fullmodel,tvec,SIR.twopatch.births.varC.sero,parms.fullmodel,Ct.fn=Ct.fn,Cc.fn=Cc.fn))
  IAR.fullmodel[,jj] <- apply(out.temp[names(out.temp)=="cum.inc"],2,max)
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
save(epsilon.vec,
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
     file = './fig_2_tempSIRvarC.RData')
