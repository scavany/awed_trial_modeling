# set working directory
setwd('~/Documents/awed_trial_modeling/')

# clear existing workspace
rm(list = ls())

# install necessary packages
if(!require(pacman)){install.packages('pacman'); library(pacman)}
p_load(deSolve,data.table)

# load necessary functions
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

## plot cluster wMel frequencies to check make sense
pdf("wMel_frequency.pdf",width=14)
par(mfrow=c(1,2))
plot(wMel.mean~Date,data=Ct.mean,
     col="red",lwd=3,ylim=c(0,1),type='b',
     main="Treatment",xlab="Date",ylab="wMel frequency")
for (cluster in Ct.clusters) {
    lines(wMel~Date,data=Ct[Cluster==cluster],col=adjustcolor("grey",0.8))
}
## lines(seq(min(Ct.mean$Date),max(Ct.mean$Date),"1 day"),
##       Ct.fn(seq(0,length.out=length(seq(min(Ct.mean$Date),max(Ct.mean$Date),"1 day")),by=1)),
##       col="green",lwd=2)
plot(wMel.mean~Date,data=Cc.mean,
     col="red",lwd=3,ylim=c(0,1),type='b',
     main="Control",xlab="Date",ylab="wMel frequency")
for (cluster in Cc.clusters) {
    lines(wMel~Date,data=Cc[Cluster==cluster],col=adjustcolor("grey",0.8))
}
dev.off()

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
R0 = (log(Sf) - log(S0)) / (Sf - S0)

## Other parms
gamma = 1/7 #recovery rate
mu <- 1/(life.expectancy*365.25) #mortality rate
tvec.full <- seq(-1,366*11,1)
day0 <- as.Date("2006-01-01")
dvec <- tvec.full + day0
start.date <- as.Date("2006-01-01"); end.date <- as.Date("2007-12-31")
tvec <- tvec.full[which(dvec >= start.date & dvec <= end.date)]
load("./tempSIR_optimout.RData",verbose=TRUE)
beta.amp <- model.fit$par[["beta.amp"]]
offset <- model.fit$par[["offset"]]
I0 <- model.fit$par[["I0"]]

## Analysis 1 - How does efficacy change with epsilon? ## 
epsilon.vec = seq(0,1,by=0.005)
Nt = Nc = 1
rho.tt = rho.cc = rho_tt_checker(b = 60, d = 1e3)
rho.tc = 1 - rho.tt
rho.ct = 1 - rho.cc
efficacy.epsilon.bestcase = numeric(length(epsilon.vec))
efficacy.epsilon.fullmodel = numeric(length(epsilon.vec))
for(jj in 1:length(epsilon.vec)){
  epsilon = epsilon.vec[jj]
  state.t.bestcase <- c(S=S0,I=I0,R=1-S0-I0,cum.inc=0)
  parms.t.bestcase <- c(R0=R0,epsilon=epsilon,C=1,N=1,gamma=gamma,beta.amp=beta.amp,mu=mu,offset=offset)
  IAR.t.bestcase <- max(as.data.frame(ode(state.t.bestcase,tvec,SIR.onepatch.births,parms.t.bestcase))$cum.inc)
  state.c.bestcase <- c(S=S0,I=I0,R=1-S0-I0,cum.inc=0)
  parms.c.bestcase <- c(R0=R0,epsilon=0,C=0,N=1,gamma=gamma,beta.amp=beta.amp,mu=mu,offset=offset)
  IAR.c.bestcase <- max(as.data.frame(ode(state.c.bestcase,tvec,SIR.onepatch.births,parms.c.bestcase))$cum.inc)
  
  ## S.f.t.bestcase <- 1 - (IAR.t.bestcase + (1 - S0))
  ## S.f.c.bestcase <- 1 - (IAR.c.bestcase + (1 - S0))
  ## efficacy.epsilon.bestcase[jj] = 1 - IAR.t.bestcase / IAR.c.bestcase #ARR
  ## efficacy.epsilon.bestcase[jj] = 1 - ((IAR.t.bestcase / IAR.c.bestcase) * (S.f.c.bestcase / S.f.t.bestcase)) #OR method 1
  efficacy.epsilon.bestcase[jj] = 1 - (IAR.t.bestcase / IAR.c.bestcase) * (1 - IAR.c.bestcase) / (1 - IAR.t.bestcase) #OR method 2
  
  state.fullmodel <- rep(c(S=S0,I=I0,R=1-S0-I0,cum.inc=0),each=2)
  parms.fullmodel <- list(R0=R0,epsilon=epsilon,N=1,gamma=gamma,beta.amp=beta.amp,
                          mu=mu,offset=offset,
                          rho.ij=matrix(c(rho.tt,rho.tc,rho.ct,rho.cc),ncol = 2, byrow = T))
  out.temp <- as.data.frame(ode(state.fullmodel,tvec,SIR.twopatch.births.varC,parms.fullmodel,Ct.fn=Ct.fn,Cc.fn=Cc.fn))
  IAR.fullmodel <- apply(out.temp[names(out.temp)=="cum.inc"],2,max)
  ## S.f.t.fullmodel <- 1 - (IAR.fullmodel[1] + (1-S0))
  ## S.f.c.fullmodel <- 1 - (IAR.fullmodel[2] + (1-S0))
  
  ## efficacy.epsilon.fullmodel[jj] = 1 - IAR.fullmodel[1] / IAR.fullmodel[2] #ARR
  ## efficacy.epsilon.fullmodel[jj] = 1 - ((IAR.fullmodel[1] / IAR.fullmodel[2]) * (S.f.c.fullmodel / S.f.t.fullmodel)) #OR method 1
  efficacy.epsilon.fullmodel[jj] = 1 - (IAR.fullmodel[1] / IAR.fullmodel[2]) * (1 - IAR.fullmodel[2]) / (1 - IAR.fullmodel[1]) #OR method 2
}

plot(epsilon.vec, efficacy.epsilon.bestcase, type = 'l',ylim=c(0,1))
lines(epsilon.vec, efficacy.epsilon.fullmodel, lty = 2, col = 'green')
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

    state.fullmodel <- rep(c(S=S0,I=I0,R=1-S0-I0,cum.inc=0),each=2)
    parms.fullmodel <- list(R0=R0,epsilon=epsilon,N=1,gamma=gamma,beta.amp=beta.amp,
                            mu=mu,offset=offset,
                            rho.ij=matrix(c(rho.tt,rho.tc,rho.ct,rho.cc),ncol = 2, byrow = T))
    out.temp <- as.data.frame(ode(state.fullmodel,tvec,SIR.twopatch.births.varC,parms.fullmodel,Ct.fn=Ct.fn,Cc.fn=Cc.fn))
    IAR <- apply(out.temp[names(out.temp)=="cum.inc"],2,max)
    
    ## S.f.t <- 1 - (IAR[1] + (1-S0))
    ## S.f.c <- 1 - (IAR[2] + (1-S0))
    ## efficacy.rho[ii,jj] = 1 - IAR.t.human / IAR.c.human #ARR
    ## efficacy.rho[ii,jj] = 1 - ((IAR[1] / IAR[2]) * (S.f.c / S.f.t)) #OR method 1
    efficacy.rho[ii,jj] = 1 - (IAR[1] / IAR[2]) * (1 - IAR[2]) / (1 - IAR[1]) #OR method 2
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
    
    state.fullmodel <- rep(c(S=S0,I=I0,R=1-S0-I0,cum.inc=0),each=2)
    parms.fullmodel <- list(R0=R0,epsilon=epsilon,N=1,gamma=gamma,beta.amp=beta.amp,
                            mu=mu,offset=offset,
                            rho.ij=matrix(c(rho.tt,rho.tc,rho.ct,rho.cc),ncol = 2, byrow = T))
    out.temp <- as.data.frame(ode(state.fullmodel,tvec,SIR.twopatch.births.varC,parms.fullmodel,Ct.fn=Ct.fn,Cc.fn=Cc.fn))
    IAR <- apply(out.temp[names(out.temp)=="cum.inc"],2,max)
    ## S.f.t <- 1 - (IAR[1] + (1-S0))
    ## S.f.c <- 1 - (IAR[2] + (1-S0))
    ## efficacy.implied[ii,jj] = 1 - IAR[1] / IAR[2] #ARR
    ## efficacy.implied[ii,jj] = 1 - ((IAR[1] / IAR[2]) * (S.f.c / S.f.t)) #OR method 1
    efficacy.implied[ii,jj] = 1 - (IAR[1] / IAR[2]) * (1 - IAR[2]) / (1 - IAR[1]) #OR method 2
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
epsilon.b.98 = c(
  epsilon.implied.vec[which.min(abs(
    efficacy.implied[which.min(
      abs(rho.implied.vec-rho_tt_checker(60,1e3))),] -
      efficacy.trial[1]))],
  epsilon.implied.vec[which.min(abs(
    efficacy.implied[which.min(
      abs(rho.implied.vec-rho_tt_checker(60,1e3))),] -
      efficacy.trial[2]))],
  epsilon.implied.vec[which.min(abs(
    efficacy.implied[which.min(
      abs(rho.implied.vec-rho_tt_checker(60,1e3))),] -
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
b.delta.vec = 60
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
    
    state.fullmodel <- rep(c(S=S0,I=I0,R=1-S0-I0,cum.inc=0),each=2)
    parms.fullmodel <- list(R0=R0,epsilon=epsilon,N=1,gamma=gamma,beta.amp=beta.amp,
                            mu=mu,offset=offset,
                            rho.ij=matrix(c(rho.tt,rho.tc,rho.ct,rho.cc),ncol = 2, byrow = T))
    out.temp <- as.data.frame(ode(state.fullmodel,tvec,SIR.twopatch.births.varC,parms.fullmodel,Ct.fn=Ct.fn,Cc.fn=Cc.fn))
    IAR <- apply(out.temp[names(out.temp)=="cum.inc"],2,max)

    ## S.f.t <- 1 - (IAR[1] + (1-S0))
    ## S.f.c <- 1 - (IAR[2] + (1-S0))
    
    ## efficacy.delta[ii,jj] = 1 - (IAR.t.human / IAR.c.human) #ARR
    ## efficacy.delta[ii,jj] = 1 - ((IAR[1] / IAR[2]) * (S.f.c / S.f.t)) #OR method 1
    efficacy.delta[ii,jj] = 1 - (IAR[1] / IAR[2]) * (1 - IAR[2]) / (1 - IAR[1]) #OR method 2
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
     file = './fig_1_tempSIRvarC.RData')

