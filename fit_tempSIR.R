## load packages
if(!require(pacman)){install.packages("pacman");library(pacman)}
p_load(deSolve,
       BayesianTools,
       data.table)

# set working directory
setwd('~/Documents/awed_trial_modeling/')

# clear existing workspace
rm(list = ls())

# load functions
source('functions_trial_sim.R')
source('functions_immune_history.R')

## load data https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6221224/
popn.yogyakarta <- 417744
start.day <- as.Date("2006-01-01")
data <- fread("./yogyakarta_dengue_incidence.csv")
data[,day:=floor(day+0.5)]
data[,incidence:=incidence/1e5]
data[,monthly.cases:=floor(incidence * popn.yogyakarta+0.5)]
data[,date:=day + start.day]
data[,month:=format(date,"%m")]
data[,year:=format(date,"%y")]
data[,yearmonth:=paste0(year,"_",month)]
plot(incidence~date,data=data,type='b',xlab="Date",ylab="Incidence")

## Parameters
pop.age = read.csv('./pop_by_age_Indonesia.csv')
age.dist = cbind(
  rep(pop.age$AgeGrpStart,each=pop.age$AgeGrpSpan) + (0:4),
  rep(pop.age$PopTotal/5,each=pop.age$AgeGrpSpan))
age.dist[,2] = age.dist[,2] * 1000
life.expectancy <- sum(pop.age[,ncol(pop.age)] * (age.dist[,2] / sum(age.dist[,2])))
FOI = 0.0317767
inv.cross.im <- 0.5
S0 = calc.prop.susc.ci(age.dist = age.dist,
                       FOI = FOI,
                       inv.cross.im = inv.cross.im)
Sf = S0 * exp(-FOI)
R0 = (log(Sf) - log(S0)) / (Sf - S0)

gamma = 1/7#recovery rate
beta.amp = 0.3 #seasonal amplitude in beta (proportional)
tvec <- seq(-1,366*11,1)
I0 <- data$incidence[1]
mu <- 1/(life.expectancy*365.25)
underreporting <- 0.1
state.init <- c(S=S0,I=I0,R=1-S0-I0,cum.inc=0)
parms <- c(R0=R0,epsilon=0,C=0,N=1,gamma=gamma,beta.amp=beta.amp,mu=mu,offset=2.8*365.25/8)
## Run model
out <- as.data.table(ode(state.init,tvec,SIR.onepatch.births,parms))
out[,actual.incidence:=c(NA,diff(cum.inc)) * 365.25/12] # convert to monthly incidence
out[,date:=time + start.day]
out[,month:=format(date,"%m")]
out[,year:=format(date,"%y")]
out[,yearmonth:=paste0(year,"_",month)]
plot(I~date,data=out,type='l')

## Compare to data
out[,incidence:=actual.incidence * underreporting]
plot(incidence~date,data=data,xlab="Date",ylab="Incidence",
     ylim=c(0,max(out$incidence,data$incidence,na.rm=TRUE)))
lines(out$date,out$incidence,col="red")

## Compare to data - numbers of cases per month
out[,daily.infections:=actual.incidence*12/365.25*popn.yogyakarta]
out[,monthly.infections:=sum(daily.infections,na.rm=TRUE),by=yearmonth]
out.monthly <- out[,.(monthly.infections=sum(daily.infections,na.rm=TRUE)),by=yearmonth
                   ][data,on=.(yearmonth=yearmonth)]
plot(monthly.infections * underreporting ~ date,data=out.monthly,
     xlab="Date",ylab="Monthly cases",type="s")
points(monthly.cases~date,data=data,col="red")

## Fit to data
## Vary amplitude, offset, underreporting, and starting prevalence
beta.amp.bounds <- c(0.1,0.9)
offset.bounds <- c(0,365.25/2)
underreporting.bounds <- c(0.01,1)
I0.bounds <- range(data$incidence) * 12 / 365.25 / gamma
par.bounds <- data.table(beta.amp=beta.amp.bounds,
                         offset=offset.bounds,
                         underreporting=underreporting.bounds,
                         I0=I0.bounds)

## Log-likelihood function
NLL <- function(par) {
    beta.amp <- par[["beta.amp"]]
    offset <- par[["offset"]]
    underreporting <- par[["underreporting"]]
    I0 <- par[["I0"]]
    state.init <- c(S=S0,I=I0,R=1-S0-I0,cum.inc=0)
    parms <- c(R0=R0,epsilon=0,C=0,N=1,gamma=gamma,beta.amp=beta.amp,mu=mu,offset=offset)
    out <- as.data.table(ode(state.init,tvec,SIR.onepatch.births,parms))
    out[,actual.incidence:=c(NA,diff(cum.inc))]
    out[,date:=time + start.day]
    out[,month:=format(date,"%m")]
    out[,year:=format(date,"%y")]
    out[,yearmonth:=paste0(year,"_",month)]
    out[,daily.infections:=actual.incidence*popn.yogyakarta]
    out[,monthly.infections:=sum(daily.infections,na.rm=TRUE),by=yearmonth]
    out.monthly <- out[,.(monthly.infections=sum(daily.infections,na.rm=TRUE)),by=yearmonth
                   ][data,on=.(yearmonth=yearmonth)]
    return(-sum(dpois(out.monthly$monthly.cases,
                      out.monthly$monthly.infections * underreporting,
                      log=TRUE)))
}
LL <- function(par){names(par) <- names(par.bounds);-NLL(par)}

## Optim fit
model.fit <- optim(par=colMeans(par.bounds),
                   fn=NLL,
                   method="L-BFGS-B",
                   lower=par.bounds[1],
                   upper=par.bounds[2])
## save(model.fit,file="tempSIR_optimout.RData")

## Check it
load(file="tempSIR_optimout.RData")
beta.amp <- model.fit$par[["beta.amp"]]
offset <- model.fit$par[["offset"]]
underreporting <- model.fit$par[["underreporting"]]
I0 <- model.fit$par[["I0"]]

state.init <- c(S=S0,I=I0,R=1-S0-I0,cum.inc=0)
parms <- c(R0=R0,epsilon=0,C=0,N=1,gamma=gamma,beta.amp=beta.amp,mu=mu,offset=offset)
out <- as.data.table(ode(state.init,tvec,SIR.onepatch.births,parms))
out[,actual.incidence:=c(NA,diff(cum.inc)) * 365.25/12] # convert to monthly incidence
out[,date:=time + start.day]
out[,month:=format(date,"%m")]
out[,year:=format(date,"%y")]
out[,yearmonth:=paste0(year,"_",month)]
out[,incidence:=actual.incidence * underreporting]
out[,daily.infections:=actual.incidence*12/365.25*popn.yogyakarta]
out[,monthly.infections:=sum(daily.infections,na.rm=TRUE),by=yearmonth]
out.monthly <- out[,.(monthly.infections=sum(daily.infections,na.rm=TRUE)),by=yearmonth
                   ][data,on=.(yearmonth=yearmonth)]
pdf("optim_fit_tempSIR.pdf",width=14)
plot(x=out.monthly$date,y=rep(NA,nrow(out.monthly)),
     ylim=c(0,1.1*max(out.monthly$monthly.infections*underreporting,out.monthly$monthly.cases)),
     xlab="Date",ylab="Monthly cases",
     bty="n",las=1,yaxs="i")
polygon(c(out.monthly$date,rev(out.monthly$date)),
        c(rep(0,nrow(out.monthly)),rev(out.monthly$monthly.infections * underreporting)),
        col=adjustcolor("grey",0.3),
        border=FALSE)
points(monthly.cases~date,data=data,col="red",type="b")
legend("topright",legend=c("model","data"),fill=c(adjustcolor("grey",0.3),NA),
       pch=c(NA,1),col=c(NA,"red"),lty=c(NA,1),border=FALSE,bty="n")
dev.off()

## Save ICs for each season
yearly.ics <- out[yday(date)==1,.(year=year(date),S,I,R)]
save(yearly.ics,file="yearly_initial_conditions.RData")




## Check mcmc fit
load(file="tempSIR_mcmcout.RData")
plot(mcmc.out,start=5000)
samples = getSample(mcmc.out,start=5e3)
n.samples <- 100
subsamples <- samples[sample(nrow(samples),n.samples,replace=TRUE),]

beta.amp.vec <- subsamples[,"beta.amp"]
offset.vec <- subsamples[,"offset"]
underreporting.vec <- subsamples[,"underreporting"]
I0.vec <- subsamples[,"I0"]
monthly.inf.mat <- matrix(NA,nrow=length(I0.vec),ncol=nrow(out.monthly))
for (ii in 1:length(I0.vec)){
    beta.amp <- beta.amp.vec[ii]
    offset <- offset.vec[ii]
    underreporting <- underreporting.vec[ii]
    I0 <-  I0.vec[ii]
    state.init <- c(S=S0,I=I0,R=1-S0-I0,cum.inc=0)
    parms <- c(R0=R0,epsilon=0,C=0,N=1,gamma=gamma,beta.amp=beta.amp,mu=mu,offset=offset)
    out <- as.data.table(ode(state.init,tvec,SIR.onepatch.births,parms))
    out[,actual.incidence:=c(NA,diff(cum.inc)) * 365.25/12] # convert to monthly incidence
    out[,date:=time + start.day]
    out[,month:=format(date,"%m")]
    out[,year:=format(date,"%y")]
    out[,yearmonth:=paste0(year,"_",month)]
    out[,incidence:=actual.incidence * underreporting]
    out[,daily.infections:=actual.incidence*12/365.25*popn.yogyakarta]
    out[,monthly.infections:=sum(daily.infections,na.rm=TRUE),by=yearmonth]
    out.monthly <- out[,.(monthly.infections=sum(daily.infections,na.rm=TRUE)),by=yearmonth
                       ][data,on=.(yearmonth=yearmonth)]
    monthly.inf.mat[ii,] <- out.monthly$monthly.infections * underreporting
}

monthly.inf.mean <- colMeans(monthly.inf.mat)
monthly.inf.ci <- apply(monthly.inf.mat,2,function(x)quantile(x,c(0,1)))

plot(out.monthly$date,monthly.inf.mean,
     ylim=c(0,1.1*max(monthly.inf.ci)),
     xlab="Date",ylab="Monthly cases",
     bty="n",las=1,yaxs="i",type='l',lwd=1)
polygon(c(out.monthly$date,rev(out.monthly$date)),
        c(monthly.inf.ci[1,],rev(monthly.inf.ci[2,])),
        col=adjustcolor("green",0.5),
        border=FALSE)
points(monthly.cases~date,data=data,col="red",type="b")
legend("topright",legend=c("model","data"),fill=c(adjustcolor("grey",0.3),NA),
       pch=c(NA,1),col=c(NA,"red"),lty=c(NA,1),border=FALSE,bty="n")
