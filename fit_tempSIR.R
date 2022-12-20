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

## Control parms
fit.model <- FALSE
save.fit <- FALSE
regenerate.plot <- TRUE

## load data https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6221224/
popn.yogyakarta <- 417744
start.day <- as.Date("2006-01-01")
data <- fread("./yogyakarta_dengue_incidence.csv")
data[,day:=floor(day+0.5)]
data[,incidence:=incidence/1e5]
data[,date:=day + start.day]
data[,month:=format(date,"%m")]
data[,year:=format(date,"%y")]
data[,yearmonth:=paste0(year,"_",month)]
data[,average:=mean(incidence),by=month]
data[,monthly.cases:=floor(incidence * popn.yogyakarta+0.5)]
data[,avg.monthly.cases:=floor(average * popn.yogyakarta+0.5)]
plot(incidence~date,data=data,type='b',xlab="Date",ylab="Incidence")
lines(average~date,data=data,col='red',lty=2,lwd=2)

## Parameters
## load("foi_and_r0.RData")
load("foi_yogya.RData")
pop.age = read.csv('./pop_by_age_Indonesia.csv')
age.dist = cbind(
  rep(pop.age$AgeGrpStart,each=pop.age$AgeGrpSpan) + (0:4),
  rep(pop.age$PopTotal/5,each=pop.age$AgeGrpSpan))
age.dist[,2] = age.dist[,2] * 1000
life.expectancy <- pop.age[1,"LifeExpectancyAtBirth"]
FOI = foi.yogya #foi.val
inv.cross.im <- 0.5
valency = calc.prop.susc.ci.byvalency(age.dist = age.dist,
                                      FOI = FOI,
                                      inv.cross.im = inv.cross.im)
S0 <- sum(valency[c(1,3,5,7)] * c(1,.75,.5,.25))
Sf = S0 * exp(-4*FOI)
R0 = (log(Sf) - log(S0)) / (Sf - S0)
## R0 <- R0a.val
## R0 <- R0b.val

gamma = 1/7#recovery rate
gamma.wan <- gamma * (inv.cross.im/365.25) / (gamma - inv.cross.im/365.25) #CI waning rate
prop.I <- gamma.wan/(gamma + gamma.wan) # proportion CI that are in the I compartment
beta.amp = 0.3 #seasonal amplitude in beta (proportional)
tvec <- seq(-366*80,366*11,1)
## I0 <- data$incidence[1]
mu <- 1/(life.expectancy*365.25)
underreporting <- 0.1
state.init <- c(S4=valency[[1]],
                I3=valency[[2]]*prop.I,R3=valency[[2]]*(1-prop.I),S3=valency[[3]],
                I2=valency[[4]]*prop.I,R2=valency[[4]]*(1-prop.I),S2=valency[[5]],
                I1=valency[[6]]*prop.I,R1=valency[[6]]*(1-prop.I),S1=valency[[7]],
                I0=valency[[8]]*prop.I,R0=valency[[8]]*(1-prop.I) + valency[[9]],
                cum.inc=0)
parms <- c(R0=R0,epsilon=0,C=0,N=1,gamma=gamma,gamma.wan=gamma.wan,
           beta.amp=beta.amp,mu=mu,offset=2.8*365.25/8)

## Fit to data
## Vary amplitude, offset, underreporting, and starting prevalence
beta.amp.bounds <- c(0,0.2)
offset.bounds <- c(0,180)
underreporting.bounds <- c(0.01,0.2)
## R0.bounds <- c(R0,2*R0)
## I0.bounds <- range(data$incidence) * 12 / 365.25 / gamma
par.bounds <- data.table(beta.amp=beta.amp.bounds,
                         offset=offset.bounds,
                         underreporting=underreporting.bounds)#,
                         ## R0=R0.bounds)

## Log-likelihood function
NLL <- function(par) {
    ## print(par)
    beta.amp <- par[["beta.amp"]]
    offset <- par[["offset"]]
    underreporting <- par[["underreporting"]]
    ## R0 <- par[["R0"]]
    ## I0 <- par[["I0"]]
    ## state.init <- c(S=S0,I=I0,R=1-S0-I0,cum.inc=0)
    parms <- c(R0=R0,epsilon=0,C=0,N=1,gamma=gamma,gamma.wan=gamma.wan,
               beta.amp=beta.amp,mu=mu,offset=offset)
    out <- as.data.table(ode(state.init,tvec,SIR.onepatch.births.sero,parms))
    out[,actual.incidence:=c(NA,diff(cum.inc))]
    out[,date:=time + start.day]
    out[,month:=format(date,"%m")]
    out[,year:=format(date,"%y")]
    out[,yearmonth:=paste0(year,"_",month)]
    out[,daily.infections:=actual.incidence*popn.yogyakarta]
    out[,monthly.infections:=sum(daily.infections,na.rm=TRUE),by=yearmonth]
    out.monthly <- out[,.(monthly.infections=sum(daily.infections,na.rm=TRUE)),by=yearmonth
                   ][data,on=.(yearmonth=yearmonth)]
    return(-sum(dpois(out.monthly$avg.monthly.cases,
                      out.monthly$monthly.infections * underreporting,
                      log=TRUE)))
}
LL <- function(par){names(par) <- names(par.bounds);-NLL(par)}

## Optim fit
if (fit.model){
    model.fit <- optim(par=colMeans(par.bounds),
                       fn=NLL,
                       method="L-BFGS-B",
                       lower=par.bounds[1],
                       upper=par.bounds[2])
    if(save.fit) save(model.fit,file="tempSIR_optimout.RData")
} else {
    load(file="tempSIR_optimout.RData")
}

## Check it
beta.amp <- model.fit$par[["beta.amp"]]
offset <- model.fit$par[["offset"]]
underreporting <- model.fit$par[["underreporting"]]
## R0 <- model.fit$par[["R0"]]
## I0 <- model.fit$par[["I0"]]

## state.init <- c(S=S0,I=I0,R=1-S0-I0,cum.inc=0)
## parms <- c(R0=R0,epsilon=0,C=0,N=1,gamma=gamma,gamma.wan=gamma.wan,
##            beta.amp=0.063,mu=mu,offset=90)
## underreporting <- 0.058
parms <- c(R0=R0,epsilon=0,C=0,N=1,gamma=gamma,gamma.wan=gamma.wan,
           beta.amp=beta.amp,mu=mu,offset=offset)
out <- as.data.table(ode(state.init,tvec,SIR.onepatch.births.sero,parms))
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

## if (regenerate.plot) jpeg("optim_fit_tempSIR.jpg",width=14,height=7,units="in",res=600)
## plot(x=out.monthly$date,y=rep(NA,nrow(out.monthly)),
##      ylim=c(0,1.1*max(out.monthly$monthly.infections*underreporting,out.monthly$monthly.cases)),
##      xlab="Date",ylab="Monthly cases",
##      bty="n",las=1,yaxs="i")
## polygon(c(out.monthly$date,rev(out.monthly$date)),
##         c(rep(0,nrow(out.monthly)),rev(out.monthly$monthly.infections * underreporting)),
##         col=adjustcolor("grey",0.5),
##         border=FALSE)
## points(avg.monthly.cases~date,data=data,col="red",type="b",lwd=2)
## points(monthly.cases~date,data=data,col=adjustcolor("red",0.3),type="b",lwd=0.5)
## legend("topright",legend=c("model","data - averaged","data - raw"),fill=c(adjustcolor("grey",0.5),rep(NA,2)),
##        pch=c(NA,1,1),col=c(NA,"red",adjustcolor("red",0.3)),lty=c(NA,1,1),border=FALSE,bty="n",lwd=c(NA,2,0.5))
## if (regenerate.plot) dev.off()

## if (regenerate.plot) jpeg("optim_fit_tempSIR_darker.jpg",width=14,height=7,units="in",res=600)
## plot(x=out.monthly$date,y=rep(NA,nrow(out.monthly)),
##      ylim=c(0,1.1*max(out.monthly$monthly.infections*underreporting,out.monthly$monthly.cases)),
##      xlab="Date",ylab="Monthly cases",
##      bty="n",las=1,yaxs="i")
## polygon(c(out.monthly$date,rev(out.monthly$date)),
##         c(rep(0,nrow(out.monthly)),rev(out.monthly$monthly.infections * underreporting)),
##         col=adjustcolor("grey",0.5),
##         border=FALSE)
## points(avg.monthly.cases~date,data=data,col="red",type="b",lwd=2)
## points(monthly.cases~date,data=data,col=adjustcolor("red",0.7),type="b",lwd=0.5)
## legend("topright",legend=c("model","data - averaged","data - raw"),fill=c(adjustcolor("grey",0.5),rep(NA,2)),
##        pch=c(NA,1,1),col=c(NA,"red",adjustcolor("red",0.7)),lty=c(NA,1,1),border=FALSE,bty="n",lwd=c(NA,2,0.5))
## if (regenerate.plot) dev.off()

if (regenerate.plot) jpeg("optim_fit_tempSIR_lines.jpg",width=14,height=7,units="in",res=600)
plot(x=out.monthly$date,y=rep(NA,nrow(out.monthly)),
     ylim=c(0,1.1*max(out.monthly$monthly.infections*underreporting,out.monthly$monthly.cases)),
     xlab="Date",ylab="Monthly cases",
     bty="n",las=1,yaxs="i")
points(monthly.cases~date,data=data,col="darkgray",type="b",pch=1,lwd=1.5)
points(avg.monthly.cases~date,data=data,col="black",type="b",cex=1.5,lwd=1.5)
lines(out.monthly$date,out.monthly$monthly.infections * underreporting,
      col=adjustcolor("red",1),lwd=3)
legend("topright",legend=c("model","data - averaged","data - raw"),fill=NA,#c(adjustcolor("grey",0.5),rep(NA,2)),
       pch=c(NA,1,1),col=c("red","black","gray"),lty=c(1,1,1),border=FALSE,bty="n",
       lwd=c(3,2,2),pt.cex=c(NA,1.5,1))
if (regenerate.plot) dev.off()

if (regenerate.plot) jpeg("optim_fit_tempSIR_poster.jpg",width=15.4,height=7.1,units="in",res=600,pointsize=20)
plot(x=out.monthly$date,y=rep(NA,nrow(out.monthly)),
     ylim=c(0,1.1*max(out.monthly$monthly.infections*underreporting,out.monthly$monthly.cases)),
     xlab="Date",ylab="Monthly cases",
     bty="n",las=1,yaxs="i")
points(monthly.cases~date,data=data,col="darkgray",type="b",pch=1,lwd=1.5)
points(avg.monthly.cases~date,data=data,col="black",type="b",cex=1.5,lwd=1.5)
lines(out.monthly$date,out.monthly$monthly.infections * underreporting,
      col=adjustcolor("red",1),lwd=3)
legend(as.Date("2011-03-01"),280,legend=c("model","data - averaged","data - raw"),fill=NA,#c(adjustcolor("grey",0.5),rep(NA,2)),
       pch=c(NA,1,1),col=c("red","black","gray"),lty=c(1,1,1),border=FALSE,bty="n",
       lwd=c(3,2,2),pt.cex=c(NA,1.5,1))
if (regenerate.plot) dev.off()

## Save ICs for each season
yearly.ics <- out[yday(date)==1,.(year=year(date),S4,I3,R3,S3,I2,R2,S2,I1,R1,S1,I0,R0)]
if (save.fit) save(yearly.ics,file="yearly_initial_conditions.RData")


## Check mcmc fit
## load(file="tempSIR_mcmcout.RData")
## plot(mcmc.out,start=5000)
## samples = getSample(mcmc.out,start=5e3)
## n.samples <- 100
## subsamples <- samples[sample(nrow(samples),n.samples,replace=TRUE),]

## beta.amp.vec <- subsamples[,"beta.amp"]
## offset.vec <- subsamples[,"offset"]
## underreporting.vec <- subsamples[,"underreporting"]
## I0.vec <- subsamples[,"I0"]
## monthly.inf.mat <- matrix(NA,nrow=length(I0.vec),ncol=nrow(out.monthly))
## for (ii in 1:length(I0.vec)){
##     beta.amp <- beta.amp.vec[ii]
##     offset <- offset.vec[ii]
##     underreporting <- underreporting.vec[ii]
##     I0 <-  I0.vec[ii]
##     state.init <- c(S=S0,I=I0,R=1-S0-I0,cum.inc=0)
##     parms <- c(R0=R0,epsilon=0,C=0,N=1,gamma=gamma,beta.amp=beta.amp,mu=mu,offset=offset)
##     out <- as.data.table(ode(state.init,tvec,SIR.onepatch.births,parms))
##     out[,actual.incidence:=c(NA,diff(cum.inc)) * 365.25/12] # convert to monthly incidence
##     out[,date:=time + start.day]
##     out[,month:=format(date,"%m")]
##     out[,year:=format(date,"%y")]
##     out[,yearmonth:=paste0(year,"_",month)]
##     out[,incidence:=actual.incidence * underreporting]
##     out[,daily.infections:=actual.incidence*12/365.25*popn.yogyakarta]
##     out[,monthly.infections:=sum(daily.infections,na.rm=TRUE),by=yearmonth]
##     out.monthly <- out[,.(monthly.infections=sum(daily.infections,na.rm=TRUE)),by=yearmonth
##                        ][data,on=.(yearmonth=yearmonth)]
##     monthly.inf.mat[ii,] <- out.monthly$monthly.infections * underreporting
## }

## monthly.inf.mean <- colMeans(monthly.inf.mat)
## monthly.inf.ci <- apply(monthly.inf.mat,2,function(x)quantile(x,c(0,1)))

## plot(out.monthly$date,monthly.inf.mean,
##      ylim=c(0,1.1*max(monthly.inf.ci)),
##      xlab="Date",ylab="Monthly cases",
##      bty="n",las=1,yaxs="i",type='l',lwd=1)
## polygon(c(out.monthly$date,rev(out.monthly$date)),
##         c(monthly.inf.ci[1,],rev(monthly.inf.ci[2,])),
##         col=adjustcolor("green",0.5),
##         border=FALSE)
## points(monthly.cases~date,data=data,col="red",type="b")
## legend("topright",legend=c("model","data"),fill=c(adjustcolor("grey",0.3),NA),
##        pch=c(NA,1),col=c(NA,"red"),lty=c(NA,1),border=FALSE,bty="n")
