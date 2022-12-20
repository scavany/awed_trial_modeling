## load packages
if(!require(pacman)){install.packages("pacman");library(pacman)}
p_load(data.table)

# set working directory
setwd('~/Documents/awed_trial_modeling/')

# clear existing workspace
rm(list = ls())

## Load and plot data
data <- fread("./wei_movement.csv")
efficacy.trial = c(0.653,0.771,0.849)
pdf("adjusted_data_comparison.pdf")
plot((1:5)/5-0.1,data$eff_mean,ylim=range(data),xlim=c(0,1),
     cex=1.5,
     xaxt="n",xlab="wMel exposure index",
     ylab="Efficacy (%)")
axis(1,at=(1:5)/5-0.1,labels=paste0("[",data$wexp_low,", ",data$wexp_upp,c(rep(")",4),"]")))
for (ii in 2:nrow(data)){
    lines(rep(ii/5-0.1,2),c(data$eff_low[ii],data$eff_upp[ii]),lwd=1.5)
    lines(c(-0.01,0.01)+ii/5-0.1,rep(data$eff_low[ii],2),lwd=1.5)
    lines(c(-0.01,0.01)+ii/5-0.1,rep(data$eff_upp[ii],2),lwd=1.5)
}
abline(h=efficacy.trial*100,lwd=c(1,2,1),lty=c(2,1,2),col="red")
abline(h=0,lty=2)
legend("right",bty="n",legend=c("Baseline","Human-movement adjusted"),
       col=c("red","black"),pch=c(NA,1),lwd=1.5)
dev.off()

## Plot as RR
logistic.fn <- function(x,y.max,y.min,k,x.mid){
    return((y.max-y.min)/(1+exp(-k*(x-x.mid))) + y.min)
}

rr.optim.fn <- function(par) {
    rr <- logistic.fn(seq(0.1,0.9,0.2),par[1],par[2],par[3],par[4])    
    return(sum(abs(rr - 1 + data$eff_mean/100)))
}

rr.optim.out <- optim(c(1.1,0.134,-6,0.5),rr.optim.fn)
rr.optim.out$par;rr.optim.out$convergence;

jpeg("adjusted_data_comparison_full.jpg",height=7, width=7, units="in",res=600)
plot((1:5)/5-0.1,1-data$eff_mean/100,ylim=range(1-data/100),xlim=c(0,1),
     cex=1.5,
     xaxt="n",xlab="wMel exposure index",
     ylab="RR")
axis(1,at=(1:5)/5-0.1,labels=paste0("[",data$wexp_low,", ",data$wexp_upp,c(rep(")",4),"]")))
for (ii in 2:nrow(data)){
    lines(rep(ii/5-0.1,2),1-c(data$eff_low[ii],data$eff_upp[ii])/100,lwd=1.5)
    lines(c(-0.01,0.01)+ii/5-0.1,1-rep(data$eff_low[ii],2)/100,lwd=1.5)
    lines(c(-0.01,0.01)+ii/5-0.1,1-rep(data$eff_upp[ii],2)/100,lwd=1.5)
}
abline(h=1-efficacy.trial,lwd=c(1,2,1),lty=c(2,1,2),col="red")
abline(h=1,lty=2)
legend("right",bty="n",legend=c("Baseline","Human-movement adjusted"),
       col=c("red","black"),pch=c(NA,1),lwd=1.5)
lines(seq(0,1,0.01),
      logistic.fn(seq(0,1,0.01),
                  rr.optim.out$par[1],rr.optim.out$par[2],rr.optim.out$par[3],rr.optim.out$par[4]))
dev.off()

## Calculate RR and efficacy between 0 and 1
rr.full <- logistic.fn(1,rr.optim.out$par[1],rr.optim.out$par[2],rr.optim.out$par[3],rr.optim.out$par[4])/logistic.fn(0,rr.optim.out$par[1],rr.optim.out$par[2],rr.optim.out$par[3],rr.optim.out$par[4])
eff.full <- (1-rr.full)*100

save(eff.full,file="eff_full.RData")
