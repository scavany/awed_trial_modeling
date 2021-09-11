## set working directory
setwd('~/Documents/awed_trial_modeling')

## clear existing workspace
rm(list = ls())

## regenerate output, or use old version?
generate.output = FALSE
save.output = FALSE
load.output=TRUE


## install necessary packages
if(!require(pacman)){install.packages('pacman'); library(pacman)}
p_load(pomp,epiR,ggplot2)

## load necessary functions
source('functions_sensitivity_analysis.R')

## load the age distribution 
pop.age = read.csv('./pop_by_age_Indonesia.csv')
age.dist = cbind(
  rep(pop.age$AgeGrpStart,each=pop.age$AgeGrpSpan) + (0:4),
  rep(pop.age$PopTotal/5,each=pop.age$AgeGrpSpan))
age.dist[,2] = age.dist[,2] * 1000

## calculate the life expectancy as the weighted average of the life expectancies across age groups
life.expectancy <- sum(pop.age[,ncol(pop.age)] * (age.dist[,2] / sum(age.dist[,2])))

## Generate sobol sweep of parameters and empty df for output
nseq <- 1e5
Ct.range <- c(0.7,1); Ct.baseline <- 0.950
Cc.range <- c(0,0.4); Cc.baseline <- 0.170
FOI.range <- c(0,0.04); FOI.baseline <- 0.0318
R0.range <- c(3.5/1.25,3.5*1.25); R0.baseline = 3.5#2.330239
epsilon.range <- c(0,1); epsilon.baseline <- 0.849
rho.tt.range <- c(0.6,1); rho.tt.baseline <- 0.887

sweep.parms <- sobol_design(lower=c(Ct=Ct.range[1],Cc=Cc.range[1],FOI=FOI.range[1],
                                    epsilon=epsilon.range[1],R0=R0.range[1],rho.tt=rho.tt.range[1]),
                            upper=c(Ct=Ct.range[2],Cc=Cc.range[2],FOI=FOI.range[2],
                                    epsilon=epsilon.range[2],R0=R0.range[2],rho.tt=rho.tt.range[2]),
                            nseq=nseq)
sweep.out <- data.frame(eff.bestcase=rep(NA,nseq),eff.mosquito=rep(NA,nseq),
                        eff.human=rep(NA,nseq),eff.hummoz=rep(NA,nseq),
                        eff.humsupp=rep(NA,nseq),eff.fullmodel=rep(NA,nseq))

## Sweep and store output
if (generate.output){
    for (ii in 1:nseq){
    sweep.out[ii,] <- calc_efficacy(age.dist=age.dist,life.expectancy=life.expectancy,
                                    population_structure="exponential",
                                    Ct=sweep.parms[ii,"Ct"],Cc=sweep.parms[ii,"Cc"],
                                    FOI=sweep.parms[ii,"FOI"],epsilon=sweep.parms[ii,"epsilon"],
                                    R0=sweep.parms[ii,"R0"],rho.tt=sweep.parms[ii,"rho.tt"])
    }
}

if(save.output) save(sweep.parms,sweep.out,file="./sweep_output.RData")

## Calculate bias outputs
if (load.output) load("./sweep_output.RData",verbose=TRUE)
subsample <- sample(1:nrow(sweep.out),1e4,replace=FALSE) # sample to make more manageable
sweep.parms <- sweep.parms[subsample,]
sweep.out <- sweep.out[subsample,]
sweep.out$bias.mosquito <-  (rowSums(sweep.out[,c(2,4,6)]) - rowSums(sweep.out[,c(1,3,5)])) / 3
sweep.out$bias.human <- (rowSums(sweep.out[,c(3,3,4)]) - rowSums(sweep.out[,c(1,1,2)])) / 3
sweep.out$bias.suppression <- (rowSums(sweep.out[,c(5,6,6)]) - rowSums(sweep.out[,c(3,4,4)])) / 3

## plot biases
## First, get one-at-a-time at baseline
nplot <- 101
Ct.plotrange <- seq(Ct.range[1],Ct.range[2],length.out=nplot)
Cc.plotrange <- seq(Cc.range[1],Cc.range[2],length.out=nplot)
FOI.plotrange <- seq(FOI.range[1]+1e-10,FOI.range[2],length.out=nplot)
rho.tt.plotrange <- seq(rho.tt.range[1],rho.tt.range[2],length.out=nplot)
R0.plotrange <- seq(R0.range[1],R0.range[2],length.out=nplot)
epsilon.plotrange <- seq(epsilon.range[1],epsilon.range[2],length.out=nplot)
efficacy.Ct <- data.frame(eff.bestcase=rep(NA,nplot),eff.mosquito=rep(NA,nplot),
                          eff.human=rep(NA,nplot),eff.hummoz=rep(NA,nplot),
                          eff.humsupp=rep(NA,nplot),eff.fullmodel=rep(NA,nplot))
efficacy.Cc <- efficacy.Ct; efficacy.rho.tt <- efficacy.Ct; efficacy.FOI <- efficacy.Ct
efficacy.R0 <- efficacy.Ct; efficacy.epsilon <- efficacy.Ct

for (ii in 1:nplot) {
    ## print(ii)
    efficacy.Ct[ii,] <- calc_efficacy(age.dist=age.dist,life.expectancy=life.expectancy,
                                     population_structure="exponential",
                                     Ct=Ct.plotrange[ii],Cc=Cc.baseline,
                                     FOI=FOI.baseline,epsilon=epsilon.baseline,
                                     R0=R0.baseline,rho.tt=rho.tt.baseline)
    efficacy.Cc[ii,] <- calc_efficacy(age.dist=age.dist,life.expectancy=life.expectancy,
                                     population_structure="exponential",
                                     Ct=Ct.baseline,Cc=Cc.plotrange[ii],
                                     FOI=FOI.baseline,epsilon=epsilon.baseline,
                                     R0=R0.baseline,rho.tt=rho.tt.baseline)
    efficacy.FOI[ii,] <- calc_efficacy(age.dist=age.dist,life.expectancy=life.expectancy,
                                      population_structure="exponential",
                                      Ct=Ct.baseline,Cc=Cc.baseline,
                                      FOI=FOI.plotrange[ii],epsilon=epsilon.baseline,
                                      R0=R0.baseline,rho.tt=rho.tt.baseline)
    efficacy.rho.tt[ii,] <- calc_efficacy(age.dist=age.dist,life.expectancy=life.expectancy,
                                         population_structure="exponential",
                                         Ct=Ct.baseline,Cc=Cc.baseline,
                                         FOI=FOI.baseline,epsilon=epsilon.baseline,
                                         R0=R0.baseline,rho.tt=rho.tt.plotrange[ii])
    efficacy.R0[ii,] <- calc_efficacy(age.dist=age.dist,life.expectancy=life.expectancy,
                                     population_structure="exponential",
                                     Ct=Ct.baseline,Cc=Cc.baseline,
                                     FOI=FOI.baseline,epsilon=epsilon.baseline,
                                     R0=R0.plotrange[ii],rho.tt=rho.tt.baseline)
    efficacy.epsilon[ii,] <- calc_efficacy(age.dist=age.dist,life.expectancy=life.expectancy,
                                          population_structure="exponential",
                                          Ct=Ct.baseline,Cc=Cc.baseline,
                                          FOI=FOI.baseline,epsilon=epsilon.plotrange[ii],
                                          R0=R0.baseline,rho.tt=rho.tt.baseline)
}
efficacy.Ct$bias.mosquito <-  (rowSums(efficacy.Ct[,c(2,4,6)]) - rowSums(efficacy.Ct[,c(1,3,5)])) / 3
efficacy.Ct$bias.suppression <- (rowSums(efficacy.Ct[,c(5,6,6)]) - rowSums(efficacy.Ct[,c(3,4,4)])) / 3
efficacy.Ct$bias.human <- (rowSums(efficacy.Ct[,c(3,3,4)]) - rowSums(efficacy.Ct[,c(1,1,2)])) / 3
efficacy.Cc$bias.mosquito <-  (rowSums(efficacy.Cc[,c(2,4,6)]) - rowSums(efficacy.Cc[,c(1,3,5)])) / 3
efficacy.Cc$bias.suppression <- (rowSums(efficacy.Cc[,c(5,6,6)]) - rowSums(efficacy.Cc[,c(3,4,4)])) / 3
efficacy.Cc$bias.human <- (rowSums(efficacy.Cc[,c(3,3,4)]) - rowSums(efficacy.Cc[,c(1,1,2)])) / 3
efficacy.rho.tt$bias.mosquito <-  (rowSums(efficacy.rho.tt[,c(2,4,6)]) - rowSums(efficacy.rho.tt[,c(1,3,5)])) / 3
efficacy.rho.tt$bias.suppression <- (rowSums(efficacy.rho.tt[,c(5,6,6)]) - rowSums(efficacy.rho.tt[,c(3,4,4)])) / 3
efficacy.rho.tt$bias.human <- (rowSums(efficacy.rho.tt[,c(3,3,4)]) - rowSums(efficacy.rho.tt[,c(1,1,2)])) / 3
efficacy.R0$bias.mosquito <-  (rowSums(efficacy.R0[,c(2,4,6)]) - rowSums(efficacy.R0[,c(1,3,5)])) / 3
efficacy.R0$bias.suppression <- (rowSums(efficacy.R0[,c(5,6,6)]) - rowSums(efficacy.R0[,c(3,4,4)])) / 3
efficacy.R0$bias.human <- (rowSums(efficacy.R0[,c(3,3,4)]) - rowSums(efficacy.R0[,c(1,1,2)])) / 3
efficacy.epsilon$bias.mosquito <-  (rowSums(efficacy.epsilon[,c(2,4,6)]) - rowSums(efficacy.epsilon[,c(1,3,5)])) / 3
efficacy.epsilon$bias.suppression <- (rowSums(efficacy.epsilon[,c(5,6,6)]) - rowSums(efficacy.epsilon[,c(3,4,4)])) / 3
efficacy.epsilon$bias.human <- (rowSums(efficacy.epsilon[,c(3,3,4)]) - rowSums(efficacy.epsilon[,c(1,1,2)])) / 3
efficacy.FOI$bias.mosquito <-  (rowSums(efficacy.FOI[,c(2,4,6)]) - rowSums(efficacy.FOI[,c(1,3,5)])) / 3
efficacy.FOI$bias.suppression <- (rowSums(efficacy.FOI[,c(5,6,6)]) - rowSums(efficacy.FOI[,c(3,4,4)])) / 3
efficacy.FOI$bias.human <- (rowSums(efficacy.FOI[,c(3,3,4)]) - rowSums(efficacy.FOI[,c(1,1,2)])) / 3

## Do PRCC
prcc.mosquito <- cbind(var.names=names(sweep.parms),index=1:ncol(sweep.parms),
                       epi.prcc(cbind(sweep.parms,sweep.out$bias.mosquito)))
prcc.human <- cbind(var.names=names(sweep.parms),index=1:ncol(sweep.parms),
                    epi.prcc(cbind(sweep.parms,sweep.out$bias.human)))
prcc.suppression <- cbind(var.names=names(sweep.parms),index=1:ncol(sweep.parms),
                          epi.prcc(cbind(sweep.parms,sweep.out$bias.suppression)))

## plot scatters with one-at-a-time overlaid
png("bias_mosquito_scatters.png",width=480*1.5,height=480*1.5,pointsize=14)
par(mfrow=c(3,2),oma=c(0,0,0,2),mar=0.1+c(5,4,4,0))
plot(sweep.parms[,"Ct"],sweep.out[,"bias.mosquito"],
     ylab="Bias - mosquito",xlab="Ct",xlim=Ct.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(-1,1))
lines(Ct.plotrange,
      efficacy.Ct$bias.mosquito,
      col="red",lwd=3)
mtext(paste("PRCC =",round(prcc.mosquito$est[prcc.mosquito$var.names == "Ct"],2)),adj=1,cex=2/3)
abline(v=Ct.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"Cc"],sweep.out[,"bias.mosquito"],
     ylab="",xlab="Cc",xlim=Cc.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(-1,1))
lines(Cc.plotrange,
      efficacy.Cc$bias.mosquito,
      col="red",lwd=3)
mtext(paste("PRCC =",round(prcc.mosquito$est[prcc.mosquito$var.names == "Cc"],2)),adj=1,cex=2/3)
abline(v=Cc.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"FOI"],sweep.out[,"bias.mosquito"],
     ylab="Bias - mosquito",xlab="FOI",xlim=FOI.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(-1,1))
lines(FOI.plotrange,
      efficacy.FOI$bias.mosquito,
      col="red",lwd=3)
mtext(paste("PRCC =",round(prcc.mosquito$est[prcc.mosquito$var.names == "FOI"],2)),adj=1,cex=2/3)
abline(v=FOI.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"rho.tt"],sweep.out[,"bias.mosquito"],
     ylab="",xlab="rho_tt",xlim=rho.tt.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(-1,1))
lines(rho.tt.plotrange,
      efficacy.rho.tt$bias.mosquito,
      col="red",lwd=3)
mtext(paste("PRCC =",round(prcc.mosquito$est[prcc.mosquito$var.names == "rho.tt"],2)),adj=1,cex=2/3)
abline(v=rho.tt.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"R0"],sweep.out[,"bias.mosquito"],
     ylab="Bias - mosquito",xlab="R0",xlim=R0.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(-1,1))
lines(R0.plotrange,
      efficacy.R0$bias.mosquito,
      col="red",lwd=3)
mtext(paste("PRCC =",round(prcc.mosquito$est[prcc.mosquito$var.names == "R0"],2)),adj=1,cex=2/3)
abline(v=R0.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"epsilon"],sweep.out[,"bias.mosquito"],
     ylab="",xlab="epsilon",xlim=epsilon.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(-1,1))
lines(epsilon.plotrange,
      efficacy.epsilon$bias.mosquito,
      col="red",lwd=3)
mtext(paste("PRCC =",round(prcc.mosquito$est[prcc.mosquito$var.names == "epsilon"],2)),adj=1,cex=2/3)
abline(v=epsilon.baseline,lty="dashed",col="red",lwd=2)
dev.off()

png("bias_human_scatters.png",width=480*1.5,height=480*1.5,pointsize=14)
par(mfrow=c(3,2),oma=c(0,0,0,2),mar=0.1+c(5,4,4,0))
plot(sweep.parms[,"Ct"],sweep.out[,"bias.human"],
     ylab="Bias - human",xlab="Ct",xlim=Ct.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(-1,1))
lines(Ct.plotrange,
      efficacy.Ct$bias.human,
      col="red",lwd=3)
mtext(paste("PRCC =",round(prcc.human$est[prcc.human$var.names == "Ct"],2)),adj=1,cex=2/3)
abline(v=Ct.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"Cc"],sweep.out[,"bias.human"],
     ylab="",xlab="Cc",xlim=Cc.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(-1,1))
lines(Cc.plotrange,
      efficacy.Cc$bias.human,
      col="red",lwd=3)
mtext(paste("PRCC =",round(prcc.human$est[prcc.human$var.names == "Cc"],2)),adj=1,cex=2/3)
abline(v=Cc.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"FOI"],sweep.out[,"bias.human"],
     ylab="Bias - human",xlab="FOI",xlim=FOI.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(-1,1))
lines(FOI.plotrange,
      efficacy.FOI$bias.human,
      col="red",lwd=3)
mtext(paste("PRCC =",round(prcc.human$est[prcc.human$var.names == "FOI"],2)),adj=1,cex=2/3)
abline(v=FOI.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"rho.tt"],sweep.out[,"bias.human"],
     ylab="",xlab="rho_tt",xlim=rho.tt.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(-1,1))
lines(rho.tt.plotrange,
      efficacy.rho.tt$bias.human,
      col="red",lwd=3)
mtext(paste("PRCC =",round(prcc.human$est[prcc.human$var.names == "rho.tt"],2)),adj=1,cex=2/3)
abline(v=rho.tt.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"R0"],sweep.out[,"bias.human"],
     ylab="Bias - human",xlab="R0",xlim=R0.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(-1,1))
lines(R0.plotrange,
      efficacy.R0$bias.human,
      col="red",lwd=3)
mtext(paste("PRCC =",round(prcc.human$est[prcc.human$var.names == "R0"],2)),adj=1,cex=2/3)
abline(v=R0.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"epsilon"],sweep.out[,"bias.human"],
     ylab="",xlab="epsilon",xlim=epsilon.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(-1,1))
lines(epsilon.plotrange,
      efficacy.epsilon$bias.human,
      col="red",lwd=3)
mtext(paste("PRCC =",round(prcc.human$est[prcc.human$var.names == "epsilon"],2)),adj=1,cex=2/3)
abline(v=epsilon.baseline,lty="dashed",col="red",lwd=2)
dev.off()

png("bias_suppression_scatters.png",width=480*1.5,height=480*1.5,pointsize=14)
par(mfrow=c(3,2),oma=c(0,0,0,2),mar=0.1+c(5,4,4,0))
plot(sweep.parms[,"Ct"],sweep.out[,"bias.suppression"],
     ylab="Bias - coupling",xlab="Ct",xlim=Ct.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(-1,1))
lines(Ct.plotrange,
      efficacy.Ct$bias.suppression,
      col="red",lwd=3)
mtext(paste("PRCC =",round(prcc.suppression$est[prcc.suppression$var.names == "Ct"],2)),adj=1,cex=2/3)
abline(v=Ct.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"Cc"],sweep.out[,"bias.suppression"],
     ylab="",xlab="Cc",xlim=Cc.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(-1,1))
lines(Cc.plotrange,
      efficacy.Cc$bias.suppression,
      col="red",lwd=3)
mtext(paste("PRCC =",round(prcc.suppression$est[prcc.suppression$var.names == "Cc"],2)),adj=1,cex=2/3)
abline(v=Cc.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"FOI"],sweep.out[,"bias.suppression"],
     ylab="Bias - coupling",xlab="FOI",xlim=FOI.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(-1,1))
lines(FOI.plotrange,
      efficacy.FOI$bias.suppression,
      col="red",lwd=3)
mtext(paste("PRCC =",round(prcc.suppression$est[prcc.suppression$var.names == "FOI"],2)),adj=1,cex=2/3)
abline(v=FOI.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"rho.tt"],sweep.out[,"bias.suppression"],
     ylab="",xlab="rho_tt",xlim=rho.tt.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(-1,1))
lines(rho.tt.plotrange,
      efficacy.rho.tt$bias.suppression,
      col="red",lwd=3)
mtext(paste("PRCC =",round(prcc.suppression$est[prcc.suppression$var.names == "rho.tt"],2)),adj=1,cex=2/3)
abline(v=rho.tt.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"R0"],sweep.out[,"bias.suppression"],
     ylab="Bias - coupling",xlab="R0",xlim=R0.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(-1,1))
lines(R0.plotrange,
      efficacy.R0$bias.suppression,
      col="red",lwd=3)
mtext(paste("PRCC =",round(prcc.suppression$est[prcc.suppression$var.names == "R0"],2)),adj=1,cex=2/3)
abline(v=R0.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"epsilon"],sweep.out[,"bias.suppression"],
     ylab="",xlab="epsilon",xlim=epsilon.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(-1,1))
lines(epsilon.plotrange,
      efficacy.epsilon$bias.suppression,
      col="red",lwd=3)
mtext(paste("PRCC =",round(prcc.suppression$est[prcc.suppression$var.names == "epsilon"],2)),adj=1,cex=2/3)
abline(v=epsilon.baseline,lty="dashed",col="red",lwd=2)
dev.off()

### Sensitivity analysis of proportion attributable to transmission COME BACK TO THIS BIT!!!!
png("prop_suppression_scatters.png",width=480*1.5,height=480*1.5)
par(mfrow=c(3,2),oma=c(0,0,0,2),mar=0.1+c(5,4,4,0))
plot(sweep.parms[,"Ct"],
     sweep.out$bias.suppression/(sweep.out$eff.fullmodel - sweep.out$eff.bestcase),
     ylab="Proportion of bias due to suppression",xlab="Ct",xlim=Ct.range,yaxs="i",bty="n",las=1,
     pch=20,ylim=c(-1,1))
lines(Ct.plotrange,
      efficacy.Ct$bias.suppression/(efficacy.Ct$eff.fullmodel - efficacy.Ct$eff.bestcase),
      col="red",lwd=3)
abline(v=Ct.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"Cc"],
     sweep.out$bias.suppression/(sweep.out$eff.fullmodel - sweep.out$eff.bestcase),
     ylab="",xlab="Cc",xlim=Cc.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(-1,1))
lines(Cc.plotrange,
      efficacy.Cc$bias.suppression/(efficacy.Cc$eff.fullmodel - efficacy.Cc$eff.bestcase),
      col="red",lwd=3)
abline(v=Cc.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"FOI"],
     sweep.out$bias.suppression/(sweep.out$eff.fullmodel - sweep.out$eff.bestcase),
     ylab="Proportion of bias due to suppression",xlab="FOI",xlim=FOI.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(-1,1))
lines(FOI.plotrange,
      efficacy.FOI$bias.suppression/(efficacy.FOI$eff.fullmodel - efficacy.FOI$eff.bestcase),
      col="red",lwd=3)
abline(v=FOI.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"rho.tt"],
     sweep.out$bias.suppression/(sweep.out$eff.fullmodel - sweep.out$eff.bestcase),
     ylab="",xlab="rho_tt",xlim=rho.tt.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(-1,1))
lines(rho.tt.plotrange,
      efficacy.rho.tt$bias.suppression/(efficacy.rho.tt$eff.fullmodel - efficacy.rho.tt$eff.bestcase),
      col="red",lwd=3)
abline(v=rho.tt.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"R0"],
     sweep.out$bias.suppression/(sweep.out$eff.fullmodel - sweep.out$eff.bestcase),
     ylab="Proportion of bias due to suppression",xlab="R0",xlim=R0.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(-1,1))
lines(R0.plotrange,
      efficacy.R0$bias.suppression/(efficacy.R0$eff.fullmodel - efficacy.R0$eff.bestcase),
      col="red",lwd=3)
abline(v=R0.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"epsilon"],sweep.out$bias.suppression/(sweep.out$eff.fullmodel - sweep.out$eff.bestcase),
     ylab="",xlab="epsilon",xlim=epsilon.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(-1,1))
lines(epsilon.plotrange,
      efficacy.epsilon$bias.suppression/(efficacy.epsilon$eff.fullmodel - efficacy.epsilon$eff.bestcase),
      col="red",lwd=3)
abline(v=epsilon.baseline,lty="dashed",col="red",lwd=2)
dev.off()

## PRCC
prcc.prop.suppression <- cbind(var.names=names(sweep.parms),index=1:ncol(sweep.parms),
                               epi.prcc(cbind(sweep.parms,1 - sweep.out[,"bias.human"]/sweep.out[,"bias.suppression"])))

png("prcc_prop_suppression_forest.png",width=480*1.5,height=480*1.5)
ggplot(data=prcc.prop.suppression, aes(y=index, x=est, xmin=est, xmax=est)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  scale_y_continuous(breaks=1:nrow(prcc.suppression), labels=prcc.suppression$var.names) +
  labs(title='PRCC - suppression', x='Partial rank correlation coefficient', y = 'Parameter') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_classic()
dev.off()
