## set working directory
setwd('~/Documents/awed_trial_modeling')

## clear existing workspace
rm(list = ls())

## install necessary packages
if(!require(pacman)){install.packages('pacman'); library(pacman)}
p_load(pomp,epiR,ggplot2,mgcv)

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
Ct.range <- c(0.5,1); Ct.baseline <- 0.958
Cc.range <- c(0,0.5); Cc.baseline <- 0.15
FOI.range <- c(0,0.06); FOI.baseline <- 0.0318
R0.range <- c(2.33/2,2.33*2); R0.baseline = 2.330239
epsilon.range <- c(0,1); epsilon.baseline <- 0.857
rho.tt.range <- c(0.5,1); rho.tt.baseline <- 0.9

sweep.parms <- sobol_design(lower=c(Ct=Ct.range[1],Cc=Cc.range[1],FOI=FOI.range[1],
                                    epsilon=epsilon.range[1],R0=R0.range[1],rho.tt=rho.tt.range[1]),
                            upper=c(Ct=Ct.range[2],Cc=Cc.range[2],FOI=FOI.range[2],
                                    epsilon=epsilon.range[2],R0=R0.range[2],rho.tt=rho.tt.range[2]),
                            nseq=nseq)
sweep.out <- data.frame(eff.bestcase=rep(NA,nseq),eff.mosquito=rep(NA,nseq),
                        eff.human=rep(NA,nseq),eff.suppression=rep(NA,nseq))

## Sweep and store output
generate.output = FALSE
save.output = FALSE
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
load.output=TRUE
if (load.output) load("./sweep_output.RData",verbose=TRUE)
subsample <- sample(1:nrow(sweep.out),1e4,replace=FALSE) # sample to make more manageable
sweep.parms <- sweep.parms[subsample,]
sweep.out <- sweep.out[subsample,]
sweep.out$bias.mosquito <- 1 - sweep.out$eff.mosquito/sweep.out$eff.bestcase
sweep.out$bias.human <- 1 - sweep.out$eff.human/sweep.out$eff.bestcase
sweep.out$bias.suppression <- 1 - sweep.out$eff.suppression/sweep.out$eff.bestcase

## plot efficacies
pdf("efficacy_bestcase_scatters.pdf")
par(mfrow=c(3,2),oma=c(0,0,0,2),mar=0.1+c(5,4,4,0))
plot(sweep.parms[,"Ct"],sweep.out[,"eff.bestcase"],
     ylab="Efficacy - bestcase",xlab="Ct",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"Cc"],sweep.out[,"eff.bestcase"],
     ylab="",xlab="Cc",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"FOI"],sweep.out[,"eff.bestcase"],
     ylab="Efficacy - bestcase",xlab="FOI",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"rho.tt"],sweep.out[,"eff.bestcase"],
     ylab="",xlab="rho_tt",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"R0"],sweep.out[,"eff.bestcase"],
     ylab="Efficacy - bestcase",xlab="R0",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"epsilon"],sweep.out[,"eff.bestcase"],
     ylab="",xlab="epsilon",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
dev.off()

pdf("efficacy_mosquito_scatters.pdf")
par(mfrow=c(3,2),oma=c(0,0,0,2),mar=0.1+c(5,4,4,0))
plot(sweep.parms[,"Ct"],sweep.out[,"eff.mosquito"],
     ylab="Efficacy - mosquito",xlab="Ct",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"Cc"],sweep.out[,"eff.mosquito"],
     ylab="",xlab="Cc",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"FOI"],sweep.out[,"eff.mosquito"],
     ylab="Efficacy - mosquito",xlab="FOI",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"rho.tt"],sweep.out[,"eff.mosquito"],
     ylab="",xlab="rho_tt",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"R0"],sweep.out[,"eff.mosquito"],
     ylab="Efficacy - mosquito",xlab="R0",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"epsilon"],sweep.out[,"eff.mosquito"],
     ylab="",xlab="epsilon",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
dev.off()

pdf("efficacy_human_scatters.pdf")
par(mfrow=c(3,2),oma=c(0,0,0,2),mar=0.1+c(5,4,4,0))
plot(sweep.parms[,"Ct"],sweep.out[,"eff.human"],
     ylab="Efficacy - human",xlab="Ct",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"Cc"],sweep.out[,"eff.human"],
     ylab="",xlab="Cc",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"FOI"],sweep.out[,"eff.human"],
     ylab="Efficacy - human",xlab="FOI",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"rho.tt"],sweep.out[,"eff.human"],
     ylab="",xlab="rho_tt",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"R0"],sweep.out[,"eff.human"],
     ylab="Efficacy - human",xlab="R0",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"epsilon"],sweep.out[,"eff.human"],
     ylab="",xlab="epsilon",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
dev.off()

pdf("efficacy_suppression_scatters.pdf")
par(mfrow=c(3,2),oma=c(0,0,0,2),mar=0.1+c(5,4,4,0))
plot(sweep.parms[,"Ct"],sweep.out[,"eff.suppression"],
     ylab="Efficacy - suppression",xlab="Ct",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"Cc"],sweep.out[,"eff.suppression"],
     ylab="",xlab="Cc",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"FOI"],sweep.out[,"eff.suppression"],
     ylab="Efficacy - suppression",xlab="FOI",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"rho.tt"],sweep.out[,"eff.suppression"],
     ylab="",xlab="rho_tt",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"R0"],sweep.out[,"eff.suppression"],
     ylab="Efficacy - suppression",xlab="R0",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"epsilon"],sweep.out[,"eff.suppression"],
     ylab="",xlab="epsilon",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
dev.off()

## plot biases
## First, get one-at-a-time at baseline
nplot <- 101
Ct.plotrange <- seq(Ct.range[1],Ct.range[2],length.out=nplot)
Cc.plotrange <- seq(Cc.range[1],Cc.range[2],length.out=nplot)
FOI.plotrange <- seq(FOI.range[1],FOI.range[2],length.out=nplot)
rho.tt.plotrange <- seq(rho.tt.range[1],rho.tt.range[2],length.out=nplot)
R0.plotrange <- seq(R0.range[1],R0.range[2],length.out=nplot)
epsilon.plotrange <- seq(epsilon.range[1],epsilon.range[2],length.out=nplot)
efficacy.Ct <- data.frame(eff.bestcase=rep(NA,nplot),eff.mosquito=rep(NA,nplot),
                          eff.human=rep(NA,nplot),eff.suppression=rep(NA,nplot))
efficacy.Cc <- efficacy.Ct; efficacy.rho.tt <- efficacy.Ct; efficacy.FOI <- efficacy.Ct
efficacy.R0 <- efficacy.Ct; efficacy.epsilon <- efficacy.Ct

for (ii in 1:nplot) {
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

## plot scatters with one-at-a-time overlaid
pdf("bias_mosquito_scatters.pdf")
par(mfrow=c(3,2),oma=c(0,0,0,2),mar=0.1+c(5,4,4,0))
plot(sweep.parms[,"Ct"],sweep.out[,"bias.mosquito"],
     ylab="Bias - mosquito",xlab="Ct",xlim=Ct.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(0,1))
lines(Ct.plotrange,
      1 - efficacy.Ct$eff.mosquito / efficacy.Ct$eff.bestcase,
      col="red",lwd=3)
abline(v=Ct.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"Cc"],sweep.out[,"bias.mosquito"],
     ylab="",xlab="Cc",xlim=Cc.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(0,1))
lines(Cc.plotrange,
      1 - efficacy.Cc$eff.mosquito / efficacy.Cc$eff.bestcase,
      col="red",lwd=3)
abline(v=Cc.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"FOI"],sweep.out[,"bias.mosquito"],
     ylab="Bias - mosquito",xlab="FOI",xlim=FOI.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(0,1))
lines(FOI.plotrange,
      1 - efficacy.FOI$eff.mosquito / efficacy.FOI$eff.bestcase,
      col="red",lwd=3)
abline(v=FOI.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"rho.tt"],sweep.out[,"bias.mosquito"],
     ylab="",xlab="rho_tt",xlim=rho.tt.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(0,1))
lines(rho.tt.plotrange,
      1 - efficacy.rho.tt$eff.mosquito / efficacy.rho.tt$eff.bestcase,
      col="red",lwd=3)
abline(v=rho.tt.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"R0"],sweep.out[,"bias.mosquito"],
     ylab="Bias - mosquito",xlab="R0",xlim=R0.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(0,1))
lines(R0.plotrange,
      1 - efficacy.R0$eff.mosquito / efficacy.R0$eff.bestcase,
      col="red",lwd=3)
abline(v=R0.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"epsilon"],sweep.out[,"bias.mosquito"],
     ylab="",xlab="R0",xlim=epsilon.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(0,1))
lines(epsilon.plotrange,
      1 - efficacy.epsilon$eff.mosquito / efficacy.epsilon$eff.bestcase,
      col="red",lwd=3)
abline(v=epsilon.baseline,lty="dashed",col="red",lwd=2)
dev.off()

pdf("bias_human_scatters.pdf")
par(mfrow=c(3,2),oma=c(0,0,0,2),mar=0.1+c(5,4,4,0))
plot(sweep.parms[,"Ct"],sweep.out[,"bias.human"],
     ylab="Bias - human",xlab="Ct",xlim=Ct.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(0,1))
lines(Ct.plotrange,
      1 - efficacy.Ct$eff.human / efficacy.Ct$eff.bestcase,
      col="red",lwd=3)
abline(v=Ct.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"Cc"],sweep.out[,"bias.human"],
     ylab="",xlab="Cc",xlim=Cc.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(0,1))
lines(Cc.plotrange,
      1 - efficacy.Cc$eff.human / efficacy.Cc$eff.bestcase,
      col="red",lwd=3)
abline(v=Cc.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"FOI"],sweep.out[,"bias.human"],
     ylab="Bias - human",xlab="FOI",xlim=FOI.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(0,1))
lines(FOI.plotrange,
      1 - efficacy.FOI$eff.human / efficacy.FOI$eff.bestcase,
      col="red",lwd=3)
abline(v=FOI.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"rho.tt"],sweep.out[,"bias.human"],
     ylab="",xlab="rho_tt",xlim=rho.tt.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(0,1))
lines(rho.tt.plotrange,
      1 - efficacy.rho.tt$eff.human / efficacy.rho.tt$eff.bestcase,
      col="red",lwd=3)
abline(v=rho.tt.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"R0"],sweep.out[,"bias.human"],
     ylab="Bias - human",xlab="R0",xlim=R0.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(0,1))
lines(R0.plotrange,
      1 - efficacy.R0$eff.human / efficacy.R0$eff.bestcase,
      col="red",lwd=3)
abline(v=R0.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"epsilon"],sweep.out[,"bias.human"],
     ylab="",xlab="R0",xlim=epsilon.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(0,1))
lines(epsilon.plotrange,
      1 - efficacy.epsilon$eff.human / efficacy.epsilon$eff.bestcase,
      col="red",lwd=3)
abline(v=epsilon.baseline,lty="dashed",col="red",lwd=2)
dev.off()

pdf("bias_suppression_scatters.pdf")
par(mfrow=c(3,2),oma=c(0,0,0,2),mar=0.1+c(5,4,4,0))
plot(sweep.parms[,"Ct"],sweep.out[,"bias.suppression"],
     ylab="Bias - suppression",xlab="Ct",xlim=Ct.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(0,1))
lines(Ct.plotrange,
      1 - efficacy.Ct$eff.suppression / efficacy.Ct$eff.bestcase,
      col="red",lwd=3)
abline(v=Ct.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"Cc"],sweep.out[,"bias.suppression"],
     ylab="",xlab="Cc",xlim=Cc.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(0,1))
lines(Cc.plotrange,
      1 - efficacy.Cc$eff.suppression / efficacy.Cc$eff.bestcase,
      col="red",lwd=3)
abline(v=Cc.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"FOI"],sweep.out[,"bias.suppression"],
     ylab="Bias - suppression",xlab="FOI",xlim=FOI.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(0,1))
lines(FOI.plotrange,
      1 - efficacy.FOI$eff.suppression / efficacy.FOI$eff.bestcase,
      col="red",lwd=3)
abline(v=FOI.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"rho.tt"],sweep.out[,"bias.suppression"],
     ylab="",xlab="rho_tt",xlim=rho.tt.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(0,1))
lines(rho.tt.plotrange,
      1 - efficacy.rho.tt$eff.suppression / efficacy.rho.tt$eff.bestcase,
      col="red",lwd=3)
abline(v=rho.tt.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"R0"],sweep.out[,"bias.suppression"],
     ylab="Bias - suppression",xlab="R0",xlim=R0.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(0,1))
lines(R0.plotrange,
      1 - efficacy.R0$eff.suppression / efficacy.R0$eff.bestcase,
      col="red",lwd=3)
abline(v=R0.baseline,lty="dashed",col="red",lwd=2)
plot(sweep.parms[,"epsilon"],sweep.out[,"bias.suppression"],
     ylab="",xlab="R0",xlim=epsilon.range,yaxs="i",bty="n",las=1,pch=20,ylim=c(0,1))
lines(epsilon.plotrange,
      1 - efficacy.epsilon$eff.suppression / efficacy.epsilon$eff.bestcase,
      col="red",lwd=3)
abline(v=epsilon.baseline,lty="dashed",col="red",lwd=2)
dev.off()


## GAM Plots
pdf("bias_mosquito_gam.pdf")
gam.out <- gam(qlogis(bias.mosquito)~s(Ct)+s(Cc)+s(FOI)+s(rho.tt)+s(epsilon)+s(R0),
               data=cbind(sweep.parms,bias.mosquito=sweep.out$bias.mosquito))
plot(gam.out, pages=1, rug=FALSE, lwd=3, residuals=TRUE, shift=coef(gam.out)[1],
     trans = plogis)
dev.off()

pdf("bias_human_gam.pdf")
gam.out <- gam(qlogis(bias.human)~s(Ct)+s(Cc)+s(FOI)+s(rho.tt)+s(epsilon)+s(R0),
               data=cbind(sweep.parms,bias.human=sweep.out$bias.human))
plot(gam.out, pages=1, rug=FALSE, lwd=3, residuals=TRUE, shift=coef(gam.out)[1],
     trans = plogis)
dev.off()

pdf("bias_suppression_gam.pdf")
gam.out <- gam(qlogis(bias.suppression)~s(Ct)+s(Cc)+s(FOI)+s(rho.tt)+s(epsilon)+s(R0),
               data=cbind(sweep.parms,bias.suppression=sweep.out$bias.suppression))
plot(gam.out, pages=1, rug=FALSE, lwd=3, residuals=TRUE, shift=coef(gam.out)[1],
     trans = plogis)
dev.off()

## Do PRCC
prcc.mosquito <- cbind(var.names=names(sweep.parms),index=1:ncol(sweep.parms),
                       epi.prcc(cbind(sweep.parms,sweep.out$bias.mosquito)))
prcc.human <- cbind(var.names=names(sweep.parms),index=1:ncol(sweep.parms),
                    epi.prcc(cbind(sweep.parms,sweep.out$bias.human)))
prcc.suppression <- cbind(var.names=names(sweep.parms),index=1:ncol(sweep.parms),
                          epi.prcc(cbind(sweep.parms,sweep.out$bias.suppression)))

## Plot PRCC
pdf("prcc_mosquito_forest.pdf")
ggplot(data=prcc.mosquito, aes(y=index, x=est, xmin=est, xmax=est)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  scale_y_continuous(breaks=1:nrow(prcc.mosquito), labels=prcc.mosquito$var.names) +
  labs(title='PRCC - mosquito', x='Partial rank correlation coefficient', y = 'Parameter') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_classic()
dev.off()

pdf("prcc_human_forest.pdf")
ggplot(data=prcc.human, aes(y=index, x=est, xmin=est, xmax=est)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  scale_y_continuous(breaks=1:nrow(prcc.human), labels=prcc.human$var.names) +
  labs(title='PRCC - human', x='Partial rank correlation coefficient', y = 'Parameter') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_classic()
dev.off()

pdf("prcc_suppression_forest.pdf")
ggplot(data=prcc.suppression, aes(y=index, x=est, xmin=est, xmax=est)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  scale_y_continuous(breaks=1:nrow(prcc.suppression), labels=prcc.suppression$var.names) +
  labs(title='PRCC - suppression', x='Partial rank correlation coefficient', y = 'Parameter') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_classic()
dev.off()
