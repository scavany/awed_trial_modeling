## set working directory
setwd('~/Documents/awed_trial_modeling')

## clear existing workspace
rm(list = ls())

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
sweep.parms <- sobol_design(lower=c(Ct=0.5,Cc=0,FOI=0,epsilon=0,rho.tt=0.5),
                            upper=c(Ct=1,Cc=0.5,FOI=0.06,epsilon=1,rho.tt=1),
                            nseq=nseq)
sweep.out <- data.frame(eff.bestcase=rep(NA,nseq),eff.mosquito=rep(NA,nseq),
                        eff.human=rep(NA,nseq),eff.suppression=rep(NA,nseq))

## Sweep and store output
save.output = FALSE
for (ii in 1:nseq){
    sweep.out[ii,] <- calc_efficacy(age.dist=age.dist,life.expectancy=life.expectancy,
                                    population_structure="exponential",
                                    Ct=sweep.parms[ii,"Ct"],Cc=sweep.parms[ii,"Cc"],
                                    FOI=sweep.parms[ii,"FOI"],epsilon=sweep.parms[ii,"epsilon"],
                                    rho.tt=sweep.parms[ii,"rho.tt"])
}
if(save.output) save(sweep.parms,sweep.out,file="./sweep_output.RData")

## Calculate bias outputs
load.output=FALSE
if (load.output) load("./sweep_output.RData",verbose=TRUE)
subsample <- sample(1:nrow(sweep.out),nseq,replace=FALSE) # sample to make more manageable
sweep.parms <- sweep.parms[subsample,]
sweep.out <- sweep.out[subsample,]
sweep.out$bias.mosquito <- 1 - sweep.out$eff.mosquito/sweep.out$eff.bestcase
sweep.out$bias.human <- 1 - sweep.out$eff.human/sweep.out$eff.bestcase
sweep.out$bias.suppression <- 1 - sweep.out$eff.suppression/sweep.out$eff.bestcase

## plot efficacies
pdf("efficacy_bestcase_scatters.pdf")
par(mfrow=c(2,2),oma=c(0,0,0,2),mar=0.1+c(5,4,4,0))
plot(sweep.parms[,"Ct"],sweep.out[,"eff.bestcase"],
     ylab="Efficacy - bestcase",xlab="Ct",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"Cc"],sweep.out[,"eff.bestcase"],
     ylab="",xlab="Cc",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"FOI"],sweep.out[,"eff.bestcase"],
     ylab="Efficacy - bestcase",xlab="FOI",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"rho.tt"],sweep.out[,"eff.bestcase"],
     ylab="",xlab="rho_tt",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
dev.off()

pdf("efficacy_mosquito_scatters.pdf")
par(mfrow=c(2,2),oma=c(0,0,0,2),mar=0.1+c(5,4,4,0))
plot(sweep.parms[,"Ct"],sweep.out[,"eff.mosquito"],
     ylab="Efficacy - mosquito",xlab="Ct",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"Cc"],sweep.out[,"eff.mosquito"],
     ylab="",xlab="Cc",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"FOI"],sweep.out[,"eff.mosquito"],
     ylab="Efficacy - mosquito",xlab="FOI",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"rho.tt"],sweep.out[,"eff.mosquito"],
     ylab="",xlab="rho_tt",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
dev.off()

pdf("efficacy_human_scatters.pdf")
par(mfrow=c(2,2),oma=c(0,0,0,2),mar=0.1+c(5,4,4,0))
plot(sweep.parms[,"Ct"],sweep.out[,"eff.human"],
     ylab="Efficacy - human",xlab="Ct",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"Cc"],sweep.out[,"eff.human"],
     ylab="",xlab="Cc",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"FOI"],sweep.out[,"eff.human"],
     ylab="Efficacy - human",xlab="FOI",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"rho.tt"],sweep.out[,"eff.human"],
     ylab="",xlab="rho_tt",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
dev.off()

pdf("efficacy_suppression_scatters.pdf")
par(mfrow=c(2,2),oma=c(0,0,0,2),mar=0.1+c(5,4,4,0))
plot(sweep.parms[,"Ct"],sweep.out[,"eff.suppression"],
     ylab="Efficacy - suppression",xlab="Ct",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"Cc"],sweep.out[,"eff.suppression"],
     ylab="",xlab="Cc",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"FOI"],sweep.out[,"eff.suppression"],
     ylab="Efficacy - suppression",xlab="FOI",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"rho.tt"],sweep.out[,"eff.suppression"],
     ylab="",xlab="rho_tt",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
dev.off()

## plot biases
pdf("bias_mosquito_scatters.pdf")
par(mfrow=c(2,2),oma=c(0,0,0,2),mar=0.1+c(5,4,4,0))
plot(sweep.parms[,"Ct"],sweep.out[,"bias.mosquito"],
     ylab="Bias - mosquito",xlab="Ct",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"Cc"],sweep.out[,"bias.mosquito"],
     ylab="",xlab="Cc",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"FOI"],sweep.out[,"bias.mosquito"],
     ylab="Bias - mosquito",xlab="FOI",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"rho.tt"],sweep.out[,"bias.mosquito"],
     ylab="",xlab="rho_tt",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
dev.off()

pdf("bias_human_scatters.pdf")
par(mfrow=c(2,2),oma=c(0,0,0,2),mar=0.1+c(5,4,4,0))
plot(sweep.parms[,"Ct"],sweep.out[,"bias.human"],
     ylab="Bias - human",xlab="Ct",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"Cc"],sweep.out[,"bias.human"],
     ylab="",xlab="Cc",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"FOI"],sweep.out[,"bias.human"],
     ylab="Bias - human",xlab="FOI",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"rho.tt"],sweep.out[,"bias.human"],
     ylab="",xlab="rho_tt",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
dev.off()

pdf("bias_suppression_scatters.pdf")
par(mfrow=c(2,2),oma=c(0,0,0,2),mar=0.1+c(5,4,4,0))
plot(sweep.parms[,"Ct"],sweep.out[,"bias.suppression"],
     ylab="Bias - suppression",xlab="Ct",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"Cc"],sweep.out[,"bias.suppression"],
     ylab="",xlab="Cc",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"FOI"],sweep.out[,"bias.suppression"],
     ylab="Bias - suppression",xlab="FOI",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
plot(sweep.parms[,"rho.tt"],sweep.out[,"bias.suppression"],
     ylab="",xlab="rho_tt",xaxs="i",yaxs="i",bty="n",las=1,pch=20)
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
ggplot(data=prcc.mosquito, aes(y=index, x=est, xmin=lower, xmax=upper)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  scale_y_continuous(breaks=1:nrow(prcc.mosquito), labels=prcc.mosquito$var.names) +
  labs(title='PRCC - mosquito', x='Partial rank correlation coefficient', y = 'Parameter') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_classic()
dev.off()

pdf("prcc_human_forest.pdf")
ggplot(data=prcc.human, aes(y=index, x=est, xmin=lower, xmax=upper)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  scale_y_continuous(breaks=1:nrow(prcc.human), labels=prcc.human$var.names) +
  labs(title='PRCC - human', x='Partial rank correlation coefficient', y = 'Parameter') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_classic()
dev.off()

pdf("prcc_suppression_forest.pdf")
ggplot(data=prcc.suppression, aes(y=index, x=est, xmin=lower, xmax=upper)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  scale_y_continuous(breaks=1:nrow(prcc.suppression), labels=prcc.suppression$var.names) +
  labs(title='PRCC - suppression', x='Partial rank correlation coefficient', y = 'Parameter') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_classic()
dev.off()
