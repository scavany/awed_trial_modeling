## set working directory
setwd('~/Documents/awed_trial_modeling')

## clear existing workspace
rm(list = ls())

## load necessary functions
source('functions_immune_history.R')
source('functions_trial_sim.R')

## load the age distribution 
pop.age = read.csv('./pop_by_age_Indonesia.csv')
age.dist = cbind(
  rep(pop.age$AgeGrpStart,each=pop.age$AgeGrpSpan) + (0:4),
  rep(pop.age$PopTotal/5,each=pop.age$AgeGrpSpan))
age.dist[,2] = age.dist[,2] * 1000

## calculate the life expectancy as the weighted average of the life expectancies across age groups
life.expectancy <- sum(pop.age[,ncol(pop.age)] * (age.dist[,2] / sum(age.dist[,2])))

## Panel A - efficacy in the mosquito model for values of Ct and Cc (at baseline epsilon?)
S0 <- 1
R0 <- 3.5
epsilon <- 0.75
Ct <- seq(0.5,1,0.01)
Cc <- rev(1-Ct)
