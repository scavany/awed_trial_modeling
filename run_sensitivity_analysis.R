# set working directory
setwd('~/Documents/AWED_Trial/code/')

# clear existing workspace
rm(list = ls())

# install necessary packages
if(!require(pomp)){install.packages('pomp'); library(pomp)}

# load necessary functions
source('functions_sensitivity_analysis.R')

# load the age distribution 
pop.age = read.csv('../data/pop_by_age_Indonesia.csv')
age.dist = cbind(
  rep(pop.age$AgeGrpStart,each=pop.age$AgeGrpSpan) + (0:4),
  rep(pop.age$PopTotal/5,each=pop.age$AgeGrpSpan))
age.dist[,2] = age.dist[,2] * 1000

# calculate the life expectancy as the weighted average of the life expectancies across age groups
life.expectancy <- sum(pop.age[,ncol(pop.age)] * (age.dist[,2] / sum(age.dist[,2])))

