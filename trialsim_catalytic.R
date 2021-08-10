
# efficacy estimates from the trial
efficacy.trial = c(0.653,0.771,0.849)

# specify forces of infection
# https://stm.sciencemag.org/content/12/528/eaax4144
FOI = 0.0317767



# ideal catalytic
FOI.c = FOI
FOI.t = function(epsilon){(1-epsilon)*FOI}

# get epsilons that would have resulted in observed efficacies
PE = function(epsilon){
  1 - (1-exp(-2*FOI.t(epsilon))) / (1-exp(-2*FOI.c))
}
epsilon.trial = numeric(length(efficacy.trial))
epsilon.trial[1] =
  optimize(f=function(par)abs(PE(par)-efficacy.trial[1]),interval=c(0,1))$minimum
epsilon.trial[2] =
  optimize(f=function(par)abs(PE(par)-efficacy.trial[2]),interval=c(0,1))$minimum
epsilon.trial[3] =
  optimize(f=function(par)abs(PE(par)-efficacy.trial[3]),interval=c(0,1))$minimum
epsilon.trial
# 0.6601513 0.7765421 0.8530148



# mosquito contamination
Ct = 0.958
Cc = 0.15
efficacy.trial.mosquito =
  1 - (1-exp(-2*FOI*(1-Ct*epsilon.trial))) / (1-exp(-2*FOI*(1-Cc*epsilon.trial)))
# get epsilons that would have resulted in observed efficacies
PE = function(epsilon){
  1 - (1-exp(-2*FOI*(1-Ct*epsilon))) / (1-exp(-2*FOI*(1-Cc*epsilon)))
}
epsilon.trial.mosquito = numeric(length(efficacy.trial))
epsilon.trial.mosquito[1] =
  optimize(f=function(par)abs(PE(par)-efficacy.trial[1]),interval=c(0,1))$minimum
epsilon.trial.mosquito[2] =
  optimize(f=function(par)abs(PE(par)-efficacy.trial[2]),interval=c(0,1))$minimum
epsilon.trial.mosquito[3] =
  optimize(f=function(par)abs(PE(par)-efficacy.trial[3]),interval=c(0,1))$minimum
epsilon.trial.mosquito
# 0.7270582 0.8393322 0.9108813



# human contamination
Ct = 0.958
Cc = 0.15
rho.cc = rho.tt = 0.95
rho.ct = rho.tc = 1 - rho.cc
efficacy.trial.human = 1 -
  (1-exp(-2*FOI*(rho.tt*(1-Ct*epsilon.trial)+rho.tc*(1-Cc*epsilon.trial)))) /
  (1-exp(-2*FOI*(rho.cc*(1-Cc*epsilon.trial)+rho.ct*(1-Ct*epsilon.trial))))
# get epsilons that would have resulted in observed efficacies
PE = function(epsilon){
  1 - (1-exp(-2*FOI*(rho.tt*(1-Ct*epsilon)+rho.tc*(1-Cc*epsilon)))) /
      (1-exp(-2*FOI*(rho.cc*(1-Cc*epsilon)+rho.ct*(1-Ct*epsilon))))
}
epsilon.trial.human = numeric(length(efficacy.trial))
epsilon.trial.human[1] =
  optimize(f=function(par)abs(PE(par)-efficacy.trial[1]),interval=c(0,1))$minimum
epsilon.trial.human[2] =
  optimize(f=function(par)abs(PE(par)-efficacy.trial[2]),interval=c(0,1))$minimum
epsilon.trial.human[3] =
  optimize(f=function(par)abs(PE(par)-efficacy.trial[3]),interval=c(0,1))$minimum
epsilon.trial.human
# 0.7729897 0.8865528 0.9582026

