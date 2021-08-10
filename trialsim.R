library(bbmle)
library(deSolve)


# gamma distribution with parameters shape=2.59 and rate=0.0131 

f = function(b,delta){
  b/2/delta*(1-exp(-delta/b))
}
g = function(b,delta,Delta){
  b/2/delta*(exp(-Delta/b)-2*exp(-(delta+Delta)/b)+exp(-(2*delta+Delta)/b))
}
g.checkerboard = function(b,delta,n){
  b/2/delta*(exp(-(n-1)*delta/b)-2*exp(-n*delta/b)+exp(-(n+1)*delta/b))
}

x = read.csv('infector_infectee_pair_distances.csv')
x = subset(x, distance > 0)
displacement = function(c){
  -sum(log(2 * pi * x$distance * c * exp(-c * x$distance)))
}
m = mle(displacement,start=list(c=1/mean(x$distance)))
coef(m)
# 0.005039196


A0 = 1-sum(sapply(1:10,function(nn)2*g.checkerboard(1/(5e-3),1e3,nn)))




FOI = 0.0317767
pop.age = read.csv('WPP2019_PopulationByAgeSex_Medium.csv')
pop.age = subset(pop.age, Time==2019 & Location=='Indonesia')
age.dist = cbind(
  rep(pop.age$AgeGrpStart,each=pop.age$AgeGrpSpan) + (0:4),
  rep(pop.age$PopTotal/5,each=pop.age$AgeGrpSpan))
age.dist[,2] = age.dist[,2] / sum(age.dist[,2])
age.dist = cbind(
  age.dist,
  exp(-FOI * age.dist[,1]),
  c(0, sapply(1:max(age.dist[,1]),function(a)
    (1-exp(-FOI))*exp(-3/4*FOI*(a-1)) +
    sum(sapply(1:a,function(aa)
      exp(-FOI*aa)*(1-exp(-FOI))*exp(-3/4*FOI*(a-aa)))))))
S0 = sum((age.dist[,3]+age.dist[,4])*age.dist[,2])
N = 1

R0 = 3.582622
gamma = 0.1
beta = R0 * gamma



loss.one = function(pi){
  abs(pi - 1 + S0 * exp(-R0 * pi))
}

loss.two = function(pi.t,pi.c){
  abs(pi.t - 1 + S0 * exp(-(
    pi.t * 1/(1+N) * rho.tt * (1 - Ct * epsilon) * R0 +
    pi.c * N/(1+N) * rho.tc * (1 - Cc * epsilon) * R0))) +
  abs(pi.c - 1 + S0 * exp(-(
    pi.t * 1/(1+N) * rho.ct * (1 - Ct * epsilon) * R0 +
    pi.c * N/(1+N) * rho.cc * (1 - Cc * epsilon) * R0)))
}

IAR = optimize(loss.one,interval=c(0,1))$minimum
IAR = optim(c(0.5,0.5),function(par)
  loss(par[1],par[2]),lower=c(0,0),upper=c(1,1),method='L-BFGS-B')$par

VE.true = epsilon = 0.95

R0 = 3.582622
S0 = sum((age.dist[,3]+age.dist[,4])*age.dist[,2])
Ct = 1
Cc = 0
rho.tt = rho.cc = 0.683
rho.tc = 1 - rho.tt
rho.ct = 1 - rho.cc
IAR = optim(c(0.5,0.5),function(par)
  loss.two(par[1],par[2]),lower=c(0,0),upper=c(1,1),method='L-BFGS-B')$par - (1 - S0)
VE.23 = 1 - IAR[1] / IAR[2]

R0 = 3.582622
S0 = sum((age.dist[,3]+age.dist[,4])*age.dist[,2])
Ct = 0.958
Cc = 0.15
rho.tt = rho.cc = 0.683
rho.tc = 1 - rho.tt
rho.ct = 1 - rho.cc
IAR = optim(c(0.5,0.5),function(par)
  loss.two(par[1],par[2]),lower=c(0,0),upper=c(1,1),method='L-BFGS-B')$par - (1 - S0)
VE.123 = 1 - IAR[1] / IAR[2]

R0.orig = 3.582622
S0 = sum((age.dist[,3]+age.dist[,4])*age.dist[,2])
Ct = 0.958
Cc = 0.15
R0 = (1 - Ct * epsilon) * R0.orig
IAR.t = optimize(loss.one,interval=c(0,1))$minimum - (1 - S0)
R0 = (1 - Cc * epsilon) * R0.orig
IAR.c = optimize(loss.one,interval=c(0,1))$minimum - (1 - S0)
VE.1 = 1 - IAR.t / IAR.c

c(VE.true,VE.1,VE.23,VE.123)







model = function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # forces of infection
    FOIt = beta * (1-eff) * rho_tt * It + beta * (1 - rho_tt) * Ic / N
    FOIc = beta * rho_cc * Ic / N + beta * (1-eff) * (1 - rho_cc) * It
    
    # rate of change
    dSt = - FOIt
    dIt = FOIt * St - gamma * It
    dSc = - FOIc
    dIc = FOIc * Sc - gamma * Ic

    # return the rate of change
    list(c(dSt, dIt, dSc, dIc))
  }) # end with(as.list ...
}



df = expand.grid(
  eff = seq(0.1,0.9,0.2))
df$iar_t = NA
df$iar_c = NA

for(ii in 1:nrow(df)){
  parameters = c(
    beta = 6/18,
    gamma = 1/18,
    rho_tt = 1,
    rho_cc = 1,
    eff = df$eff[ii],
    N = 3)
  
  state = c(
    St = 0.2,
    It = 1e-3,
    Sc = 3 * 0.2,
    Ic = 3 * 1e-3)
  
  times = seq(0, 1e3, by = 1)
  
  out = ode(y = state, times = times, func = model, parms = parameters)
  df$iar_t[ii] = out[1,'St'] - out[nrow(out),'St']
  df$iar_c[ii] = out[1,'Sc'] - out[nrow(out),'Sc']
}
