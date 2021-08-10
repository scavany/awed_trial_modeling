# loss functions for estimating IAR in one or two patch scenarios
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3506030/
loss.one = function(pi){
  abs(pi - 1 + S0 * exp(-R0 * pi))
}
loss.two = function(pi.t,pi.c){
  abs(pi.t - 1 + S0 * exp(-(
    pi.t * Nt/(Nt+Nc) * rho.tt * (1 - Ct * epsilon) * R0 +
    pi.c * Nc/(Nt+Nc) * rho.tc * (1 - Cc * epsilon) * R0))) ^ 2 +
  abs(pi.c - 1 + S0 * exp(-(
    pi.t * Nt/(Nt+Nc) * rho.ct * (1 - Ct * epsilon) * R0 +
    pi.c * Nc/(1+Nc) * rho.cc * (1 - Cc * epsilon) * R0))) ^ 2
}

# calculate susceptible population with Indonesia demography and FOI
# https://stm.sciencemag.org/content/12/528/eaax4144
FOI = 0.0317767
# https://population.un.org/wpp/Download/Standard/Population/
# pop.age = read.csv('WPP2019_PopulationByAgeSex_Medium.csv')
# pop.age = subset(pop.age, Time==2019 & Location=='Indonesia')
# pop.age = write.csv(pop.age,file='pop_by_age_Indonesia.csv')
pop.age = read.csv('pop_by_age_Indonesia.csv')
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



# look at relationship between epsilon and efficacy
R0 = 3.582622
Ct = 0.958
Cc = 0.15
epsilon.vec = seq(0,1,by=0.005)
S0 = sum((age.dist[,3]+age.dist[,4])*age.dist[,2])
Nt = Nc = 1
rho.tt = rho.cc = 1
rho.tc = 1 - rho.tt
rho.ct = 1 - rho.cc
efficacy = numeric(length(epsilon.vec))
for(jj in 1:length(epsilon.vec)){
  epsilon = epsilon.vec[jj]
  IAR = optim(c(0.9,0.9),function(par)
    loss.two(par[1],par[2]),lower=c(0,0),upper=c(1,1),method='L-BFGS-B',
    control = list(reltol=1e-12))$par - (1 - S0)
  efficacy[jj] = 1 - IAR[1] / IAR[2]
}

epsilon.fn = approxfun(efficacy,epsilon.vec)
efficacy.trial = c(0.653,0.771,0.849)
epsilon.trial = epsilon.fn(efficacy.trial)

plot(
  epsilon.vec,efficacy,
  las=1,type='l',xaxs='i',yaxs='i',xlim=c(0,1),ylim=c(0,1),
  xlab=expression('Reduction in '*R[0]*' ('*epsilon*')'),
  ylab='Efficacy observed in trial')
segments(
  c(1,1,1),efficacy.trial,epsilon.trial,efficacy.trial,
  lty=c(2,1,2))
segments(
  epsilon.trial,c(1,1,1),epsilon.trial,efficacy.trial,
  lty=c(2,1,2))
axis(4,at=efficacy.trial,labels=round(efficacy.trial*100)/100,las=1)
axis(3,at=epsilon.trial,labels=c('','',''),las=0)
text(x = epsilon.trial,
     y = 1.16,
     labels = round(epsilon.trial*100)/100,
     xpd = NA,
     srt = 270)


# look at efficacy as a function of rho at a fixed epsilon
R0 = 3.582622
Ct = 0.958
Cc = 0.15
S0 = sum((age.dist[,3]+age.dist[,4])*age.dist[,2])
Nt = Nc = 1
rho.vec = seq(0.5,1,by=0.001)
efficacy = matrix(0,length(rho.vec),3)
for(jj in 1:3){
epsilon = epsilon.trial[jj]
  for(ii in 1:length(rho.vec)){
    rho.tt = rho.cc = rho.vec[ii]
    rho.tc = 1 - rho.tt
    rho.ct = 1 - rho.cc
    IAR = optim(c(0.9,0.9),function(par)
      loss.two(par[1],par[2]),lower=c(0,0),upper=c(1,1),method='L-BFGS-B',
      control = list(reltol=1e-12))$par - (1 - S0)
    efficacy[ii,jj] = 1 - IAR[1] / IAR[2]
  }
}

par(oma=c(0,0,0,1))
plot(
  rho.vec,efficacy[,2],
  type='l',las=1,xaxs='i',yaxs='i',ylim=c(0,0.86),
  xlab=expression('Time at risk in home patch ('*rho*')'),
  ylab='Efficacy observed in trial')
lines(rho.vec,efficacy[,1],lty=2)
lines(rho.vec,efficacy[,3],lty=2)
axis(4,at=efficacy.trial,labels=round(efficacy.trial*100)/100,las=1)



# look at true epsilon implied by inferred epsilon as a function of rho
R0 = 3.582622
Ct = 0.958
Cc = 0.15
epsilon.vec = seq(0.6,1,by=0.005)
S0 = sum((age.dist[,3]+age.dist[,4])*age.dist[,2])
Nt = Nc = 1
rho.vec = seq(0.8,1,by=0.001)
efficacy = matrix(0,length(rho.vec),length(epsilon.vec))
for(jj in 1:length(epsilon.vec)){
  epsilon = epsilon.vec[jj]
  for(ii in 1:length(rho.vec)){
    rho.tt = rho.cc = rho.vec[ii]
    rho.tc = 1 - rho.tt
    rho.ct = 1 - rho.cc
    IAR = optim(c(0.9,0.9),function(par)
      loss.two(par[1],par[2]),lower=c(0,0),upper=c(1,1),method='L-BFGS-B',
      control = list(reltol=1e-12))$par - (1 - S0)
    efficacy[ii,jj] = 1 - IAR[1] / IAR[2]
  }
}

epsilon.fn = approxfun(efficacy[nrow(efficacy),],epsilon.vec)
epsilon.mat = matrix(epsilon.fn(efficacy),nrow(efficacy),ncol(efficacy))
epsilon.trial = epsilon.fn(c(0.653,0.771,0.849))

par(oma=c(0,0.5,0,1))
contour(
  rho.vec,epsilon.vec,epsilon.mat,
  levels=epsilon.trial,
  lty=c(2,1,2),lwd=c(1,2,1),
  las=1,xaxs='i',yaxs='i',drawlabels=F,
  xlab=expression('Time at risk in home patch ('*rho*')'),
  ylab=expression(epsilon*' required for observed efficacy'))
axis(4,at=epsilon.trial,labels=round(epsilon.trial*100)/100,las=1)



# spatial landscapes

# island
A.d1d2 = function(b,d1,d2){
  b / (2 * d1) * (1 - exp(-d1/b) - exp(-d2/b) + exp(-(d1+d2)/b))
}
A.d = function(b,d){
  1 - b / d * (1 - exp(-d/b))
}
rho_tt_island = function(b,d,D){
  A.d(b,d)^2 /
    (A.d(b,d)^2 + 4*A.d(b,d)*A.d1d2(b,d,D) + 4*A.d1d2(b,d,D)^2)
}
rho_tc_island = function(b,d,D){
  (4*A.d(b,d)*A.d1d2(b,d,D) + 4*A.d1d2(b,d,D)^2) /
    (A.d(b,d)^2 + 4*A.d(b,d)*A.d1d2(b,d,D) + 4*A.d1d2(b,d,D)^2)
}
rho_cc_island = function(b,d,D){
  (4*A.d(b,D)^2 + 4*A.d(b,d)*A.d(b,D)) /
    (4*A.d(b,D)^2 + 4*A.d(b,d)*A.d(b,D) + 4*A.d1d2(b,D,d)^2 + 4*A.d1d2(b,D,d)*A.d(b,d))
}
rho_ct_island = function(b,d,D){
  (4*A.d1d2(b,D,d)^2 + 4*A.d1d2(b,D,d)*A.d(b,d)) /
    (4*A.d(b,D)^2 + 4*A.d(b,d)*A.d(b,D) + 4*A.d1d2(b,D,d)^2 + 4*A.d1d2(b,D,d)*A.d(b,d))
}
Nc_island = function(d,D){
  (4 * D ^ 2 + 4 * D * d) / (d ^ 2)
}

# archipelago
A.d1d2d3 = function(b,d1,d2,d3){
  b / (2 * d1) * (exp(-d2/b) - exp(-(d1+d2)/b) - exp(-(d2+d3)/b) + exp(-(d1+d2+d3)/b))
}
B.fn = function(b,d,D){
  A.d(b,d)^2 + 4*A.d1d2d3(b,d,D,d)*A.d(b,d) + 4*A.d1d2d3(b,d,D,d)^2 + 4*A.d1d2d3(b,d,2*D+d,d)*A.d(b,d) +
  8*A.d1d2d3(b,d,2*D+d,d)*A.d1d2d3(b,d,D,d) + 4*A.d1d2d3(b,d,2*D+d,d)^2
}
C.fn = function(b,d,D){
  4*A.d1d2(b,d,D)^2 + 4*A.d1d2(b,d,D)*A.d(b,d) + 8*A.d1d2(b,d,D)*A.d1d2d3(b,d,D,d) +
  8*A.d1d2d3(b,d,D+d,D)*A.d1d2(b,d,D) + 4*A.d1d2d3(b,d,D+d,D)^2 +
  8*A.d1d2d3(b,d,2*D+d,d)*A.d1d2(b,d,D) + 8*A.d1d2d3(b,d,2*D+d,d)^2
}
D.fn = function(b,d,D){
  4*A.d1d2(b,D,d)^2 + 8*A.d1d2d3(b,D,d+D,d)*A.d1d2(b,D,d) + 4*A.d1d2d3(b,D,d+D,d)^2
}
E.fn = function(b,d,D){
  A.d(b,D)^2 +
  4*A.d1d2(b,D,d)*A.d(b,D) + 4*A.d1d2d3(b,D,d,D)*A.d(b,D) + 4*A.d1d2d3(b,D,d,D)^2 +
  8*A.d1d2d3(b,D,d,D)*A.d1d2(b,D,d) + 4*A.d1d2d3(b,D,d+D,d)*A.d(b,D) +
  8*A.d1d2d3(b,D,d+D,d)*A.d1d2d3(b,D,d,D)
}
F.fn = function(b,d,D){
  2*A.d1d2(b,D,d)*A.d(b,d) + 4*A.d1d2(b,D,d)*A.d1d2d3(b,d,D,d) +
  4*A.d1d2(b,D,d)*A.d1d2d3(b,d,2*D+d,d) + 2*A.d1d2d3(b,D,d+D,d)*A.d(b,d) +
  4*A.d1d2d3(b,D,d+D,d)*A.d1d2d3(b,d,D,d) + 4*A.d1d2d3(b,D,d+D,d)*A.d1d2d3(b,d,d+2*D,d)
}
G.fn = function(b,d,D){
  A.d(b,D)*A.d(b,d) + 2*A.d1d2(b,d,D)*A.d(b,D) + 2*A.d1d2d3(b,d,D,d)*A.d(b,D) + 2*A.d1d2d3(b,d,D+d,D)*A.d(b,D) +
  2*A.d1d2d3(b,d,2*D+d,d)*A.d(b,D) + 4*A.d1d2(b,d,D)*A.d1d2(b,D,d) + 4*A.d1d2d3(b,d,D+d,D)*A.d1d2(b,D,d) +
  2*A.d1d2d3(b,D,d,D)*A.d(b,d) + 4*A.d1d2d3(b,D,d,D)*A.d1d2(b,d,D) + 4*A.d1d2d3(b,D,d,D)*A.d1d2d3(b,d,D,d) +
  4*A.d1d2d3(b,D,d,D)*A.d1d2d3(b,d,D+d,D) + 4*A.d1d2d3(b,D,d,D)*A.d1d2d3(b,d,2*D+d,d) +
  4*A.d1d2d3(b,D,d+D,d)*A.d1d2(b,d,D) + 4*A.d1d2d3(b,D,d+D,d)*A.d1d2d3(b,d,D+d,D)
}
rho_tt_archipelago = function(b,d,D){
  B.fn(b,d,D) / (B.fn(b,d,D) + C.fn(b,d,D))
}
rho_tc_archipelago = function(b,d,D){
  C.fn(b,d,D) / (B.fn(b,d,D) + C.fn(b,d,D))
}
rho_ct_archipelago = function(b,d,D){
  (D*D.fn(b,d,D) + 2*d*F.fn(b,d,D)) /
    (D*(D.fn(b,d,D)+E.fn(b,d,D)) + 2*d*(F.fn(b,d,D)+G.fn(b,d,D)))
}
rho_cc_archipelago = function(b,d,D){
  (D*E.fn(b,d,D) + 2*d*G.fn(b,d,D)) /
    (D*(D.fn(b,d,D)+E.fn(b,d,D)) + 2*d*(F.fn(b,d,D)+G.fn(b,d,D)))
}
Nc_archipelago = function(d,D){
  (D^2 + 2*D*d) / (d ^ 2)
}

# checkerboard
H.fn = function(b,d){
  A.d(b,d)^2 + 4*A.d1d2(b,d,d)^2 + 4*A.d1d2d3(b,d,d,d)*A.d(b,d) + 4*A.d1d2d3(b,d,d,d)^2 +
  8*A.d1d2d3(b,d,2*d,d)*A.d1d2(b,d,d) + 4*A.d1d2d3(b,d,2*d,d)^2
}
I.fn = function(b,d){
  4*A.d1d2(b,d,d)*A.d(b,d) + 8*A.d1d2d3(b,d,d,d)*A.d1d2(b,d,d) +
  4*A.d1d2d3(b,d,2*d,d)*A.d(b,d) + 8*A.d1d2d3(b,d,2*d,d)^2
}
rho_tt_checker = rho_cc_checker = function(b,d){
  H.fn(b,d) / (H.fn(b,d) + I.fn(b,d))
}
rho_tc_checker = rho_ct_checker = function(b,d){
  I.fn(b,d) / (H.fn(b,d) + I.fn(b,d))
}
Nc_checker = function(d){
  1
}



# plot things about the scale of movement and time at risk
b.vec = seq(1,1e3,1)
rho.vec = numeric(length(b.vec))
for(ii in 1:length(b.vec)){
  rho.vec[ii] = rho_tt_checker(b.vec[ii],1e3)
}
plot(
  b.vec,rho.vec,type='l',xaxs='i',yaxs='i',
  xlim=c(0,1e3),ylim=c(0.5,1),
  las=1,xlab='Mean distance of time at risk (b)',
  ylab=expression('Time at risk in home patch ('*rho*')'))
b.95 = which.min(abs(rho.vec-0.95))
# b.95 = 25
b.90 = which.min(abs(rho.vec-0.90))
# b.90 = 50



# what would epsilon have been under these different values of b?
rho.vec = seq(0.8,1,by=0.001)
epsilon.b.95 = epsilon.vec[which.min(abs(
  efficacy[which.min(
    abs(rho.vec-rho_tt_checker(b.95,1e3))),] -
  efficacy.trial[2]))]
epsilon.b.90 = epsilon.vec[which.min(abs(
  efficacy[which.min(
    abs(rho.vec-rho_tt_checker(b.90,1e3))),] -
    efficacy.trial[2]))]



# what if the dimensions of the checkerboard were bigger or smaller?
R0 = 3.582622
Ct = 0.958
Cc = 0.15
epsilon.vec = epsilon.b.90
b.vec = b.90
S0 = sum((age.dist[,3]+age.dist[,4])*age.dist[,2])
delta.vec = seq(1e2,1e4,length.out=200)
efficacy = matrix(0,3,length(delta.vec))
for(ii in 1:length(epsilon.vec)){
  epsilon = epsilon.vec[ii]
  for(jj in 1:length(delta.vec)){
    delta = delta.vec[jj]
    Nt = 1
    Nc = Nc_checker(delta)
    rho.tt = rho_tt_checker(b.vec[ii],delta)
    rho.cc = rho_cc_checker(b.vec[ii],delta)
    rho.tc = 1 - rho.tt
    rho.ct = 1 - rho.cc
    IAR = optim(c(0.9,0.9),function(par)
      loss.two(par[1],par[2]),lower=c(0,0),upper=c(1,1),method='L-BFGS-B',
      control = list(reltol=1e-12))$par - (1 - S0)
    efficacy[ii,jj] = 1 - IAR[1] / IAR[2]
  }
}
plot(
  delta.vec,efficacy[1,],type='l',
  ylim=c(0,1),xaxs='i',yaxs='i',las=1,
  xlab=expression('Width of treatment patch in m ('*delta*')'),
  ylab='Efficacy that would be observed')
lines(delta.vec,efficacy[2,])
lines(delta.vec,efficacy[3,])
abline(v=1e3,lty=2)



# what about an equivalently sized island?
R0 = 3.582622
Ct = 0.958
Cc = 0.15
epsilon.vec = c(epsilon.b.95,epsilon.b.90)
b.vec = c(b.95,b.90)
S0 = sum((age.dist[,3]+age.dist[,4])*age.dist[,2])
Delta = 5e2
delta.vec = seq(1e2,1e4,length.out=200)
efficacy = matrix(0,3,length(delta.vec))
for(ii in 1:length(epsilon.vec)){
  epsilon = epsilon.vec[ii]
  for(jj in 1:length(delta.vec)){
    delta = delta.vec[jj]
    Nt = 1
    Nc = Nc_island(delta,Delta)
    rho.tt = rho_tt_island(b.vec[ii],delta,Delta)
    rho.cc = rho_cc_island(b.vec[ii],delta,Delta)
    rho.tc = rho_tc_island(b.vec[ii],delta,Delta)
    rho.ct = rho_ct_island(b.vec[ii],delta,Delta)
    IAR = optim(c(0.9,0.9),function(par)
      loss.two(par[1],par[2]),lower=c(0,0),upper=c(1,1),method='L-BFGS-B',
      control = list(reltol=1e-12))$par - (1 - S0)
    efficacy[ii,jj] = 1 - IAR[1] / IAR[2]
  }
}
plot(
  delta.vec,efficacy[1,],type='l',
  ylim=c(0,1),xaxs='i',yaxs='i',las=1,
  xlab=expression('Width of intervention patch in m ('*delta*')'),
  ylab='Efficacy that would be observed')
# lines(delta.vec,efficacy[2,],lty=2)
lines(delta.vec,efficacy[3,])
abline(v=1e3,lty=2)
# legend(
#   'bottom',legend=c('25 m','50 m'),lty=c(1,2),
#   title='Mean dist. (b)',bty='n',horiz=T)



# what about an equivalently sized archipelago?
R0 = 3.582622
Ct = 0.958
Cc = 0.15
epsilon.vec = c(epsilon.b.95,epsilon.b.90)
b.vec = c(b.95,b.90)
S0 = sum((age.dist[,3]+age.dist[,4])*age.dist[,2])
Delta = 1e3
delta.vec = seq(1e1,1e3,length.out=200)
efficacy = matrix(0,3,length(delta.vec))
for(ii in 1:length(epsilon.vec)){
  epsilon = epsilon.vec[ii]
  for(jj in 1:length(delta.vec)){
    delta = delta.vec[jj]
    Nt = 1
    Nc = Nc_archipelago(delta,Delta)
    rho.tt = rho_tt_archipelago(b.vec[ii],delta,Delta)
    rho.cc = rho_cc_archipelago(b.vec[ii],delta,Delta)
    rho.tc = rho_tc_archipelago(b.vec[ii],delta,Delta)
    rho.ct = rho_ct_archipelago(b.vec[ii],delta,Delta)
    IAR = optim(c(0.9,0.9),function(par)
      loss.two(par[1],par[2]),lower=c(0,0),upper=c(1,1),method='L-BFGS-B',
      control = list(reltol=1e-12))$par - (1 - S0)
    efficacy[ii,jj] = 1 - IAR[1] / IAR[2]
  }
}
plot(
  delta.vec,efficacy[1,],type='l',
  ylim=c(0,1),xaxs='i',yaxs='i',las=1,
  xlab=expression('Width of intervention patch in m ('*delta*')'),
  ylab='Efficacy that would be observed')
# lines(delta.vec,efficacy[2,],lty=2)
# legend(
#   'bottom',legend=c('25 m','50 m'),lty=c(1,2),
#   title='Mean dist. (b)',bty='n',horiz=T)



# get FOI for catalytic model equivalent to IAR under SIR model
R0 = 3.582622
Ct = 1
Cc = 0
S0 = sum((age.dist[,3]+age.dist[,4])*age.dist[,2])
Nt = 1
Nc = Nc_checker(delta)
rho.tt = 1
rho.cc = 1
rho.tc = 1 - rho.tt
rho.ct = 1 - rho.cc
epsilon = 0
IAR.max = (optim(c(0.9,0.9),function(par)
  loss.two(par[1],par[2]),lower=c(0,0),upper=c(1,1),method='L-BFGS-B',
  control = list(reltol=1e-12))$par - (1 - S0))[1]
FOI.max = -log(1-IAR.max)



# look at effect on IAR
R0 = 3.582622
Ct = 0.958
Cc = 0.15
epsilon.vec = seq(0,1,by=0.005)
S0 = sum((age.dist[,3]+age.dist[,4])*age.dist[,2])
Nt = 1
Nc = Nc_checker(delta)
rho.tt = rho_tt_checker(b.vec[ii],delta)
rho.cc = rho_cc_checker(b.vec[ii],delta)
rho.tc = 1 - rho.tt
rho.ct = 1 - rho.cc
IAR.mosquito = IAR.human = IAR = matrix(0,2,length(epsilon.vec))
for(jj in 1:length(epsilon.vec)){
  epsilon = epsilon.vec[jj]
  IAR.mosquito[,jj] = c(
    1-exp(-FOI.max*(1-Ct*epsilon)),
    1-exp(-FOI.max*(1-Cc*epsilon)))
  IAR.human[,jj] = c(
    1-exp(-FOI.max*(rho.tt*(1-Ct*epsilon)+rho.tc*(1-Cc*epsilon))),
    1-exp(-FOI.max*(rho.cc*(1-Cc*epsilon)+rho.ct*(1-Ct*epsilon))))
  IAR[,jj] = optim(c(0.9,0.9),function(par)
    loss.two(par[1],par[2]),lower=c(0,0),upper=c(1,1),method='L-BFGS-B',
    control = list(reltol=1e-12))$par - (1 - S0)
}
matplot(t(IAR.mosquito),type='l',lty=1,ylim=c(0,1))
matplot(t(IAR.human),type='l',lty=1,ylim=c(0,1))
matplot(t(IAR),type='l',lty=1,ylim=c(0,1))
matplot(
  epsilon.vec,cbind(
  IAR.mosquito[2,],
  IAR.human[2,],
  IAR[2,]),type='l',lty=1,ylim=c(0,1),
  xlab=expression(epsilon),
  ylab='IAR in control arm')
legend('topright',legend=1:3,col=1:3,lty=1)
matplot(
  epsilon.vec,cbind(
    IAR.mosquito[1,],
    IAR.human[1,],
    IAR[1,]),type='l',lty=1,ylim=c(0,1),
  xlab=expression(epsilon),
  ylab='IAR in treatment arm')
legend('topright',legend=1:3,col=1:3,lty=1)
matplot(
  epsilon.vec,cbind(
    1-IAR.mosquito[1,]/IAR.mosquito[2,],
    1-IAR.human[1,]/IAR.human[2,],
    1-IAR[1,]/IAR[2,]),type='l',lty=1,ylim=c(0,1),
  xlab=expression(epsilon),
  ylab='Efficacy')
legend('bottomright',legend=1:3,col=1:3,lty=1)


