# loss functions for estimating IAR in one or two patch scenarios
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3506030/
## Miller versions
## loss.one = function(pi, S0, R0){
##   abs(pi - 1 + S0 * exp(-R0 * pi))
## }
## loss.two = function(pi.t, pi.c, rho.tt, rho.tc, rho.cc, rho.ct, Ct, Cc, epsilon, S0, R0){
##   abs(pi.t - 1 + S0 * exp(-(
##     pi.t * rho.tt * (1 - Ct * epsilon) * R0 +
##     pi.c * rho.tc * (1 - Cc * epsilon) * R0))) ^ 2 +
##   abs(pi.c - 1 + S0 * exp(-(
##     pi.t * rho.ct * (1 - Ct * epsilon) * R0 +
##     pi.c * rho.cc * (1 - Cc * epsilon) * R0))) ^ 2
## }
## Quan versions
loss.one = function(pi, S0, R0){
  abs(1 - S0 / pi + S0 * exp(-R0 * pi) / pi)
}
loss.two = function(pi.t, pi.c, rho.tt, rho.tc, rho.cc, rho.ct, Ct, Cc, epsilon, S0, R0){
  abs(1 - S0 / pi.t + S0 / pi.t * exp(-(
    pi.t * rho.tt * (1 - Ct * epsilon) * R0 +
    pi.c * rho.tc * (1 - Cc * epsilon) * R0))) ^ 2 +
  abs(1 - S0 / pi.c + S0 / pi.c * exp(-(
    pi.t * rho.ct * (1 - Ct * epsilon) * R0 +
    pi.c * rho.cc * (1 - Cc * epsilon) * R0))) ^ 2
}

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
    4*A.d1d2d3(b,d,2*d,d)*A.d(b,d) + 8*A.d1d2d3(b,d,2*d,d)*A.d1d2d3(b,d,d,d)
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


## SIR models
## Temperature-dependent beta
beta.varT <- function(time,b.mean,b.amp,offset=2.8*365.25/8) {
    return(b.mean * (1 + b.amp * sin((offset + time) * 2 * pi / 365.25)))
}

## single patch without births
SIR.onepatch <- function(t, state, parms) {
    with(as.list(parms), {

        beta <- beta.varT(t,b.mean=R0*gamma,b.amp=beta.amp)
        S <- state[names(state)=="S"]
        I <- state[names(state)=="I"]
        R <- state[names(state)=="R"]
        dS <- -S * (1 - C * epsilon) * beta * I / N
        dI <- S * (1 - C * epsilon) * beta * I / N - gamma*I
        dR <- gamma*I

        return(list(c(dS, dI, dR)))
    })
}

## two patch model without births
SIR.twopatch <- function(t, state, parms) {
    with(as.list(parms), {

        beta <- beta.varT(t,b.mean=R0*gamma,b.amp=beta.amp)
        S <- state[names(state)=="S"]
        I <- state[names(state)=="I"]
        R <- state[names(state)=="R"]
        dS <- -S*(rho.ij %*% ((1-C*epsilon)*beta*I/N))
        dI <- S*(rho.ij %*% ((1-C*epsilon)*beta*I/N)) - gamma*I
        dR <- gamma*I

        return(list(c(dS, dI, dR)))
    })
}

## single patch
SIR.onepatch.births <- function(t, state, parms) {
    with(as.list(parms), {

        beta <- beta.varT(t,b.mean=R0*gamma,b.amp=beta.amp,offset=offset)
        S <- state[names(state)=="S"]
        I <- state[names(state)=="I"]
        R <- state[names(state)=="R"]
        cum.inc <- state[names(state)=="cum.inc"]
        dS <- -S * (1 - C * epsilon) * beta * I / N + mu * (I + R)
        dI <- S * (1 - C * epsilon) * beta * I / N - gamma*I - mu*I
        dR <- gamma*I - mu*R
        dcum.inc <- S * (1 - C * epsilon) * beta * I / N

        return(list(c(dS, dI, dR, dcum.inc)))
    })
}

## two patch model
SIR.twopatch.births <- function(t, state, parms) {
    with(as.list(parms), {

        beta <- beta.varT(t,b.mean=R0*gamma,b.amp=beta.amp,offset=offset)
        S <- state[names(state)=="S"]
        I <- state[names(state)=="I"]
        R <- state[names(state)=="R"]
        cum.inc <- state[names(state)=="cum.inc"]
        dS <- -S*(rho.ij %*% ((1-C*epsilon)*beta*I/N)) + mu * (I + R)
        dI <- S*(rho.ij %*% ((1-C*epsilon)*beta*I/N)) - gamma*I - mu*I
        dR <- gamma*I - mu*R
        dcum.inc <- S*(rho.ij %*% ((1-C*epsilon)*beta*I/N))

        return(list(c(dS, dI, dR, dcum.inc)))
    })
}

## single patch with variable wMel
SIR.onepatch.births.varC <- function(t, state, parms,C.fn) {
    with(as.list(parms), {

        beta <- beta.varT(t,b.mean=R0*gamma,b.amp=beta.amp,offset=offset)
        C <- C.fn(t)
        S <- state[names(state)=="S"]
        I <- state[names(state)=="I"]
        R <- state[names(state)=="R"]
        cum.inc <- state[names(state)=="cum.inc"]
        dS <- -S * (1 - C * epsilon) * beta * I / N + mu * (I + R)
        dI <- S * (1 - C * epsilon) * beta * I / N - gamma*I - mu*I
        dR <- gamma*I - mu*R
        dcum.inc <- S * (1 - C * epsilon) * beta * I / N

        return(list(c(dS, dI, dR, dcum.inc)))
    })
}

## two patch model with variable wMel freq
SIR.twopatch.births.varC <- function(t, state, parms, Ct.fn, Cc.fn) {
    with(as.list(parms), {

        beta <- beta.varT(t,b.mean=R0*gamma,b.amp=beta.amp,offset=offset)
        C <- c(Ct.fn(t),Cc.fn(t))
        S <- state[names(state)=="S"]
        I <- state[names(state)=="I"]
        R <- state[names(state)=="R"]
        cum.inc <- state[names(state)=="cum.inc"]
        dS <- -S*(rho.ij %*% ((1-C*epsilon)*beta*I/N)) + mu * (I + R)
        dI <- S*(rho.ij %*% ((1-C*epsilon)*beta*I/N)) - gamma*I - mu*I
        dR <- gamma*I - mu*R
        dcum.inc <- S*(rho.ij %*% ((1-C*epsilon)*beta*I/N))

        return(list(c(dS, dI, dR, dcum.inc)))
    })
}

## one patch model with serotypes
SIR.onepatch.births.sero <- function(t, state, parms) {
    with(as.list(parms), {

        beta <- beta.varT(t,b.mean=R0*gamma,b.amp=beta.amp,offset=offset)
        S4 <- state[names(state)=="S4"]
        I3 <- state[names(state)=="I3"]
        R3 <- state[names(state)=="R3"]
        S3 <- state[names(state)=="S3"]
        I2 <- state[names(state)=="I2"]
        R2 <- state[names(state)=="R2"]
        S2 <- state[names(state)=="S2"]
        I1 <- state[names(state)=="I1"]
        R1 <- state[names(state)=="R1"]
        S1 <- state[names(state)=="S1"]
        I0 <- state[names(state)=="I0"]
        R0 <- state[names(state)=="R0"]
        cum.inc <- state[names(state)=="cum.inc"]
        I <- I3 + I2 + I1 + I0 
        dS4 <- mu * (1 - S4) - S4 * (1 - C * epsilon) * beta * I / N * 4/4
        dI3 <- S4 * (1 - C * epsilon) * beta * I / N * 4/4 - (gamma + mu)*I3
        dR3 <- gamma*I3 - (gamma.wan + mu)*R3
        dS3 <- gamma.wan*R3 - S3 * (1 - C * epsilon) * beta * I / N * 3/4 - mu*S3
        dI2 <- S3 * (1 - C * epsilon) * beta * I / N * 3/4 - (gamma + mu)*I2
        dR2 <- gamma*I2 - (gamma.wan + mu)*R2
        dS2 <- gamma.wan*R2 - S2 * (1 - C * epsilon) * beta * I / N * 2/4 - mu*S2 
        dI1 <- S2 * (1 - C * epsilon) * beta * I / N * 2/4 - (gamma + mu)*I1
        dR1 <- gamma*I1 - (gamma.wan + mu)*R1
        dS1 <- gamma.wan*R1 - S1 * (1 - C * epsilon) * beta * I / N * 1/4 - mu*S1 
        dI0 <- S1 * (1 - C * epsilon) * beta * I / N * 1/4 - (gamma + mu)*I0
        dR0 <- gamma*I0 - mu*R0
        dcum.inc <- (S4 * 4/4 + S3 * 3/4 + S2 * 2/4 + S1 * 1/4) * (1 - C * epsilon) * beta * I / N

        return(list(c(dS4, dI3, dR3, dS3, dI2, dR2, dS2, dI1, dR1, dS1, dI0, dR0, dcum.inc)))
    })
}

## two patch model with serotypes 
SIR.twopatch.births.sero <- function(t, state, parms) {
    with(as.list(parms), {

        beta <- beta.varT(t,b.mean=R0*gamma,b.amp=beta.amp,offset=offset)
        S4 <- state[names(state)=="S4"]
        I3 <- state[names(state)=="I3"]
        R3 <- state[names(state)=="R3"]
        S3 <- state[names(state)=="S3"]
        I2 <- state[names(state)=="I2"]
        R2 <- state[names(state)=="R2"]
        S2 <- state[names(state)=="S2"]
        I1 <- state[names(state)=="I1"]
        R1 <- state[names(state)=="R1"]
        S1 <- state[names(state)=="S1"]
        I0 <- state[names(state)=="I0"]
        R0 <- state[names(state)=="R0"]
        cum.inc <- state[names(state)=="cum.inc"]
        I <- I3 + I2 + I1 + I0
        dS4 <-  mu * (1 - S4) - S4 * (rho.ij %*% ((1-C*epsilon)*beta*I/N)) * 4/4
        dI3 <- S4 * (rho.ij %*% ((1-C*epsilon)*beta*I/N)) * 4/4 - (gamma + mu)*I3
        dR3 <- gamma*I3 - (gamma.wan + mu)*R3
        dS3 <- gamma.wan*R3 - S3 * (rho.ij %*% ((1-C*epsilon)*beta*I/N)) * 3/4 - mu*S3
        dI2 <- S3 * (rho.ij %*% ((1-C*epsilon)*beta*I/N)) * 3/4 - (gamma + mu)*I2
        dR2 <- gamma*I2 - (gamma.wan + mu)*R2
        dS2 <- gamma.wan*R2 - S2 * (rho.ij %*% ((1-C*epsilon)*beta*I/N)) * 2/4 - mu*S2 
        dI1 <- S2 * (rho.ij %*% ((1-C*epsilon)*beta*I/N)) * 2/4 - (gamma + mu)*I1
        dR1 <- gamma*I1 - (gamma.wan + mu)*R1
        dS1 <- gamma.wan*R1 - S1 * (rho.ij %*% ((1-C*epsilon)*beta*I/N)) * 1/4 - mu*S1 
        dI0 <- S1 * (rho.ij %*% ((1-C*epsilon)*beta*I/N)) * 1/4 - (gamma + mu)*I0
        dR0 <- gamma*I0 - mu*R0
        dcum.inc <- (S4 * 4/4 + S3 * 3/4 + S2 * 2/4 + S1 * 1/4) * (rho.ij %*% ((1-C*epsilon)*beta*I/N))

        return(list(c(dS4, dI3, dR3, dS3, dI2, dR2, dS2, dI1, dR1, dS1, dI0, dR0, dcum.inc)))
    })
}


## single patch with serotypes and variable wMel
SIR.onepatch.births.varC.sero <- function(t, state, parms,C.fn) {
    with(as.list(parms), {
        beta <- beta.varT(t,b.mean=R0*gamma,b.amp=beta.amp,offset=offset)
        C <- C.fn(t)
        S4 <- state[names(state)=="S4"]
        I3 <- state[names(state)=="I3"]
        R3 <- state[names(state)=="R3"]
        S3 <- state[names(state)=="S3"]
        I2 <- state[names(state)=="I2"]
        R2 <- state[names(state)=="R2"]
        S2 <- state[names(state)=="S2"]
        I1 <- state[names(state)=="I1"]
        R1 <- state[names(state)=="R1"]
        S1 <- state[names(state)=="S1"]
        I0 <- state[names(state)=="I0"]
        R0 <- state[names(state)=="R0"]
        cum.inc <- state[names(state)=="cum.inc"]
        I <- I3 + I2 + I1 + I0 
        dS4 <-  mu * (1 - S4) - S4 * (1 - C * epsilon) * beta * I / N * 4/4
        dI3 <- S4 * (1 - C * epsilon) * beta * I / N * 4/4 - (gamma + mu)*I3
        dR3 <- gamma*I3 - (gamma.wan + mu)*R3
        dS3 <- gamma.wan*R3 - S3 * (1 - C * epsilon) * beta * I / N * 3/4 - mu*S3
        dI2 <- S3 * (1 - C * epsilon) * beta * I / N * 3/4 - (gamma + mu)*I2
        dR2 <- gamma*I2 - (gamma.wan + mu)*R2
        dS2 <- gamma.wan*R2 - S2 * (1 - C * epsilon) * beta * I / N * 2/4 - mu*S2 
        dI1 <- S2 * (1 - C * epsilon) * beta * I / N * 2/4 - (gamma + mu)*I1
        dR1 <- gamma*I1 - (gamma.wan + mu)*R1
        dS1 <- gamma.wan*R1 - S1 * (1 - C * epsilon) * beta * I / N * 1/4 - mu*S1 
        dI0 <- S1 * (1 - C * epsilon) * beta * I / N * 1/4 - (gamma + mu)*I0
        dR0 <- gamma*I0 - mu*R0
        dcum.inc <- (S4 * 4/4 + S3 * 3/4 + S2 * 2/4 + S1 * 1/4) * (1 - C * epsilon) * beta * I / N

        return(list(c(dS4, dI3, dR3, dS3, dI2, dR2, dS2, dI1, dR1, dS1, dI0, dR0, dcum.inc)))
    })
}

## two patch model  with serotypes and variable wMel
SIR.twopatch.births.varC.sero <- function(t, state, parms, Ct.fn, Cc.fn) {
    with(as.list(parms), {

        beta <- beta.varT(t,b.mean=R0*gamma,b.amp=beta.amp,offset=offset)
        C <- c(Ct.fn(t),Cc.fn(t))
        S4 <- state[names(state)=="S4"]
        I3 <- state[names(state)=="I3"]
        R3 <- state[names(state)=="R3"]
        S3 <- state[names(state)=="S3"]
        I2 <- state[names(state)=="I2"]
        R2 <- state[names(state)=="R2"]
        S2 <- state[names(state)=="S2"]
        I1 <- state[names(state)=="I1"]
        R1 <- state[names(state)=="R1"]
        S1 <- state[names(state)=="S1"]
        I0 <- state[names(state)=="I0"]
        R0 <- state[names(state)=="R0"]
        cum.inc <- state[names(state)=="cum.inc"]
        I <- I3 + I2 + I1 + I0
        dS4 <-  mu * (1 - S4) - S4 * (rho.ij %*% ((1-C*epsilon)*beta*I/N)) * 4/4
        dI3 <- S4 * (rho.ij %*% ((1-C*epsilon)*beta*I/N)) * 4/4 - (gamma + mu)*I3
        dR3 <- gamma*I3 - (gamma.wan + mu)*R3
        dS3 <- gamma.wan*R3 - S3 * (rho.ij %*% ((1-C*epsilon)*beta*I/N)) * 3/4 - mu*S3
        dI2 <- S3 * (rho.ij %*% ((1-C*epsilon)*beta*I/N)) * 3/4 - (gamma + mu)*I2
        dR2 <- gamma*I2 - (gamma.wan + mu)*R2
        dS2 <- gamma.wan*R2 - S2 * (rho.ij %*% ((1-C*epsilon)*beta*I/N)) * 2/4 - mu*S2 
        dI1 <- S2 * (rho.ij %*% ((1-C*epsilon)*beta*I/N)) * 2/4 - (gamma + mu)*I1
        dR1 <- gamma*I1 - (gamma.wan + mu)*R1
        dS1 <- gamma.wan*R1 - S1 * (rho.ij %*% ((1-C*epsilon)*beta*I/N)) * 1/4 - mu*S1 
        dI0 <- S1 * (rho.ij %*% ((1-C*epsilon)*beta*I/N)) * 1/4 - (gamma + mu)*I0
        dR0 <- gamma*I0 - mu*R0
        dcum.inc <- (S4 * 4/4 + S3 * 3/4 + S2 * 2/4 + S1 * 1/4) * (rho.ij %*% ((1-C*epsilon)*beta*I/N))

        return(list(c(dS4, dI3, dR3, dS3, dI2, dR2, dS2, dI1, dR1, dS1, dI0, dR0, dcum.inc)))
    })
}
