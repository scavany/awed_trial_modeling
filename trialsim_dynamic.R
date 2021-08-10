model = function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    rho_tc = 1 - rho_tt
    rho_ct = 1 - rho_cc
    
    # forces of infection
    FOIt = beta * (1-eff*Ct) * rho_tt * It + beta * (1-eff*Cc) * rho_tc * Ic
    FOIc = beta * (1-eff*Ct) * rho_ct * It + beta * (1-eff*Cc) * rho_cc * Ic
    
    # rate of change
    dSt = - FOIt
    dIt = FOIt * St - gamma * It
    dSc = - FOIc
    dIc = FOIc * Sc - gamma * Ic
    
    # return the rate of change
    list(c(dSt, dIt, dSc, dIc))
  }) # end with(as.list ...
}



parameters = c(
  beta = R0/14,
  gamma = 1/14,
  rho_tt = 1,
  rho_cc = 1,
  Ct = 0.958,
  Cc = 0.15,
  eff = epsilon.trial[2],
  N = 1)

state = c(
  St = S0,
  It = 1e-5,
  Sc = S0,
  Ic = 1e-5)

times = seq(0, 365, by = 1)

out = ode(y = state, times = times, func = model, parms = parameters)
plot(out)
