# loss functions for estimating IAR in one or two patch scenarios
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3506030/
loss.one = function(pi, S0, R0){
  abs(pi - 1 + S0 * exp(-R0 * pi))
}
loss.two = function(pi.t, pi.c, rho.tt, rho.tc, rho.cc, rho.ct, Ct, Cc, epsilon, S0, R0){
  abs(pi.t - 1 + S0 * exp(-(
    pi.t * rho.tt * (1 - Ct * epsilon) * R0 +
    pi.c * rho.tc * (1 - Cc * epsilon) * R0))) ^ 2 +
  abs(pi.c - 1 + S0 * exp(-(
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
