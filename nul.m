function nul=nul(T)
%NUL(T)     Kinematic viscosity of the liquid [m2/s].
%
%  Calls MUL, RHO.

nul = mul(T)./rho(T);
