function q = qevap(m,T,x)
%QEVAP(m,T,x) Heat flux for complete evaporation.

% xr= h'dT+lam*r(T)*q/m-h0; h0=q0/m
q = m*((lambda(T,x)-1).*r(T) + q_m(T,x));
