function dpk = dpkdT(T)
%DPKDT(T)   Derivative of the vapor pressure at a curved interface [Pa/K].
%
%  dpk/dT = (pk/ps)*(dps/dT + psMv'/RT*curv(sig/T + dsig/dT)
%
%  Calls CURV, DPSDT, DSIGDT, KELV, MOLM, PS, RHO, SIG.

[R M]=molm;
dpk = kelv(T).*(dpsdT(T)+ps(T)*M*curv*(sig(T)./T-dsigdT(T))./(R*rho(T).*T));
