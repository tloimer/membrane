function rkelv = rkelv(T)
%RKELV(T)   Evaporation enthalpy at a curved interface [J/kg].
%  Correction to evaporation enthalpy:
%  rkelv = r + (pk -ps)*(dh''/dp - dh'/dp) + (dh'/dp)*pcap
%
%  Calls CURV, DHDP, DPSDT, DRHODT, KELV, MOLM, PS, R, RHO, SIG.

[R M]=molm;
rh = rho(T);
dhldp = (1+T.*drhodT(T)./rh)./rh;
rkelv = r(T) + (kelv(T)-1).*ps(T).*(dhdp(T)-dhldp) + dhldp*curv.*sig(T);
%rs = r(T) + curv*sig(T)./rh
