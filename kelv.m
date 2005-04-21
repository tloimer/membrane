function kelv=kelv(T)
%KELV(T)    Vapor pressure reduction, KELV = pk/ps.
%
%  Returns the vapor pressure reduction, KELV = PK/PS, where PK is the
%  reduced vapor pressure.
%
%  Calls MOLM, CURV, SIG, RHO.

% to disable this:
%kelv=1;

[R M]=molm;
kelv=exp(-curv*sig(T)*M./(R*rho(T).*T));
