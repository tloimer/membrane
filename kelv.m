function kelv=kelv(T)
%KELV(T)    Vapor pressure reduction, p/ps.
%
%  Returns the vapor pressure reduction, KELV = P/PS, where P is the
%  reduced vapor pressure.
%
%  Calls MOLM, PCAP, RHO.

% to disable this:
%kelv=1;

[R M]=molm;
kelv=exp(-curv*sig(T)*M./(R*rho(T).*T));
