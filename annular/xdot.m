function xdot=xdot(T,a)
%XDOT(T,A)  Mobility, annular.
%  Mobility = Vapor mass flux over total mass flux.
%  XDOT = 1- nu*(1-a)^2/nul.
%
%  Calls NU, MUL, RHO.

xdot = 1 - nu(T,a).*(1-a).^2.*rho(T)./mul(T);
