function fl = intl(m,T0,p0,q0,zrange)
%INTL       Integrate the liquid-flow region in the membrane.
%  FL = INTL(M,T0,P0,Q0,ZRANGE) with ZRANGE = [Z0 Z1] integrates the
%  liquid-flow region in the membrane from Z0 to Z1 for given mass flux
%  M and initial conditions T0, P0 and heat flux Q0. Returns a structure
%  FL that can be used with DEVAL to evaluate the solution. INTL does
%  not terminate if the pressure drops below the saturation pressure.
%
%  Calls NUL, K, KAPPA, KELV, CPL, PS, ODE45.
%  Called from ADDINT.
%
%  See also DEVAL.

options=odeset('InitialStep',5e-7); %,'Events',@terml);
fl = ode45(@odeg,zrange,[T0 p0 q0],options,m);

%-----------------------------------------------------------------------
function dy = odeg(z,y,m)
T=y(1);
p=y(2);
q=y(3);
dp=-m*nul(T)/kappa;
dT=-q/k(T,0);
dq=-m*cpl(T)*dT; % + dhdp(T,p)*dp );
dy=[dT;dp;dq];

%-----------------------------------------------------------------------
function [val,isterm,direction] = terml(z,y,m)
isterm=1;
direction=1;
% OK: ps(T)-p < 0
% Termination: ps(T)-p>=0, only when increasing
%  without direction, integration already terminates at start
val=kelv(y(1)).*ps(y(1))-y(2);
