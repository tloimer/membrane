function fg = intg(m,T0,p0,q0,zrange)
%INTG       Integrate the gas-flow region in the membrane.
%  FL = INTG(M,T0,P0,Q0,ZRANGE) with ZRANGE = [Z0 Z1] integrates the
%  gas-flow region in the membrane from Z0 to Z1 for given mass flux M
%  and initial conditions T0, P0 and heat flux Q0. Returns a structure
%  FL that can be used with DEVAL to evaluate the solution. INTG
%  terminates it the pressure rises above the saturation pressure.
%
%  Calls NUG, KAPPA, K, CPG, DHDP, PS, ODE45.
%  Called from FLOWBACK.
%
%  See also DEVAL.

options=odeset('InitialStep',5e-7,'Events',@termg);
fg = ode45(@odeg,zrange,[T0 p0 q0],options,m);

%-----------------------------------------------------------------------
function dy = odeg(z,y,m)
T=y(1);
p=y(2);
q=y(3);
dp=-m*nug(T,p)/kappa;
dT=-q/k(T,1);
dq=-m*( cpg(T,p)*dT + dhdp(T)*dp );
dy=[dT;dp;dq];

%-----------------------------------------------------------------------
function [val,isterm,direction] = termg(z,y,m)
isterm=1;
direction=-1;
% OK: ps(T)-p > 0
% Termination: ps(T)-p<=0, only when falling
%  without direction, integration already terminates at start
val=ps(y(1))-y(2);
