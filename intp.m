function fp = intp(m,T0,p0,q0,zrange,cp,k)
%INTP       Integrate Darcy's law in the liquid-flow region in the membrane.
%  FL = INTP(M,T0,P0,Q0,ZRANGE,CP,K) with ZRANGE = [Z0 Z1] integrates
%  darcy's law in the liquid-flow region in the membrane from Z0 to Z1
%  for given mass flux M, initial conditions T0, P0 and heat flux Q0.
%  The temperature distribution is obtained by calling TEMP with
%  material properties CP and K. INTP returns the pressure distribution
%  in a structure FP that can be used with DEVAL to evaluate the
%  solution. INTP terminates if the pressure plus a pressure difference
%  due to capillary pressure drops below the saturation pressure.
%% INTP does not terminate if the pressure drops below the
%% saturation pressure but detects such an event.
%
%  Calls KAPPA, K, NUL, ODE45, PCAP, PS, TEMP.
%  Called from FLOWBACK.
%
%  See also DEVAL.

options=odeset('InitialStep',5e-7,'Events',@termd);
fp = ode45(@darcy,zrange,p0,options,m,T0,q0,zrange(1),cp,k);

%-----------------------------------------------------------------------
function dp = darcy(z,p,m,T0,q0,z0,cp,k)
dp=-m*nul(temp(z,m,T0,q0,z0,cp,k))/kappa;

%-----------------------------------------------------------------------
function [val,isterm,direction] = termd(z,p,m,T0,q0,z0,cp,k)
isterm=1; %terminate
%isterm=0; %do not terminate
direction=1;
% OK: ps(T)-p < 0
% Termination: ps(T)-p>=0, only when increasing
%  without direction, integration already terminates at start
%val=ps(temp(z,m,T0,q0,z0,cp,k))-p;
T = temp(z,m,T0,q0,z0,cp,k);
val = ps(T)-p-pcap(T); % also pcap subtracted
