function f2ph = int2ph(m,T0,a0,q0,zrange)
%INT2PH     Integrate the 2ph-flow region in the membrane.
%  FL = INT2PH(M,T0,A0,Q0,ZRANGE) with ZRANGE = [Z0 Z1] integrates the
%  2ph-region in the membrane from Z0 to Z1 for given mass flux M and
%  initial conditions T0, vapor volume fraction A0 and heat flux Q0.
%  The initial jump in vapor content is calculated. Vapor content is
%  checked to be within the range [0 1].
%
%  Calls NU, K, XDOT, KAPPA, Q_M, INTCPL, R, DPSDT, ODE23T.
%  Called from FLOWBACK.
%
%  See also DEVAL.

% initial enthalpie (flux - canceled)
% for gas-2ph this is: q0=0,a0=1,xdot(T0,1)=1
h0 = q0/m + xdot(T0,a0)*r(T0);

% calculate initial vapor content a3
a3=fzero(@ares,[0 1],optimset('fzero'),T0,T0,h0);
  
% and initial temperature gradient
dTdz1=-m*q_m(T0,a3)/k(T0,a3);

% ode15s  stiff, low to medium accuracy
% ode23t  moderately stiff, low accuracy, solution without numerical damping
options=odeset('Events',@term2ph,'Mass',[1 0;0 0],'MStateDependence','none',...
  'MassSingular','yes','InitialSlope',[dTdz1;0],'InitialStep',5e-7);
f2ph = ode23t(@ode2ph,zrange,[T0 a3],options,m,T0,h0);

%-----------------------------------------------------------------------
function ar = ares(a,T,T0,h0)

% ar= h'dT+xdot*r(T)*q/m-h0;
ar=intcpl(T0,T) + xdot(T,a)*r(T) + q_m(T,a) - h0;

%-----------------------------------------------------------------------
function dy = ode2ph(z,y,m,T0,h0)
T=y(1);
a=y(2);
dT=-m*nu(T,a)/(dpsdT(T)*kappa);
da=ares(a,T,T0,h0);
dy=[dT;da];

%-----------------------------------------------------------------------
function [val,isterm,direction] = term2ph(z,y,m,T0,h0)
%TERM2PH    Check vapor content to be between [1 0]. A direction of
%  zero-crossing must be given, otherwise integration stops at start.
a=y(2);
isterm=[1;1];
direction=-ones(2,1);
val=[1-a;a];
