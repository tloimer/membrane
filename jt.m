function jt = jt(T,p,deltap)
%JT         Joule-Thomson coefficient [K/Pa].
%  JT(T,P) returns the differential Joule-Thomson coefficient,
%  -(dh/dp)/cp.
%
%  JT(T,P,DELTAP) returns temperature difference for initial conditions
%  T, P and pressure difference DELTAP.
%
%  Calls DHDP, CPG, ODE45.

if (nargin==2)
  jt = jtloc(p,T);
elseif (nargin==3)
  [pi ti] = ode45(@jtloc,[p p-deltap],T);
  jt = T-ti(end);
end

function jtloc = jtloc(p,T)
  jtloc = -dhdp(T)./cpg(T,p);
