function flsetup = flowsetup(T2,Tmax,theta,s,mem,f) %----------------- flowsetup
%FLOWSETUP  Setup of flow properties.
% FS = FLOWSETUP(T2,Tmax,THETA,S,MEM,F) returns a struct FS that contains flow
% properties.
%  Fields:
%    FS.curv              Curvature.
%    FS.kelv(T,sigma,rho) pK/psat
%    FS.pkelv(T)          pK
%    FS.dpkdT(T)          [dpK/dT pk]
%    FS.hgK(T)            Specific enthalpy of the vapor, h(T,pk(T)).
%    FS.hvapK(T)          Enthalpy of vaporization [J/kg].
%    FS.hvapKraw(T,...)   Enthalpy of vaporization [J/kg].
%    FS.intdhdpdpsatdT(T) Int_T2^Tmax dh/dp dpsat/dT dT.
%    FS.nuapp(T,p)        Apparant vapor viscosity (viscous + molecular flow).
%    FS.knudsen(T,p)      Knudsen number.
%    FS.nu2ph(T,pk,a)     Apparent 2ph viscosity, using app. vapor viscosity.
%    FS.kmgas(T)          Heat conductivity of vapor-filled membrane.
%    FS.kmliq(T)          Heat conductivity of liquid-filled membrane.
%    FS.k2ph(T,a)         Heat conductivity of two-phase filled membrane.
%    FS.xdot(T,pk,a)      Vapor mass flow fraction.
%    FS.odemaxstep(range,delta) Minimum number of integration steps.

flsetup = struct('curv',[],'kelv',[],'pkelv',[],'dpkdT',[],'hgK',[],...
  'hvapK',[],'hvapKraw',[],'intdhdpdpsatdT',[],'nuapp',[],'knudsen',[],...
  'nu2ph',[],'kmgas',[],'kmliq',[],'k2ph',[],'xdot',[],'odemaxstep',[]);

flsetup.curv = mem.fcurv(cos(theta*pi/180));
flsetup.kelv = @(T,sigma,rho) exp(-flsetup.curv.*sigma./(s.R.*rho.*T));

% Apparent vapor viscosity
%   nuapp = nug/(1 + beta Kn),  Kn = 3 nug sqrt(pi/(8*R*T)) / dia,
if mem.beta
  kn_nu = @(T) 3*sqrt(pi/(8*s.R)) / (sqrt(T)*mem.dia);
  flsetup.knudsen = @(T,p) kn_nu(T)*s.nug(T,p);
  flsetup.nuapp = @(T,p) 1 / (1/s.nug(T,p) + mem.beta*kn_nu(T));
  flsetup.nu2ph = @nu2phapp;
else % mem.beta = 0
  % Viscous flow
  flsetup.knudsen = @(T,p) 0;
  flsetup.nuapp = @(T,p) s.nug(T,p);
  flsetup.nu2ph = @(T,pk,a) f.nu2ph(a,s.v(T,pk),1/s.rho(T),s.mug(T),s.mul(T));
end

% Flow model and effective two-phase properties.
flsetup.kmgas = @(T) f.kmgas(mem.epsilon,mem.km,s.kg(T));
flsetup.kmliq = @(T) f.kmliq(mem.epsilon,mem.km,s.kl(T));
flsetup.k2ph = @(T,a) f.k2ph(a,mem.epsilon,mem.km,s.kg(T),s.kl(T));
flsetup.xdot = @(T,pk,a) f.xdot(a,s.v(T,pk),1/s.rho(T));
%flsetup.x = @(T,a) f.x(a,s.v(T,pk(T)),1/s.rho(T));

% Step size calculation for ode45, ode23t.
% minsteps (= 2): minimum number of integration steps;
% odemaxstep: number of (power of 2: 2, 4, 8, ...) steps such that T or p does
% not change more than maxTperstep (maxpperstep) within one step.
flsetup.odemaxstep = ... % 1/max(minsteps,...
  @(range,maxstep) 1/max(2,pow2(ceil(log2(abs(range)/maxstep))));

% Calculate hgK, intdhdpsdpsatdT and pkelv. Check if T > Tc, i.e., ps(T) == Inf.
% Return Inf or NaN if the fluid is anywhere supercritical.
% Might be changed with moderate effort to return functions that give exact
% values.
if ~isfinite(s.ps(T2))
  flsetup.pkelv = @(T) Inf;
  flsetup.hgK = @(T) NaN;
  flsetup.intdhdpdpsatdT = @(T) NaN;
  flsetup.hvapK = @(T) 0;
  flsetup.hvapKraw = @(T,prad,psat,pcap,rho,drho) 0;
  return
end

flsetup.dpkdT = @dpkdT;
flsetup.pkelv = @(T) s.ps(T)*flsetup.kelv(T,s.sigma(T),s.rho(T));
%flsetup.hgK = @(T) deval(solhgK,T);
%flsetup.intdhdpdpsatdT = @(T) deval(soldhdps,T);
flsetup.hvapK = @hvapK;
flsetup.hvapKraw = @hvapKraw;

% inttol: relative tolerance for integration
inttol = 1e-6;
stepT = 0.4;

% do an integer number of steps
%intTrange = ceil(intTrange/stepT)*stepT;
intTrange = ceil((Tmax-T2)/stepT)*stepT;
if intTrange > 4
  intTrange = [T2 T2+4*intTrange]; %DEBUG use 2 for normal use
else
  intTrange = [T2 T2+10];
end
%intTrange = [295 305];

options = odeset('RelTol',inttol,'Refine',1,...
  'InitialStep',stepT,'MaxStep',stepT);
%solhgK = ode45(@inthgK,[Trange],initial condition,options);
solhgK = ode45(@inthgK,intTrange,0,options);
function dy = inthgK(T,y)
  % Calculate the enthalpy of the vapor along (T, pk(T)). h = 0 at T = T2.
  % With dh = (dh/dp) dp + (dh/dT) dT, dp = (dpk/dT)dT, (dh/dT) = cpg, follows
  %   dh/dT = (dh/dp) dpK/dT + cpg
  % One equation: Scaling (making dimensionless) is not necessary.
  [dpk pk] = dpkdT(T);
  [dhdp cpg] = s.dhcpg(T,pk);
  dy = dpk*dhdp + cpg;
end
% Quadrature, with fun(T) = dpk*dhdp + cpg and delhgK = quad(fun,T1,T2) would
% also have been possible. But quadrature would have to be evaluated each time.
% Instead, solhgK returns h(T).

% options stay the same as before
soldhdps = ode45(@intdhdpdpsdT,intTrange,0,options);
function dh = intdhdpdpsdT(T,h)
  % calculate the term
  %   T2_int^T dh/dp dpsat/dT' dT'
  % in the liquid phase, hence dh/dp = v - Tdv/dT = 1/rho + (T/rho^2) drho/dT.
  [drho rho] = s.drho(T); [ps dps] = s.ps(T);
  dh = dps*(1+T*drho)/rho;
end

flsetup.hgK = @(T) deval(solhgK,T);
flsetup.intdhdpdpsatdT = @(T) deval(soldhdps,T);

function [dpk pk] = dpkdT(T) %-------------------------------------------- dpkdT
%DPKDT      Derivative of the equilibrium pressure at a curved meniscus, dpk/dT.
%
% [DPK PK] = DPKDT(T) returns the pk and dpk/dT. A copy from HVAPK, for
% convenience.
%
% See below.
[psat dps] = s.ps(T);
[dsig sigma] = s.dsig(T);
[drho rho] = s.drho(T);
pk_ps = flsetup.kelv(T,sigma,rho);
pk = pk_ps*psat;
dpk = pk_ps * (dps + psat*flsetup.curv*sigma*(1/T-dsig+drho)/(s.R*rho*T));
end %----------------------------------------------------------------- end dpkdT

function [pk dpk hvapK dpcap pcap] = hvapK(T) %--------------------------- hvapK
%HVAPK      Enthalpy of vaporization at a curved interface.
%
% [PK DPK HVK DPCAP PCAP] = HVAPK(T) returns the enthalpy of vaporization,
% the vapor pressure, the derivatives of the vapor pressure and the capillary
% pressure with respect to temperature at a curved interface as well as the
% capillary pressure. [PK] = Pa, [DPK] = Pa/K, [HVK] = J/kg [DPCAP] = Pa/K,
% [PCAP] = Pa.
%
% DPK: See Eq. (11) in Loimer (2007).
%
% HVK: See Eq. (14) in [Loimer, Proc. STAMM 2004; Wang, Hutter (eds.), Trends in
% Applications of Mathematics to Mechanics, 247-256, Shaker (2005)].
% With dh/dp = v - T dv/dT.
[psat dps] = s.ps(T);
[dsig sigma] = s.dsig(T);
[drho rho] = s.drho(T);
pk_ps = flsetup.kelv(T,sigma,rho);
pk = pk_ps*psat;
pcap = flsetup.curv*sigma;
dpcap = pcap*dsig;
dpk = pk_ps * (dps + psat*pcap*(1/T-dsig+drho)/(s.R*rho*T));
%hvapK = s.hvap(T) + (pk-psat)*(s.dhdp(T,pk)-dhldp) + dhldp*curv*sigma;
hvapK = flsetup.hvapKraw(T,pk,psat,pcap,rho,drho);
end %----------------------------------------------------------------- end hvapK

function hvK = hvapKraw(T,prad,psat,pcap,rho,drho) %------------------- hvapKraw
  dhldp = (1 + T*drho)/rho;
  hvK = s.hvap(T) + (prad-psat)*(s.dhdp(T,prad)-dhldp) + dhldp*pcap;
end %-------------------------------------------------------------- end hvapKraw

function nu2app = nu2phapp(T,pk,a) %----------------------------------- nu2phapp
%NU2PHAPP   Apparent viscosity (Knudsen + viscous flow) of the 2ph-mixture.
%
% NU2PHAPP(T,P,A) returns the apparent kinematic viscosity [m2/s] that results
% from assuming an apparent viscosity for the vapor flow,
%   nuapp = nug / (1 + beta*Kn).
vol = s.v(T,pk);
nu2app = f.nu2ph(a,vol,1/s.rho(T),flsetup.nuapp(T,pk)/vol,s.mul(T));
end %-------------------------------------------------------------- end nu2phapp
end %------------------------------------------------------------- end flowsetup