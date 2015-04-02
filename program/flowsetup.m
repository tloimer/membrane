function flsetup = flowsetup(T2,Tmax,theta,s,mem,f)
%FLOWSETUP  Setup flow properties.
%  FLOWSETUP(T2,TMAX,THETA,S,MEM,F) returns a struct FS that contains flow
%  properties for a combination of contact angle THETA (in degrees), a
%  substance S, a homogeneous membrane MEM and a fmodel F. The integral of
%  dh/dp over temperature T is calculated between T2 and at least TMAX.
%
%  FREEFLOW = FLOWSETUP(S) sets the flowstruct FREEFLOW to properties
%  appropriate for flow through free space. Uses the homogeneous flow
%  model, FMODEL('plug'). Sets FREEFLOW.hgK and FREEFLOW.intdhdpsatdT to
%  empty matrices.
%
%  Fields:
%    FS.curv              Curvature.
%    FS.kelv(T,sigma,rho) pK/psat
%    FS.pkps(T)           pK/psat
%    FS.pkelv(T)          pK
%    FS.pkpcap(T)         pk and capillary pressure pcap, returns [pk pcap]
%    FS.dpkdT(T)          [dpK/dT pk]
%    FS.hgK(T)            Specific enthalpy of the vapor, h(T,pk(T))
%    FS.hvapK(T)          Enthalpy of vaporization [J/kg]
%    FS.hvapKraw(T,...)   Enthalpy of vaporization [J/kg]
%    FS.q2ph(m,T,a)       Heat flux, and two-phase pressure in two-phase flow
%    FS.qminqmax(m,T)     Minimum and maximum heat flux for vapor and liquid
%    FS.intdhdpdpsatdT(T) Int_T2^Tmax dh/dp dpsat/dT dT
%    FS.nuapp(T,p)        Apparant vapor viscosity (viscous + molecular flow)
%    FS.knudsen(T,p)      Knudsen numbe.
%    FS.nu2ph(T,pk,a)     Apparent 2ph viscosity, using app. vapor viscosity
%    FS.kmgas(T)          Thermal conductivity of vapor-filled membran.
%    FS.kmliq(T)          Thermal conductivity of liquid-filled membrane
%    FS.k2ph(T,a)         Thermal conductivity of two-phase filled membrane
%    FS.xdot(T,pk,a)      Vapor mass flow fraction
%    FS.odemaxstep(range,delta)   Minimum number of integration steps
%
%  See also SUBSTANCE, MEMBRANE, MSTACKSTRUCT, FMODEL.

flsetup = struct('curv',[],'kelv',[],'pkps',[],'pkelv',[],'pkpcap',[],...
  'dpkdT',[],'hgK',[],'hvapK',[],'hvapKraw',[],'q2ph',[],'qminqmax',[],...
  'intdhdpdpsatdT',[],'nuapp',[],'knudsen',[],'nu2ph',[],...
  'kmgas',[],'kmliq',[],'k2ph',[],'xdot',[],'odemaxstep',[]);

% set for free space
if nargin == 1
  s = T2;
  theta = 90;
  mem = membrane('free');
  f = fmodel('plug');
  T2 = [];
  isfreespace = true;
  issupercritical = false; % No error checks here for free space.
else
  issupercritical = ~isfinite(s.ps(T2));
  isfreespace = strcmp(mem.tname,'free');
end

% First set what is set in any case.
% Return early further below.

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

if theta == 90 || isfreespace || issupercritical
  flsetup.curv = 0;
  flsetup.kelv = @(T,sigma,rho) 1;
  if issupercritical
    flsetup.pkps = @(T) NaN;
    flsetup.pkelv = @(T) Inf;
    flsetup.pkpcap = @Infzero;
    flsetup.dpkdT = @(T) NaN; % should return 2 arguments;
    flsetup.hgK = @(T) NaN;
    flsetup.intdhdpdpsatdT = @(T) NaN;
    flsetup.hvapK = @(T) 0;
    flsetup.hvapKraw = @(T,prad,psat,pcap,rho,drho) 0;
    % this should be sufficient, to trigger an error if necessary
    return
  elseif isfreespace
    flsetup.pkps = @(T) 1;
    flsetup.pkelv = @(T) s.ps(T);
    flsetup.pkpcap = @psatpcap;
    flsetup.dpkdT = @dpsatdT;
    flsetup.hvapK = @hvapsat;
    flsetup.hvapKraw = @hvapKraw;
    % simply do not set hgK and intdhdp... for free space
    flsetup.qminqmax = @(m,T) qminqmaxfree(T);
    return
  end
  % do not return for theta = 90
else
  flsetup.curv = mem.fcurv(cos(theta*pi/180));
  flsetup.kelv = @(T,sigma,rho) exp(-flsetup.curv.*sigma./(s.R.*rho.*T));
end

flsetup.pkps = @(T) flsetup.kelv(T,s.sigma(T),s.rho(T));
flsetup.pkelv = @(T) s.ps(T)*flsetup.kelv(T,s.sigma(T),s.rho(T));
flsetup.pkpcap = @pkpcap;
flsetup.dpkdT = @dpkdT;
flsetup.hvapK = @hvapK;
flsetup.hvapKraw = @hvapKraw;
flsetup.q2ph = @q2ph; % Also works in free space or for theta = 90.
flsetup.qminqmax = @qminqmax;
% solhhK needs to exist, to be 'deval'uated
% set further below
%flsetup.hgK = @(T) deval(solhgK,T);
%flsetup.intdhdpdpsatdT = @(T) deval(soldhdps,T);

% Calculate hgK and intdhdpsdpsatdT.
% inttol: relative tolerance for integration
inttol = 1e-6;
stepT = 0.4;

% do an integer number of steps
intTrange = ceil((Tmax-T2)/stepT)*stepT;
% the large interval is probably needed for the shooting-iteration
if intTrange > 4
  intTrange = [T2 T2+2*intTrange]; %DEBUG use 2 for normal use
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

%%% NESTED FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NESTED FUNCTIONS %%%

function [pk, pcap] = pkpcap(T) %---------------------------------------- pkpcap
%PKPCAP     Vapor pressure at a curved meniscus and capillary pressure.
% [PK, PCAP] = PKPCAP(T) returns the equilibrium pressure at a curved meniscus,
% PK, and the capillary pressure due to Young's-Laplace equation.
sigma = s.sigma(T);
pk = s.ps(T)*flsetup.kelv(T,sigma,s.rho(T));
pcap = flsetup.curv*sigma;
end %---------------------------------------------------------------- end pkpcap

function [psat, pcap] = psatpcap(T) %---------------------------------- psatpcap
psat = s.ps(T);
pcap = 0;
end %-------------------------------------------------------------- end psatpcap

function [infinity, nix] = Infzero(T) %--------------------------------- Infzero
infinity = Inf;
nix = 0;
end %--------------------------------------------------------------- end Infzero

function [dpk pk] = dpkdT(T) %-------------------------------------------- dpkdT
%DPKDT      Derivative of the equilibrium pressure at a curved meniscus, dpk/dT.
% [DPK PK] = DPKDT(T) returns pk and dpk/dT. A copy from HVAPK, for convenience.

% See below.
[psat dps] = s.ps(T);
[dsig sigma] = s.dsig(T);
[drho rho] = s.drho(T);
pk_ps = flsetup.kelv(T,sigma,rho);
pk = pk_ps*psat;
dpk = pk_ps * (dps + psat*flsetup.curv*sigma*(1/T-dsig+drho)/(s.R*rho*T));
end %----------------------------------------------------------------- end dpkdT

function [dps, psat] = dpsatdT(T) %------------------------------------- dpsatdT
[psat, dps] = s.ps(T);
end %--------------------------------------------------------------- end dpsatdT

function [pk dpk hvapK dpcap pcap] = hvapK(T) %--------------------------- hvapK
%HVAPK      Enthalpy of vaporization at a curved interface.
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

function [psat dps hvapK dpcap pcap] = hvapsat(T) %--------------------- hvapsat
[psat dps] = s.ps(T);
pcap = 0;
dpcap = 0;
hvapK = s.hvap(T);
end %--------------------------------------------------------------- end hvapsat

function hvK = hvapKraw(T,prad,psat,pcap,rho,drho) %------------------- hvapKraw
  dhldp = (1 + T*drho)/rho;
  hvK = s.hvap(T) + (prad-psat)*(s.dhdp(T,prad)-dhldp) + dhldp*pcap;
end %-------------------------------------------------------------- end hvapKraw

function [q,p2ph,pk,pcap] = q2ph(m,T,a) %---------------------------------- q2ph
%Q2PH       Heat flux in two-phase flow.
%  [Q,P2PH,PK,PCAP] = Q2PH(M,T,A) Return heat flux Q [W/m2], two-phase pressure
%  P2PH = PK - (1-A)*PCAP, PK and PCAP, PCAP = pgas - pliq.

%  The bulk of this is copied from hvapK.
%  8<--  % copy from hvapk
[psat dps] = s.ps(T);
[dsig sigma] = s.dsig(T);
[drho rho] = s.drho(T);
pk_ps = flsetup.kelv(T,sigma,rho);
pk = pk_ps*psat;
pcap = flsetup.curv*sigma;
dpcap = pcap*dsig;
dpk = pk_ps * (dps + psat*pcap*(1/T-dsig+drho)/(s.R*rho*T));
%  8<--
q = m*flsetup.nu2ph(T,pk,a)*flsetup.k2ph(T,a)/((dpk-(1-a)*dpcap)*mem.kappa);
p2ph = pk - (1-a)*pcap;
end %------------------------------------------------------------------ end q2ph

function [qmin qmax hvapK dpk dpcap] = qminqmax(m,T) %----------------- qminqmax
%QMINQMAX   Minimum and maximum heat flux for vapor and liquid flow, respectively.
% [QMIN QMAX HVAPK DPK DPCAP] = QMINMAX(T) returns the minimum heat flux, such
% that a vapor remains a vapor flow, and the maximum heat flux allowed for a
% liquid to remain a liquid flow. Both liquid and vapor are at their states in
% equilibrium with the other phase at a curved meniscus.

% The minimum heat flux, such that a vapor remains a vapor (does not
% become too cold upstream of 1) is found from
%   m = -(kappa/nu) dp/dz,
%   q = -k dT/dz.
% A vapor stays marginally a vapor if p follows pk(T), dp/dT = dpk/dT,
%   m = -(kappa/nu) dpk/dT dT/dz,  m = (kappa/nu)*(dpk/dT)*q/k,
%   qmin = m*nu*k/(kappa*dpk/dT).
% Analoguous, a liquid stays marginally a liquid if p follows pk - pcap,
%   qmax = m*nu*k/(kappa*(dpk/dT-dpcap/dT)).
[pk dpk hvapK dpcap pcap] = flsetup.hvapK(T);
qmin = m*flsetup.kmgas(T)*flsetup.nuapp(T,pk)/(mem.kappa*dpk);
qmax = m*flsetup.kmliq(T)*s.nul(T)/(mem.kappa*(dpk-dpcap));
end %-------------------------------------------------------------- end qminqmax

function [qmin qmax hvapK dpk dpcap] = qminqmaxfree(T) %----------- qminqmaxfree
%QMINQMAXFREE QMINQMAX for free space.
% See also FLOWSETUP>QMINQMAX.

[pk dpk hvapK dpcap pcap] = flsetup.hvapK(T);
% Two-phase flow in free space cannot tolerate any heat flux. A vapor must not
% become colder upstream of its saturation state in the free space, hence qmin =
% 0. A liquid must not become warmer, hence qmax = 0.
qmin = 0;
qmax = 0;
end %---------------------------------------------------------- end qminqmaxfree

function nu2app = nu2phapp(T,pk,a) %----------------------------------- nu2phapp
%NU2PHAPP   Apparent viscosity (Knudsen + viscous flow) of the 2ph-mixture.
% NU2PHAPP(T,P,A) returns the apparent kinematic viscosity [m2/s] that results
% from assuming an apparent viscosity for the vapor flow,
%   nuapp = nug / (1 + beta*Kn).
vol = s.v(T,pk);
nu2app = f.nu2ph(a,vol,1/s.rho(T),flsetup.nuapp(T,pk)/vol,s.mul(T));
end %-------------------------------------------------------------- end nu2phapp

end %%% END FLOWSETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END FLOWSETUP %%%
